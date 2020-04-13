
simulate_single_cpt <- function(setup = init_data_setup(n = 200, p = 5, rho = 0.8, band = 2, locations = 80, durations = 120),
                                method = "lr", seed = NULL, nb_struct = "true",
                                est_band = 2, b = 1) {
  if (!is.null(seed)) set.seed(seed)
  x <- simulate_cor(n                 = setup$n,
                    p                 = setup$p,
                    mu                = setup$mu,
                    Sigma             = setup$Sigma,
                    locations         = setup$locations,
                    durations         = setup$durations,
                    proportions       = setup$proportions,
                    change_type       = setup$change_type)

  diff_x <- x[2:nrow(x), ] - x[1:(nrow(x) - 1), ]
  if (nb_struct == "true")
    Q_hat <- 2 * estimate_precision_mat(diff_x, setup$Sigma_inv)
  else if (nb_struct == "banded")
    Q_hat <- 2 * estimate_precision_mat(diff_x, adjacency_mat(banded_neighbours(est_band, setup$p)))

  if (method == "inspect") cpt <- single_cor_inspect(t(x), Q_hat)$changepoint
  else if (method == "mvlrt") cpt <- single_mvnormal_changepoint(x, Q_hat, b = b)$cpt
  cpt
}

sim_cpt_distr <- function(method, change_setup, nb_struct, est_band, root_seeds, b) {
  add_setup_info <- function(res_dt) {
    res_dt$method <- rep(method, nrow(res_dt))
    res_dt$rho <- rep(change_setup$rho, nrow(res_dt))
    res_dt$band <- rep(change_setup$band, nrow(res_dt))
    res_dt$cor_mat_type <- rep(change_setup$cor_mat_type, nrow(res_dt))
    res_dt$n <- rep(change_setup$n, nrow(res_dt))
    res_dt$p <- rep(change_setup$p, nrow(res_dt))
    res_dt$proportion <- rep(change_setup$proportions, nrow(res_dt))
    res_dt$change_type <- rep(change_setup$change_type, nrow(res_dt))
    res_dt$mu <- rep(change_setup$mu[1], nrow(res_dt))
    res_dt$vartheta <- rep(res_dt$mu[1] * sqrt(round(change_setup$proportions * change_setup$p)), nrow(res_dt))
    res_dt$location <- rep(change_setup$locations[1], nrow(res_dt))
    res_dt$est_nb_struct <- rep(nb_struct, nrow(res_dt))
    if (method == "inspect") res_dt$b <- rep(NA, nrow(res_dt))
    else if (method == "mvlrt") res_dt$b <- rep(b, nrow(res_dt))
    if (nb_struct == "true") res_dt$est_band <- rep(NA, nrow(res_dt))
    else if (nb_struct == "banded") res_dt$est_band <- rep(est_band, nrow(res_dt))
    res_dt
  }

  cpts <- do.call("c", lapply(root_seeds, function(seed) {
    simulate_single_cpt(change_setup, method, seed, nb_struct, est_band, b)
  }))
  res <- data.table("cpt" = cpts)
  add_setup_info(res)
}

run_single_cpt_sim <- function(out_file, n = 500, p = 20, rhos = c(0.5, 0.7, 0.9, 0.99), band = 2,
                           proportions = c(1/p, round(p^{-5/8}, 2), 1), varthetas = 1,
                           locations = 0.4 * n, durations = n - locations,
                           cor_mat_type = "banded", change_type = "adjacent",
                           nb_struct = "true", est_band = 2, b = 1,
                           n_sim = 100, tol = 10, init_seed = 16) {
  methods <- c('inspect', 'mvlrt')
  set.seed(init_seed)
  start_time <- proc.time()[3]
  res <- do.call("rbind", lapply(varthetas, function(vartheta) {
    do.call("rbind", lapply(proportions, function(prop) {
      do.call("rbind", lapply(rhos, function(rho) {
        start_seed <- sample(1:10^6, 1)
        root_seeds <- (start_seed + 1):(start_seed + n_sim)
        do.call("rbind", lapply(methods, function(method) {
          print(paste0("Signal strength = ", vartheta, ", prop = ", prop,  ', rho = ', rho, ', method = ', method))
          print(paste0("Time elapsed: ", round((proc.time()[3] - start_time) / 60, 1), " min."))
          change_setup <- init_data_setup(n = n,
                                          p = p,
                                          proportions = prop,
                                          mu = vartheta / sqrt(round(prop * p)),
                                          locations = locations,
                                          durations = durations,
                                          rho = rho,
                                          band = band,
                                          change_type = change_type,
                                          cor_mat_type = cor_mat_type)
          sim_cpt_distr(method, change_setup, nb_struct, est_band, root_seeds, b)
        }))
      }))
    }))
  }))
  fwrite(res, file = out_file, append = TRUE)
  res
}

single_cpt_runs <- function() {
  res <- run_mvcapa_sim("single_cpt_test2.csv", n = 200, p = 20, rhos = c(0.5, 0.7, 0.9, 0.99), proportions = c(1/20, 4/20, 1), varthetas = seq(0.25, 1.5, 0.25), n_sim = 200, init_seed = 20, b = 1)
  res <- run_mvcapa_sim("single_cpt_n400p100.csv", n = 400, p = 100, rhos = c(0.5, 0.7, 0.9, 0.99), proportions = c(1/100, 5/100, 1), varthetas = seq(0.25, 1.5, 0.25), n_sim = 200, init_seed = 40, b = 1)
}

cpt_distr <- function(res, prop = 0.2, ss = 1) {
  sub_res <- res[vartheta == ss][proportion == prop]
  # sub_res <- res
  title <- paste0("pr=", round(prop, 2), ", ss=", ss)
  ggplot2::ggplot(sub_res, ggplot2::aes(cpt, colour = method)) +
    ggplot2::stat_bin(alpha = 0.8, ggplot2::aes(y=..density..), geom = "step", position = "identity", binwidth = 1) +
    # ggplot2::geom_histogram(alpha = 0.5, ggplot2::aes(y = ..density..), position = 'identity', binwidth = 1) +
    ggplot2::ggtitle(title) +
    # ggplot2::scale_fill_manual() +
    ggplot2::scale_x_continuous(name = "Estimated cpt",
                                limits = c(0, res$n[1]),
                                breaks = res$location[1] - 1)
}

#' @export
prop_mean_distr <- function(file_name, rho_ = 0.9, b_ = 1, varthetas = NULL, props = NULL,
                            n_obs = 200, n_var = 20, cmt = "banded",
                            band = 2, duration = NA, location = NA,
                            ct = "adjacent", nb_struct = "banded") {
  res <- fread(file_name)
  res <- res[rho == rho_][b == b_ | is.na(b)]
  if (cmt == "lattice") res <- adjust_res_lattice(res)

  if (is.null(varthetas)) varthetas <- unique(res$vartheta)
  if (is.null(props)) props <- unique(res$proportion)

  mean_prop_combs <- expand.grid(varthetas, props)
  plots <- Map(cpt_distr, prop = mean_prop_combs[, 2], ss = mean_prop_combs[, 1],
               MoreArgs = list("res" = res))
  ggpubr::ggarrange(plotlist = plots, nrow = length(props), ncol = length(varthetas),
                    common.legend = TRUE, legend = "right")
}
