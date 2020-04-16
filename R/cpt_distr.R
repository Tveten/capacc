
#' @export
sim_cpt_distr <- function(out_file,
                          data = init_data(n = 200, p = 10, rho = 0.8, band = 2,
                                           locations = 80, durations = 120),
                          params = method_params(cost = "mvlrt"),
                          n_sim = 100, seed = NA) {
  add_setup_info <- function(res) {
    res <- cbind(res,
                 as.data.table(data[!(grepl("Sigma", names(data)) |
                                        names(data) == "changing_vars")]),
                 as.data.table(params[names(data) != "maxsl"]))
    res$vartheta <- res$mu * sqrt(round(data$proportions * data$p))
    res$n_sim <- n_sim
    res$seed <-  seed
    res
  }

  message(paste0("Estimating changepoint distribution for n=", data$n,
                 ", p=", data$p,
                 ", cost=", params$cost,
                 ", precision=", data$precision_type,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", prop=", data$proportions,
                 ", mu=", round(data$mu, 1),
                 ", location=", data$locations,
                 ", precision_est_struct=", params$precision_est_struct,
                 ", est_band=", params$est_band,
                 "."))

  if (!is.na(seed)) set.seed(seed)
  res <- data.table(cpt = replicate(n_sim, simulate_mvcpt(data, params)$cpt))
  fwrite(add_setup_info(res), paste0("./results/", out_file), append = TRUE)
}

# list("locations" = c(40, 50, 60), "vartheta" = c(0.5, 1, 2))
# Map(sim_cpt_distr,
#     data     = data_list,
#     params   = params_list,
#     MoreArgs = list(out_file = out_file,
#                     seed     = seed,
#                     n_sim    = n_sim))

#' @export
many_cpt_distr <- function(out_file,
                           data = init_data(n = 200, p = 10, rho = 0.8, band = 2,
                                            locations = 80, durations = 120),
                           params = method_params(), n_sim = 100,
                           costs = c("inspect", "mvlrt"), bands = 2,
                           rhos = 0.9, props = 0.1, varthetas = 1,
                           locations = round(data$n / 2),
                           precision_est_structs = "correct", est_bands = NA) {
  lapply(locations, function(location) {
    lapply(varthetas, function(vartheta) {
      lapply(props, function(prop) {
        lapply(rhos, function(rho) {
          lapply(bands, function(band) {
            lapply(precision_est_structs, function(precision_est_struct) {
              lapply(est_bands, function(est_band) {
                seed <- sample(1:10^6, 1)
                lapply(costs, function(cost) {
                  data$locations <- location
                  data$durations <- data$n - location
                  data$mu <- vartheta / sqrt(round(prop * data$p))
                  data$proportions <- prop
                  data$rho <- rho
                  data$band <- band
                  params$cost <- cost
                  params$precision_est_struct <- precision_est_struct
                  if (precision_est_struct == "correct") params$est_band <- NA
                  else params$est_band <- est_band
                  sim_cpt_distr(out_file, data, params, n_sim, seed)
                })
              })
            })
          })
        })
      })
    })
  })
  NULL
}

#' @export
cpt_runs <- function(n_sim = 100) {
  out_file <- "cpt_distr.csv"
  many_cpt_distr(out_file,
                 init_data(n = 100, p = 10),
                 method_params(b = 1),
                 n_sim = n_sim,
                 props = c(0.1, 1),
                 varthetas = c(0.5, 1, 2),
                 rhos = c(0.7, 0.9, 0.99),
                 locations = c(10, 50),
                 precision_est_structs = c("correct", "banded"),
                 est_bands = 0)
  many_cpt_distr(out_file,
                 init_data(n = 100, p = 10),
                 method_params(b = 0.3),
                 n_sim = n_sim,
                 props = c(0.1, 1),
                 varthetas = c(0.5, 1, 2),
                 rhos = c(0.7, 0.9, 0.99),
                 locations = c(10, 50),
                 precision_est_structs = c("correct", "banded"),
                 est_bands = 0)
  many_cpt_distr(out_file,
                 init_data(n = 200, p = 100),
                 method_params(b = 1),
                 n_sim = n_sim,
                 props = c(0.01, 1),
                 varthetas = c(0.5, 1, 2),
                 rhos = c(0.7, 0.9, 0.99),
                 locations = c(10, 50, 100),
                 precision_est_structs = c("correct", "banded"),
                 est_bands = 0)
}

read_single_cpt_distr <- function(file_name, data = init_data(),
                                  params = method_params(), .n_sim = 100) {
  res <- fread(paste0("./results/", file_name))
  res <- res[n == data$n &
               p == data$p &
               rho == data$rho &
               precision_type == data$precision_type &
               block_size == data$block_size &
               is_equal(proportions, data$proportions) &
               locations == data$locations &
               durations == data$durations &
               change_type == data$change_type &
               minsl == params$minsl &
               maxsl == params$maxsl &
               n_sim == .n_sim]

  if (data$precision_type %in% c("banded", "block_banded"))
    res <- res[band == data$band]

  if (params$cost == "cor") {
    res <- res[precision_est_struct == params$precision_est_struct]
    if (is.na(params$est_band)) res <- res[is.na(est_band)]
    else res <- res[est_band == params$est_band]
  }
  res
}

plot_cpt_distr <- function(file_name,
                           data = init_data(n = 100, p = 10, rho = 0.9, band = 2,
                                            locations = 50, durations = 50),
                           params = method_params(),
                           n_sim = 100) {
  get_title <- function() {
      if (data$precision_type == "lattice")
        title <- "lattice"
      if (data$precision_type == "banded")
        title <- paste0(data$band, "-banded")
      if (data$precision_type == "block_banded")
        title <- paste0(data$band, "-banded, m=", data$block_size)

      title <- paste0(title,
                      " rho=", data$rho,
                      ", p=", data$p,
                      ", n=", data$n,
                      ", cpt=", data$locations,
                      ", prop=", round(data$proportions, 2))
  }

  res <- read_single_cpt_distr(file_name, data, params, n_sim)
  print(res)
  # sub_res <- res
  ggplot2::ggplot(res, ggplot2::aes(cpt, colour = cost)) +
    ggplot2::stat_bin(alpha = 0.8, ggplot2::aes(y=..density..), geom = "step",
                      position = "identity", binwidth = 1) +
    # ggplot2::geom_histogram(alpha = 0.5, ggplot2::aes(y = ..density..), position = 'identity', binwidth = 1) +
    ggplot2::ggtitle(get_title()) +
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
  plots <- Map(plot_cpt_distr, prop = mean_prop_combs[, 2], ss = mean_prop_combs[, 1],
               MoreArgs = list("res" = res))
  ggpubr::ggarrange(plotlist = plots, nrow = length(props), ncol = length(varthetas),
                    common.legend = TRUE, legend = "right")
}
