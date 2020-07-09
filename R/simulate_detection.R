
#' @export
init_data <- function(n = 100, p = 10, proportions = round(sqrt(p))/p,
                      vartheta = 1, shape = 0, change_seed = NA,
                      locations = 50, durations = 10,
                      change_type = 'adjacent', changing_vars = NA,
                      point_locations = NA, point_proportions = NA,
                      point_mu = NA, precision_type = 'banded', rho = 0.9,
                      band = 2, block_size = p, min_nbs = 1, max_nbs = 3) {
  get_Sigma <- function(precision_type) {
    if (precision_type == 'iid') {
      return(list('mat' = diag(1, p), 'inverse' = diag(1, p)))
    } else if (precision_type == 'ar1') {
      return(list('mat' = ar_cor_mat(p, rho), 'inverse' = ar_precision_mat(p, rho)))
    } else if (precision_type == "global_const") {
      Sigma <- constant_cor_mat(p, rho)
      Q <- Matrix::Matrix(solve(Sigma), sparse = TRUE)
      return(list("mat" = Sigma, "inverse" = Q))
    } else if (precision_type == "Wishart") {
      precision_mat <- Wishart_precision(p = p,
                                         n = 20 * p,
                                         rho = rho,
                                         precision_type = "banded",
                                         band = band,
                                         min_nbs = min_nbs,
                                         max_nbs = max_nbs)
      return(list("mat" = solve(precision_mat), "inverse" = precision_mat))
    } else {
      precision_mat <- block_precision_mat(p                 = p,
                                           m                 = block_size,
                                           rho               = rho,
                                           within_block_type = precision_type,
                                           band              = band,
                                           min_nbs           = min_nbs,
                                           max_nbs           = max_nbs)
      return(list('mat' = solve(precision_mat), 'inverse' = precision_mat))
    }
  }

  if (any(change_type == "custom")) {
    if (!is.na(changing_vars))
      proportions <- length(changing_vars) / p
    else
      stop("If change_type is 'custom', you must provide a changing_vars vector")
  }
  # if (precision_type == "lattice" && change_type == "adjacent")
  #   change_type <- "adjacent_lattice"

  if (isTRUE(all.equal(1/p, proportions))) shape <- 0
  proportions <- ceiling(proportions * p) / p

  block_size <- min(p, block_size)
  Sigma_obj <- get_Sigma(precision_type)
  if (precision_type == "lattice") band <- band(Sigma_obj$inverse)
  else if (precision_type == "global_const") band <- p - 1

  list('n'                 = n,
       'p'                 = p,
       'mu'                = mu_from_vartheta(vartheta, p, proportions),
       'vartheta'          = vartheta,
       'shape'             = shape,
       'change_seed'       = change_seed,
       'precision_type'    = precision_type,
       'rho'               = rho,
       'band'              = band,
       'block_size'        = block_size,
       'Sigma'             = Sigma_obj$mat,
       'Sigma_inv'         = Sigma_obj$inverse,
       'locations'         = locations,
       'durations'         = durations,
       'proportions'       = proportions,
       'change_type'       = change_type,
       'changing_vars'     = changing_vars,
       'point_locations'   = point_locations,
       'point_proportions' = point_proportions,
       'point_mu'          = point_mu)
}

init_data_mc <- function(n = 1000, p = 10, vartheta = 1, shape = 6,
                         location = 300, duration = 10,
                         precision_type = 'banded', rho = 0.9,
                         point_anoms = FALSE,
                         band = 2, block_size = p) {
  locations <- 1:3 * location
  durations <- 3:1 * duration
  varthetas <- 1:3 * vartheta
  proportions <- c(1/p, round(sqrt(p)) / p, 3 * round(sqrt(p)) / p)
  change_types <- c("adjacent", "adjacent", "block_scattered")
  if (point_anoms) {
    point_locations <- round(c(5, 10, 50, 200, 350, 370, 660, 700, 800, 950) / 1000 * n)
    point_proportions <- rep(1/p, length(point_locations))
    point_mu <- rnorm(length(point_locations), 0, 4 * log(p))
  } else point_locations <- point_proportions <- point_mu <- NA

  init_data(n                 = n,
            p                 = p,
            proportions       = proportions,
            vartheta          = varthetas,
            shape             = shape,
            locations         = locations,
            durations         = durations,
            change_type       = change_types,
            point_locations   = point_locations,
            point_proportions = point_proportions,
            point_mu          = point_mu,
            precision_type    = precision_type,
            rho               = rho,
            band              = band,
            block_size        = block_size)
}

init_data_ <- function(data) {
  init_data(n                 = data$n,
            p                 = data$p,
            proportions       = data$proportions,
            vartheta          = data$vartheta,
            shape             = data$shape,
            change_seed       = data$change_seed,
            locations         = data$locations,
            durations         = data$durations,
            change_type       = data$change_type,
            changing_vars     = data$changing_vars,
            point_locations   = data$point_locations,
            point_proportions = data$point_proportions,
            point_mu          = data$point_mu,
            precision_type    = data$precision_type,
            rho               = data$rho,
            band              = data$band,
            block_size        = data$block_size)
}

#' @export
method_params <- function(cost = "cor", b = 1, minsl = 2, maxsl = 100,
                          precision_est_struct = "correct", est_band = NA,
                          size_mu = NA) {
  if (grepl("inspect", cost) && !is.na(precision_est_struct)) {
    if (is.na(est_band) || est_band != 0)
      precision_est_struct <- "correct"
  }
  if (grepl("iid", cost)) {
    precision_est_struct <- "banded"
    est_band <- 0
  } else {
    if (is.na(precision_est_struct) || precision_est_struct == "correct")
      est_band <- NA
    else if (precision_est_struct == "banded") {
      if (is.na(est_band)) est_band <- 2
    }
  }
  if (cost != "cor_dense_optimal") size_mu <- NA
  else {
    if (is.na(size_mu)) stop("size_mu must be specified if cost = cor_dense_optimal")
  }
  list("cost"                 = cost,
       "b"                    = b,
       "minsl"                = minsl,
       "maxsl"                = maxsl,
       "precision_est_struct" = precision_est_struct,
       "est_band"             = est_band,
       "size_mu"              = size_mu)
}

method_params_ <- function(method) {
  method_params(cost                 = method$cost,
                b                    = method$b,
                minsl                = method$minsl,
                maxsl                = method$maxsl,
                precision_est_struct = method$precision_est_struct,
                est_band             = method$est_band,
                size_mu              = method$size_mu)
}

#' @export
simulate_detection <- function(data = init_data(), method = method_params(),
                               seed = NULL, standardise_output = TRUE) {
  format_cpt_out <- function(cost, res) {
    if (res$value <= 0)
      return(data.table("cpt" = integer(0), "variate" = integer(0),
                        "mean1" = numeric(0), "mean2" = numeric(0)))
    else {
      if (cost == "sinspect") res$J <- which(res$proj != 0)
      mean1 <- colMeans(x$x[1:res$cpt, res$J, drop = FALSE])
      mean2 <- colMeans(x$x[(res$cpt + 1):nrow(x$x), res$J, drop = FALSE])
      return(data.table("cpt" = res$cpt, "variate" = res$J,
                        "mean1" = mean1, "mean2" = mean2))
    }
  }

  format_output <- function(cost, res) {
    if (cost == "cor_exact")
      return(list("collective" = collective_anomaliesR(res),
                  "point"      = point_anomaliesR(res)))
    else if (cost == "cor")
      return(list("collective" = collective_anomalies(list("anoms" = res)),
                  "point"      = point_anomalies(list("anoms" = res))))
    else if (cost == "iid")
      return(list("collective" = anomaly::collective_anomalies(res),
                  "point"      = anomaly::point_anomalies(res)))
    else if (cost %in% c("mvlrt", "sinspect"))
      return(format_cpt_out(cost, res))
    else if (cost == "inspect")
      return(anomalies_inspect(res, x$x))
  }

  if (!is.null(seed)) set.seed(seed)
  x <- simulate_cor_(data)
  if (method$cost == "iid") {
    x$x <- robust_scale(x$x)
    beta <- iid_penalty(data$n, data$p, method$b)
    beta_tilde <- iid_point_penalty(data$n, data$p, max(0.05, method$b))
    res <- anomaly::capa.mv(x$x,
                            beta        = beta,
                            beta_tilde  = beta_tilde,
                            min_seg_len = method$minsl,
                            max_seg_len = method$maxsl,
                            type        = "mean")
  } else {
    Q_hat <- get_Q_hat(x$x, data, method)
    if (grepl("cor", method$cost)) {
      x$x <- centralise(x$x)
      if (method$cost == "cor_exact") mvcapa_func <- mvcapa_cor_exact
      else if(method$cost == "cor")   mvcapa_func <- mvcapa_cor
      res <- mvcapa_func(x$x, Q_hat,
                         b                = method$b,
                         b_point          = max(0.05, method$b),
                         min_seg_len      = method$minsl,
                         max_seg_len      = method$maxsl)
    } else if (method$cost == "mvlrt")
      res <- single_mvnormal_changepoint(x$x, Q_hat,
                                         b = method$b,
                                         min_seg_len = method$minsl)
    else if (method$cost == "sinspect")
      res <- single_cor_inspect(t(x$x), Q_hat, method$b)
    else if (method$cost == "inspect")
      res <- cor_inspect(t(x$x), Q_hat, method$b)$changepoints
  }
  if (standardise_output) return(format_output(method$cost, res))
  else return(res)
}

#' @export
simulate_detection_known <- function(data = init_data(), method = method_params(),
                                     seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  x <- simulate_cor_(data)
  if (grepl("iid", method$cost)) {
    x$x <- robust_scale(x$x)
    x_anom <- x$x[(data$locations + 1):(data$locations + data$durations), ]
    if (method$cost == "iid") {
      penalty_vec <- cumsum(iid_penalty(data$n, data$p, method$b))
      return(optimise_mvnormal_iid_saving(x_anom, penalty_vec))
    } else if (method$cost == "iid_dense") {
      alpha <- get_penalty("dense", data$n, data$p, method$b)$alpha_const
      return(list(S_max = saving_iid(1:data$p, x_anom) - alpha))
    }
  } else if (grepl("cor", method$cost)) {
    Q_hat <- get_Q_hat(x$x, data, method)
    x$x <- centralise(x$x)
    x_anom <- x$x[(data$locations + 1):(data$locations + data$durations), ]
    if (method$cost == "cor") {
      penalty <- get_penalty('combined', data$n, data$p, method$b)
      lower_nbs <- lower_nbs(Q_hat)
      extended_nbs <- extended_lower_nbs(lower_nbs)
      return(optimise_mvnormal_saving(x_anom, Q_hat, lower_nbs, extended_nbs,
                                      penalty$alpha_const, penalty$beta,
                                      penalty$alpha_lin))
    } else if (method$cost == "cor_exact") {
      penalty <- get_penalty_vec("combined", data$n, data$p, method$b)
      return(optimise_mvnormal_saving_BF(x_anom, Q_hat, penalty, mu_MLE(Q_hat)))
    } else if (method$cost == "cor_BF") {
      penalty <- get_penalty_vec("combined", data$n, data$p, method$b)
      return(optimise_mvnormal_saving_BF(x_anom, Q_hat, penalty, mu_aMLE()))
    } else if (method$cost == "cor_dense") {
      alpha <- get_penalty("dense", data$n, data$p, method$b)$alpha_const
      mean_x <- matrix(colMeans(x_anom), nrow = data$p)
      return(list(S_max = dense_mvnormal_savings(mean_x, Q_hat, nrow(x_anom)) - alpha))
    } else if (method$cost == "cor_dense_optimal") {
      alpha <- get_penalty("dense", data$n, data$p, method$b)$alpha_const
      mu <- generate_change(method$size_mu, data$p, 0)
      return(list(S_max = optimal_mvnormal_dense_saving(x_anom, Q_hat, mu, alpha)))
    }
  } else if (method$cost == "sinspect") {
    Q_hat <- get_Q_hat(x$x, data, method)
    single_cor_inspect(t(x$x), Q_hat, method$b, data$locations, standardize.series = TRUE)
  } else if (method$cost == "mvlrt") {
    Q_hat <- get_Q_hat(x$x, data, method)
    return(optimise_mvnormal_lr(data$locations, x$x, Q_hat, method$b))
  }
}

test_optim <- function(data, method) {
  x <- simulate_cor_(data)
  Q_hat <- get_Q_hat(x$x, data, method)
  x$x <- centralise(x$x)
  x_anom <- x$x[(data$locations + 1):(data$locations + data$durations), ]
  penalty <- get_penalty('combined', data$n, data$p, method$b)
  lower_nbs <- lower_nbs(Q_hat)
  extended_nbs <- extended_lower_nbs(lower_nbs)
  optimise_mvnormal_saving(x_anom, Q_hat, lower_nbs, extended_nbs,
                           penalty$alpha_const, penalty$beta,
                           penalty$alpha_lin)
}

#' @export
test_runtime <- function(cost, b = 1, est_band = 2, times = 5) {
  data <- init_data(p = 100, n = 200, vartheta = 0, rho = 0.9)
  method <- method_params(cost, b, precision_est_struct = "banded", est_band = est_band)
  print(microbenchmark::microbenchmark(
    simulate_detection(data, method, standardise_output = TRUE),
    times = times
  ))
}

#' @export
simple_example <- function(p = 4, n = 200, vartheta = 10, method = "cor",
                           b = 2, b_point = 0.5, v = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (v == 1) {
    true_anoms <- data.frame(start = c(20, n - 50), duration = c(round(n/10), 30))
    true_anoms$end <- true_anoms$start + true_anoms$duration
    point_locs <- c(round(n / 2)) + c(-5, 0, 5)
    x <- simulate_cor(n = n, p = p,
                      locations = true_anoms$start,
                      durations = true_anoms$duration,
                      proportions = c(2/p, 1/p),
                      vartheta = vartheta, change_type = "random",
                      point_locations = point_locs, point_proportions = 1/p,
                      point_mu = c(1.5 * vartheta, - vartheta, vartheta),
                      Sigma = car_precision_mat(banded_neighbours(3, p), rho = 0.9))$x
  } else if (v == 2) {
    true_anoms <- data.frame(start = c(50, round(n / 3), n - 100),
                             duration = c(round(n / 20), round(n / 40), round(n / 10)))
    true_anoms$end <- true_anoms$start + true_anoms$duration
    point_locs <- c(200, 210, 215, 260, 400, 406, 460, 630, 700, 705, 712, 800)
    x <- simulate_cor(n = 1000, p = p,
                      locations = true_anoms$start,
                      durations = true_anoms$duration,
                      proportions = c(2/p, 1, 1/p),
                      vartheta = c(0.8, 4, 0.8), change_type = "random",
                      point_locations = point_locs,
                      point_mu = 4, point_proportions = 2/p,
                      Sigma = constant_cor_mat(p, 0.5))$x
  }
  true_anoms <- rbind(true_anoms, data.frame(start = point_locs - 1, duration = 1, end = point_locs))
  if (method == "cor") {
    x <- centralise(x)
    Q_hat <- estimate_precision_mat(x, adjacency_mat(banded_neighbours(4, p), sparse = FALSE))
    res <- mvcapa_cor(x, Q_hat, b = b, b_point = b_point)
    plot_capa(list("x" = x, "anoms" = res), true_anoms = true_anoms)
  } else if (method == "iid") {
    beta <- iid_penalty(n, p, b)
    beta_tilde <- iid_point_penalty(n, p, b_point)
    res <- anomaly::capa.mv(x,
                            beta        = beta,
                            beta_tilde  = beta_tilde,
                            type        = "mean")
    plot_capa(res, cost = "iid", true_anoms = true_anoms)
  }
}

save_simple_example <- function(p = 4, n = 200, vartheta = 10,
                                b = 2, b_point = 2, seed = NULL) {
  file_name <- paste0("simple_example",
                      "_p", p,
                      "_vartheta", vartheta,
                      "_b", b,
                      "_bpoint", b_point,
                      ".png")
  show(simple_example(p, n, vartheta, b, b_point, seed))
  ggsave(paste0("./images/", file_name), width = 7, height = 5, units = "in")
}

simple_examples_presentation <- function() {
  save_simple_example(vartheta = 0, b = 10^6, b_point = 10^6, seed = 6)
  save_simple_example(b = 10^6, b_point = 10^6, seed = 6)
  save_simple_example(seed = 6)
}

simple_examples_paper <- function(p = 10, save = FALSE) {
  plots <- list(simple_example(p = p, n = 1000, method = "cor", b = 1, b_point = 0.5,
                               vartheta = 1, v = 2, seed = 103),
                simple_example(p = p, n = 1000, method = "iid", b = 1.68, b_point = 0.5,
                               vartheta = 1, v = 2, seed = 103))
  figure <- ggpubr::ggarrange(plotlist = plots, nrow = 1, ncol = 2,
                              common.legend = TRUE,
                              legend = "bottom")
  if (save) {
    show(figure)
    ggsave(paste0("./images/simple_example", p, ".png"), width = 8, height = 5, units = "in")
  } else return(figure)
}



