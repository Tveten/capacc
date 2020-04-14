#' @export
init_data <- function(n = 200, p = 10, proportions = sqrt(p)/p, mu = 1,
                            locations = n - durations - 1, durations = 10,
                            change_type = 'adjacent', point_locations = NA,
                            point_proportions = NA, point_mu = NA,
                            precision_type = 'banded', rho = 0.8, band = 2,
                            block_size = p,
                            min_nbs = 1, max_nbs = 3) {
  get_Sigma <- function(precision_type) {
    if (precision_type == 'iid') {
      return(list('mat'     = diag(1, p),
                  'inverse' = diag(1, p)))
    } else if (precision_type == 'ar1') {
      return(list('mat'     = ar_cor_mat(p, rho),
                  'inverse' = ar_precision_mat(p, rho)))
    } else {
      if (precision_type == 'lattice')
        precision_mat <- car_precision_mat(lattice_neighbours(p), rho)
      else if (precision_type == 'banded')
        precision_mat <- car_precision_mat(banded_neighbours(band, p), rho)
      else if (precision_type == 'random')
        precision_mat <- car_precision_mat(random_neighbours(p), rho, min_nbs = min_nbs, max_nbs = max_nbs)
      else if (precision_type == "block_banded")
        precision_mat <- block_precision_mat(p, block_size, within_block_type = "banded", rho, band = band)
      else if (precision_type == "block_lattice")
        precision_mat <- block_precision_mat(p, block_size, within_block_type = "lattice", rho)
      return(list('mat'     = solve(precision_mat),
                  'inverse' = precision_mat))
    }
  }

  Sigma_obj <- get_Sigma(precision_type)
  list('n'                 = n,
       'p'                 = p,
       'mu'                = mu,
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
       'point_locations'   = point_locations,
       'point_proportions' = point_proportions,
       'point_mu'          = point_mu)
}

#' @export
mvcapa_params <- function(cost = "cor", b = 1, minsl = 2, maxsl = 100,
                          precision_est_struct = "correct", est_band = 2) {
  if (precision_est_struct == "correct") est_band <- NA
  list("cost"                 = cost,
       "b"                    = b,
       "minsl"                = minsl,
       "maxsl"                = maxsl,
       "precision_est_struct" = precision_est_struct,
       "est_band"             = est_band)
}

#' @export
simulate_mvcapa <- function(data = init_data(), params = mvcapa_params(),
                            seed = NULL, return_anom_only = FALSE) {
  get_adj_mat <- function(est_struct) {
    if (est_struct == "correct")
      return(data$Sigma_inv)
    else if (est_struct == "banded")
      return(adjacency_mat(banded_neighbours(params$est_band, data$p)))
  }

  if (!is.null(seed)) set.seed(seed)
  x <- simulate_cor(n                 = data$n,
                    p                 = data$p,
                    mu                = data$mu,
                    Sigma             = data$Sigma,
                    locations         = data$locations,
                    durations         = data$durations,
                    proportions       = data$proportions,
                    change_type       = data$change_type,
                    point_locations   = data$point_locations,
                    point_proportions = data$point_proportions,
                    point_mu          = data$point_mu)
  x <- anomaly::robustscale(x)
  if (params$cost == "cor") {
    Q_hat <- estimate_precision_mat(x, get_adj_mat(params$precision_est_struct))
    res <- mvcapa_cor(x, Q_hat,
                      b                = params$b,
                      b_point          = max(0.05, params$b),
                      min_seg_len      = params$minsl,
                      max_seg_len      = params$maxsl)
    if (!return_anom_only) return(res)
    else return(list("collective" = collective_anomalies(list("anoms" = res)),
                     "point"      = point_anomalies(list("anoms" = res))))
  } else if (params$cost == "iid") {
    beta <- iid_penalty(data$n, data$p, params$b)
    beta_tilde <- iid_point_penalty(data$n, data$p, max(0.05, params$b))
    res <- anomaly::capa.mv(x,
                            beta = beta,
                            beta_tilde = beta_tilde,
                            min_seg_len = params$minsl,
                            max_seg_len = params$maxsl,
                            type        = "mean")
    if (!return_anom_only) return(res)
    else return(list("collective" = anomaly::collective_anomalies(res),
                     "point"      = anomaly::point_anomalies(res)))
  }
}

change_sd <- function(p, prop, duration) {
  xi <- - log(prop) / log(p)
  if (xi >= 0 && xi <= 1/2) mu_boundary <- p^(- 1/2 + xi)
  else if (xi > 1/2 && xi <= 3/4) mu_boundary <- sqrt((2 * xi - 1) * log(p))
  else if (xi > 3/4 && xi <= 1) mu_boundary <- sqrt(2 * log(p)) * (1 - sqrt(1 - xi))

  # > pnorm(1/8) - pnorm(-1/8)
  # [1] 0.09947645
  # I.e., 90% of drawn changes will be greater than the detection boundary
  # (when averaging over durations).
  mu_boundary / sqrt(duration)
}

draw_anomaly_data <- function(data, mean_duration = 10,
                               anomaly_rate = 3 / data$n,
                               change_sd = log(data$p), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  locations <- list(rnbinom(1, 1, anomaly_rate) + 1)
  durations <- list(rpois(1, mean_duration))
  if (locations[[1]] + durations[[1]] >= data$n)
    locations[[1]] <- data$n - durations[[1]] - 1
  mus <- list(rnorm(1, sd = change_sd(data$p, data$proportions, mean_duration)))
  i <- 1
  while (locations[[i]] + durations[[i]] < data$n) {
    locations[[i + 1]] <- locations[[i]] + durations[[i]] + rnbinom(1, 1, anomaly_rate)
    durations[[i + 1]] <- rpois(1, mean_duration)
    mus[[i + 1]] <- rnorm(1, sd = change_sd(data$p, data$proportions, mean_duration))
    i <- i + 1
  }
  n_anoms <- length(locations) - 1
  data$locations <- unlist(locations)[1:n_anoms]
  data$durations <- unlist(durations)[1:n_anoms]
  data$mu <- unlist(mus)[1:n_anoms]
  data
}
