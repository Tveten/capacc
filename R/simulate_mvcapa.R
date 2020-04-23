
#' @export
init_data <- function(n = 100, p = 10, proportions = sqrt(p)/p, mu = NA,
                      vartheta = 1, shape = 0, change_seed = NA,
                      locations = 50, durations = 10,
                      change_type = 'adjacent', changing_vars = NA,
                      point_locations = NA, point_proportions = NA,
                      point_mu = NA, precision_type = 'banded', rho = 0.9,
                      band = 2, block_size = p, min_nbs = 1, max_nbs = 3) {
  get_Sigma <- function(precision_type) {
    if (precision_type == 'iid') {
      return(list('mat'     = diag(1, p),
                  'inverse' = diag(1, p)))
    } else if (precision_type == 'ar1') {
      return(list('mat'     = ar_cor_mat(p, rho),
                  'inverse' = ar_precision_mat(p, rho)))
    } else {
      precision_mat <- block_precision_mat(p, block_size,
                                           within_block_type = precision_type,
                                           band = band, min_nbs = min_nbs,
                                           max_nbs = max_nbs)
      return(list('mat'     = solve(precision_mat),
                  'inverse' = precision_mat))
    }
  }

  if (all(is.na(c(mu, vartheta))))
    mu <- vartheta <- 0
  else if (!is.na(mu) && is.na(vartheta))
    vartheta <- vartheta_from_mu(mu, p, proportions)
  else
    mu <- mu_from_vartheta(vartheta, p, proportions)

  if (change_type == "custom") {
    if (!is.na(changing_vars))
      proportions <- length(changing_vars) / p
    else
      stop("If change_type is 'custom', you must provide a changing_vars vector")
  }

  Sigma_obj <- get_Sigma(precision_type)
  if (precision_type == "lattice") band <- band(Sigma_obj$inverse)

  list('n'                 = n,
       'p'                 = p,
       'mu'                = mu,
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
                          precision_est_struct = "correct", est_band = NA) {
  if (grepl("iid", cost)) {
    precision_est_struct <- "banded"
    est_band <- 0
  } else {
    if (precision_est_struct == "correct") est_band <- NA
    else if (precision_est_struct == "banded") {
      if (is.na(est_band)) est_band <- 2
    }
  }
  if (grepl("inspect", cost)) b <- NA
  list("cost"                 = cost,
       "b"                    = b,
       "minsl"                = minsl,
       "maxsl"                = maxsl,
       "precision_est_struct" = precision_est_struct,
       "est_band"             = est_band)
}

method_params_ <- function(method) {
  method_params(cost                 = method$cost,
                b                    = method$b,
                minsl                = method$minsl,
                maxsl                = method$maxsl,
                precision_est_struct = method$precision_est_struct,
                est_band             = method$est_band)
}


#' @export
simulate_mvcapa <- function(data = init_data(), params = method_params(),
                            seed = NULL, return_anom_only = FALSE) {
  get_adj_mat <- function(est_struct) {
    if (est_struct == "correct")
      return(data$Sigma_inv)
    else if (est_struct == "banded")
      return(adjacency_mat(banded_neighbours(params$est_band, data$p)))
  }

  if (!is.null(seed)) set.seed(seed)
  x <- anomaly::robustscale(simulate_cor_(data))
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
