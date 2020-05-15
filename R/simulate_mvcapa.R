
#' @export
init_data <- function(n = 100, p = 10, proportions = sqrt(p)/p,
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

  if (change_type == "custom") {
    if (!is.na(changing_vars))
      proportions <- length(changing_vars) / p
    else
      stop("If change_type is 'custom', you must provide a changing_vars vector")
  }

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
  if (grepl("inspect", cost)) b <- NA
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


get_adj_mat <- function(data, method) {
  if (method$precision_est_struct == "correct")
    return(data$Sigma_inv)
  else if (method$precision_est_struct == "banded")
    return(adjacency_mat(banded_neighbours(method$est_band, data$p)))
}

get_Q_hat <- function(x, data, method) {
  if (is.na(method$precision_est_struct))
    return(data$Sigma_inv)
  else
    return(estimate_precision_mat(x, get_adj_mat(data, method)))
}

#' @export
simulate_mvcapa <- function(data = init_data(), method = method_params(),
                            seed = NULL, return_anom_only = TRUE) {

  get_anom_only_list <- function(cost, res) {
    if (cost == "cor_exact")
      return(list("collective" = collective_anomaliesR(res),
                  "point"      = point_anomaliesR(res)))
    else if (cost == "cor")
      return(list("collective" = collective_anomalies(list("anoms" = res)),
                  "point"      = point_anomalies(list("anoms" = res))))
    else if (cost == "iid")
      return(list("collective" = anomaly::collective_anomalies(res),
                  "point"      = anomaly::point_anomalies(res)))
  }

  if (!is.null(seed)) set.seed(seed)
  x <- simulate_cor_(data)
  x$x <- anomaly::robustscale(x$x)
  if (grepl("cor", method$cost)) {
    Q_hat <- get_Q_hat(x$x, data, method)
    if (method$cost == "cor_exact") mvcapa_func <- mvcapa_cor_exact
    else if(method$cost == "cor")   mvcapa_func <- mvcapa_cor
    res <- mvcapa_func(x$x, Q_hat,
                       b                = method$b,
                       b_point          = max(0.05, method$b),
                       min_seg_len      = method$minsl,
                       max_seg_len      = method$maxsl)
  } else if (method$cost == "iid") {
    beta <- iid_penalty(data$n, data$p, method$b)
    beta_tilde <- iid_point_penalty(data$n, data$p, max(0.05, method$b))
    res <- anomaly::capa.mv(x$x,
                            beta        = beta,
                            beta_tilde  = beta_tilde,
                            min_seg_len = method$minsl,
                            max_seg_len = method$maxsl,
                            type        = "mean")
  }
  if (return_anom_only) return(get_anom_only_list(method$cost, res))
  else return(res)
}

#' @export
simulate_mvcapa_known <- function(data = init_data(), method = method_params(),
                                  seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  x <- simulate_cor_(data)
  # print(x$mu)
  x$x <- anomaly::robustscale(x$x)
  x_anom <- x$x[(data$locations + 1):(data$locations + data$durations), ]
  if (method$cost == "iid") {
    penalty_vec <- cumsum(iid_penalty(data$n, data$p, method$b))
    return(optimise_mvnormal_iid_saving(x_anom, penalty_vec))
  } else if (method$cost == "iid_dense") {
    alpha <- get_penalty("dense", data$n, data$p, method$b)$alpha_const
    return(list(S_max = saving_iid(1:data$p, x_anom) - alpha))
  } else if (grepl("cor", method$cost)) {
    Q_hat <- get_Q_hat(x$x, data, method)
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
      return(optimise_mvnormal_saving_BF(x_anom, Q_hat, penalty, mu_aMLE))
    } else if (method$cost == "cor_dense") {
      alpha <- get_penalty("dense", data$n, data$p, method$b)$alpha_const
      mean_x <- matrix(colMeans(x_anom), nrow = data$p)
      return(list(S_max = dense_mvnormal_savings(mean_x, Q_hat, nrow(x_anom)) - alpha))
    } else if (method$cost == "cor_dense_optimal") {
      alpha <- get_penalty("dense", data$n, data$p, method$b)$alpha_const
      mu <- generate_change(method$size_mu, data$p, 0)
      # mu <- as.vector(x$mu)
      return(list(S_max = optimal_mvnormal_dense_saving(x_anom, Q_hat, mu, alpha)))
    }
  }
}

