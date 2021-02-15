#' @export
init_data <- function(n = 100, p = 10, proportions = round(sqrt(p))/p,
                      vartheta = 1, shape = 0, change_seed = NA,
                      locations = 50, durations = 10,
                      change_type = 'adjacent', changing_vars = NA,
                      point_locations = NA, point_proportions = NA,
                      point_mu = NA, precision_type = 'banded', rho = 0.9,
                      band = 2, block_size = p, min_nbs = 1, max_nbs = 3,
                      n_sd_changes = 0) {
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
       'point_mu'          = point_mu,
       'n_sd_changes'      = n_sd_changes)
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
            block_size        = data$block_size,
            n_sd_changes      = data$n_sd_changes)
}

init_data_mc <- function(n = 1000, p = 10, vartheta = 1, shape = 6,
                         location = 300, duration = 10,
                         precision_type = 'banded', rho = 0.9,
                         point_anoms = FALSE, band = 2, block_size = p,
                         n_sd_changes = 0) {
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
            block_size        = block_size,
            n_sd_changes      = n_sd_changes)
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
  }
  if (cost == "gflars") {
    b <- minsl <- maxsl <- precision_est_struct <- est_band <- size_mu <- NA
  }
  if (cost == "var_pgl") {
    precision_est_struct <- est_band <- size_mu <- NA
  }
  if (cost == "decor")
    precision_est_struct <- est_band <- size_mu <- NA
  if (is.na(precision_est_struct) || precision_est_struct == "correct")
    est_band <- NA
  else if (precision_est_struct == "banded") {
    if (is.na(est_band)) est_band <- 2
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
tuning_params <- function(alpha = 0.05, tol = 0.02, max_iter = 50,
                          init_b = c(0.05, 1, 10), n_sim = 200) {
  list("alpha"           = alpha,
       "alpha_tol"       = tol,
       "tuning_max_iter" = max_iter,
       "init_b"          = init_b,
       "tuning_n_sim"    = n_sim)
}

#' @export
curve_params <- function(max_dist = 0.2, max_iter = 50, n_sim = 100,
                         init_values = c(0.1, 5)) {
  list("curve_max_dist" = max_dist,
       "curve_max_iter" = max_iter,
       "curve_n_sim"    = n_sim,
       "init_values"    = init_values)
}
