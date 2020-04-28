
get_sim_seeds <- function(params_list, variables) {
  # Uses structure of output from expand.grid.
  if (names(variables)[1] != "cost")
    stop("'cost' must be the first element in 'variables' for seeds to be correct.")
  if (!is.null(variables$precision_est_struct) && names(variables)[2] != "precision_est_struct")
    stop("'precision_est_struct' must be the second element in 'variables' if present for seeds to be correct.")
  n_seeds <- length(params_list[[1]]) / (length(variables$cost) * length(variables$precision_est_struct))
  seeds <- sample(1:10^6, n_seeds)
  rep(seeds, each = length(variables$cost))
}


mu_from_vartheta <- function(vartheta, p, prop) {
  vartheta / sqrt(round(prop * p))
}

vartheta_from_mu <- function(mu, p, prop) {
  mu * sqrt(round(prop * p))
}

split_params <- function(a, grouping) {

  get_init_func <- function(name) {
    if      (name == "data") return(init_data_)
    else if (name == "method") return(method_params_)
    else                      return(identity)
  }
  split_list <- lapply(grouping, extract_nested_elements, list_of_lists = a)
  names(split_list) <- names(grouping)
  out <- list()
  for (i in 1:length(split_list)) {
    init_func <- get_init_func(names(grouping)[i])
    out[[length(out) + 1]] <- Map(init_func, split_list[[i]])
  }
  names(out) <- names(grouping)
  out
}

#
# change_sd <- function(p, prop, duration) {
#   xi <- - log(prop) / log(p)
#   if (xi >= 0 && xi <= 1/2) mu_boundary <- p^(- 1/2 + xi)
#   else if (xi > 1/2 && xi <= 3/4) mu_boundary <- sqrt((2 * xi - 1) * log(p))
#   else if (xi > 3/4 && xi <= 1) mu_boundary <- sqrt(2 * log(p)) * (1 - sqrt(1 - xi))
#
#   # > pnorm(1/8) - pnorm(-1/8)
#   # [1] 0.09947645
#   # I.e., 90% of drawn changes will be greater than the detection boundary
#   # (when averaging over durations).
#   mu_boundary / sqrt(duration)
# }
#
# draw_anomaly_data <- function(data, mean_duration = 10,
#                                anomaly_rate = 3 / data$n,
#                                change_sd = log(data$p), seed = NULL) {
#   if (!is.null(seed)) set.seed(seed)
#
#   locations <- list(rnbinom(1, 1, anomaly_rate) + 1)
#   durations <- list(rpois(1, mean_duration))
#   if (locations[[1]] + durations[[1]] >= data$n)
#     locations[[1]] <- data$n - durations[[1]] - 1
#   mus <- list(rnorm(1, sd = change_sd(data$p, data$proportions, mean_duration)))
#   i <- 1
#   while (locations[[i]] + durations[[i]] < data$n) {
#     locations[[i + 1]] <- locations[[i]] + durations[[i]] + rnbinom(1, 1, anomaly_rate)
#     durations[[i + 1]] <- rpois(1, mean_duration)
#     mus[[i + 1]] <- rnorm(1, sd = change_sd(data$p, data$proportions, mean_duration))
#     i <- i + 1
#   }
#   n_anoms <- length(locations) - 1
#   data$locations <- unlist(locations)[1:n_anoms]
#   data$durations <- unlist(durations)[1:n_anoms]
#   data$mu <- unlist(mus)[1:n_anoms]
#   data
# }

