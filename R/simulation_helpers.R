get_adj_mat <- function(data, method) {
  if (method$precision_est_struct == "correct")
    return(as.matrix(data$Sigma_inv))
  else if (method$precision_est_struct == "banded")
    return(adjacency_mat(banded_neighbours(method$est_band, data$p), sparse = FALSE))
}

get_Q_hat <- function(x, data, method) {
  if (is.na(method$precision_est_struct))
    return(data$Sigma_inv)
  else if (method$cost %in% c("mvlrt", "inspect", "sinspect")) {
    diff_x <- x[2:nrow(x), ] - x[1:(nrow(x) - 1), ]
    if (grepl("inspect", method$cost) && method$precision_est_struct == "correct")
      return(2 * solve(robust_cov_mat(diff_x)))
    else
      return(2 * robust_sparse_precision(diff_x, get_adj_mat(data, method)))
  } else
    return(robust_sparse_precision(x, get_adj_mat(data, method)))
}

robust_scale <- function(x, Q = NULL) {
  if (is.null(Q)) sigma <- Rfast::colMads(x)
  else sigma <- sqrt(diag(solve(Q)))
  med <- Rfast::colMedians(x)
  t((t(x) - med) / sigma)
}

robust_whitening <- function(x) {
  # x <- apply(x, 2, function(y) y - mean(y))
  # Sigma <- 1 / (nrow(x) - 1) * t(x) %*% x
  x <- centralise(x)
  Sigma <- robust_cov_mat(x)
  eigen_S <- eigen(Sigma, symmetric = TRUE)
  Sigma_inv_root <- eigen_S$vectors %*% diag(1 / sqrt(eigen_S$values)) %*% t(eigen_S$vectors)
  x %*% Sigma_inv_root
}

#' @export
centralise <- function(x) {
  med <- Rfast::colMedians(x)
  t(t(x) - med)
}

mu_from_vartheta <- function(vartheta, p, prop) {
  vartheta / sqrt(round(prop * p))
}

vartheta_from_mu <- function(mu, p, prop) {
  mu * sqrt(round(prop * p))
}

# get_sim_seeds <- function(params_list, variables) {
#   # Uses structure of output from expand.grid.
#   same_seed_names <- c("cost", "precision_est_struct", "est_band")
#   variable_names <- names(variables)
#   ind_names_present <- which(variable_names %in% same_seed_names)
#   if (length(ind_names_present) != max(ind_names_present))
#     stop("cost, precision_est_struct and est_band must be the first elements of 'variables' for seeds to be correct")
#   n_costs <- Reduce(`*`, vapply(variables[ind_names_present], length, numeric(1)))
#   n_settings <- length(params_list[[1]])
#   n_unique_seeds <- n_settings / n_costs
#   rep(sample(1:10^6, n_unique_seeds), each = n_costs)
# }

get_sim_seeds <- function(vars, prune = TRUE) {
  if (prune) var_grid <- cost_pruned_expand_grid(vars)
  else       var_grid <- as.data.table(expand.grid(vars, stringsAsFactors = FALSE))
  cost_vars <- c("cost", "precision_est_struct", "est_band")
  var_grid <- var_grid[, "seed" := sample(1:10^6, 1),
                       by = c(names(var_grid)[!names(var_grid) %in% cost_vars])]
  var_grid$seed
}

get_previous_seed <- function(seed, all_params, read_func, out_file) {
  # all_params <- c(data, method, tuning, curve, list(loc_tol = loc_tol))
  all_res <- read_results(out_file)
  exclude <- c("cost", "precision_est_struct", "est_band")
  # prev_res <- read_anom_class(all_res, all_params, exclude = exclude)
  prev_res <- read_func(all_res, all_params, exclude = exclude)
  prev_seed <- unique(prev_res$seed)
  if (length(prev_seed) == 0) return(seed)
  else if (length(prev_seed) >= 1) {
    if (length(prev_seed) >= 2) {
      message("WARNING: More than 1 seed used for this simulation setting. Using the first.")
      prev_res[, print(paste0("cost and est_structs for seed ", .GRP, ": ",
                              paste(unique(c(cost, precision_est_struct, est_band)), collapse = " "))),
               by = seed]
    }
    message(paste0("Overriding input seed with the same seed used for other costs: ",
                   paste(prev_seed[1], collapse = ", ")))
    return(prev_seed[1])
  }
}


expand_list <- function(a, vars, prune = TRUE) {
  if (prune) var_grid <- as.data.frame(cost_pruned_expand_grid(vars))
  else       var_grid <- expand.grid(vars, stringsAsFactors = FALSE)
  var_names <- names(var_grid)
  # print(var_grid)
  # var_grid <- var_grid[c(3, 1:2, 4:nrow(var_grid)), ]
  # print(var_grid)

  # For loop to avoid "c" being used on each row of var_grid, which might change
  # the type of variables.
  out_list <- list()
  for (i in 1:nrow(var_grid)) {
    a_copy <- a
    for (j in 1:ncol(var_grid)) {
      a_copy[[var_names[j]]] <- unname(var_grid[i, j])
    }
    out_list[[length(out_list) + 1]] <- a_copy
  }
  out_list
}

cost_pruned_expand_grid <- function(vars) {
  var_grid <- as.data.table(expand.grid(vars, stringsAsFactors = FALSE))
  if (any(names(var_grid) == "cost") && length(var_grid) == 1) return(var_grid)
  cost_vars <- c("cost", "precision_est_struct", "est_band")
  var_grid <- var_grid[, rbind(.SD[cost == "iid"][1],
                   .SD[cost == "decor"][1],
                   .SD[cost == "gflars"][1],
                   .SD[cost == "var_pgl"][1],
                   .SD[cost == "cor" & is.na(precision_est_struct)][1],
                   .SD[cost == "cor" & precision_est_struct == "correct"][1],
                   .SD[cost != "iid" & cost != "decor" & cost != "gflars" & cost != "var_pgl" &
                         !(cost == "cor" & is.na(precision_est_struct)) &
                         !(cost == "cor" & precision_est_struct == "correct")]),
                   by = c(names(var_grid)[!names(var_grid) %in% cost_vars])]
  if (any(names(var_grid) == "est_band"))
    var_grid <- var_grid[, .SD[!(cost == "cor" & precision_est_struct == "banded" & est_band == 0)],
                         by = c(names(var_grid)[!names(var_grid) %in% cost_vars])]
  var_grid <- unique(var_grid)
  var_grid[!is.na(cost)]
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

#' @export
pll_test <- function(cpus = 1) {
  test_func <- function(i) {
    Sys.sleep(0.1)
    i
  }

  if (cpus == 1) {
    out <- Map(test_func, 1:10)
  } else {
    parallelMap::parallelStart(mode = "multicore", cpus = cpus, show.info = FALSE)
    out <- parallelMap::parallelMap(
      test_func,
      1:10
    )
    parallelMap::parallelStop()
  }
  out
}

anomalies_from_cpt <- function(cpt, x, tol = 1) {
  if (length(cpt) == 0) {
    return(list("collective" = data.table(start = integer(0), end = integer(0)),
                "point"      = data.table(location = integer(0))))
  }
  n <- nrow(x)
  res <- data.table(location = c(0, cpt, n))
  res$mean_size <- 0
  for (i in 2:length(res$location)) {
    seg_mean <- colMeans(x[(res$location[i - 1] + 1):res$location[i], , drop = FALSE])
    res$mean_size[i] <- sign(sum(seg_mean)) * sqrt(sum(seg_mean^2))
  }
  starts <- integer(0)
  ends <- integer(0)
  in_anom <- FALSE
  curr_start_ind <- 0
  i <- 2
  while (i <= nrow(res)) {
    if (!in_anom && abs(res$mean_size[i]) >= tol) {
      curr_start_ind <- i - 1
      starts <- c(starts, res$location[curr_start_ind] + 1)
      in_anom <- TRUE
      i <- i + 1
    } else {
      if (in_anom) {
        end_anom <- is_in_interval(res$mean_size[i], c(-tol, tol))
        switch_anom <- (res$mean_size[curr_start_ind + 1] < 0 && res$mean_size[i] > tol) ||
          (res$mean_size[curr_start_ind + 1] > 0 && res$mean_size[i] < - tol)
        if (end_anom) {
          ends <- c(ends, res$location[i - 1])
          in_anom <- FALSE
        } else if (switch_anom) {
          ends <- c(ends, res$location[i - 1])
          curr_start_ind <- i - 1
          starts <- c(starts, res$location[curr_start_ind] + 1)
        }
      }
      i <- i + 1
    }
  }
  if (in_anom) ends <- c(ends, res[.N, location])
  if (length(starts) != length(ends)) {
    print(starts)
    print(ends)
    stop("Bug when extracting inspect anomalies. Unequal number of start and end points.")
  }
  anoms <- data.table(start = starts, end = ends)
  return(list("collective" = anoms[start != end],
              "point"      = data.table(location = anoms[start == end, start])))
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

