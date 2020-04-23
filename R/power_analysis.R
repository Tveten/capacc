#' @export
curve_params <- function(max_dist = 0.2, max_iter = 50, n_sim = 100,
                         init_values = c(0.1, 5)) {
  list("curve_max_dist" = max_dist,
       "curve_max_iter" = max_iter,
       "curve_n_sim"    = n_sim,
       "init_values"    = init_values)
}

#' @export
power_curve <- function(out_file, data = init_data(), params = method_params(),
                        tuning = tuning_params(), curve = curve_params(),
                        loc_tol = 10, seed = NA) {
  add_setup_info <- function(res) {
    which_data <- !(grepl("Sigma", names(data)) | names(data) %in% c("changing_vars", "vartheta"))
    res <- cbind(res,
                 as.data.table(data[which_data]),
                 as.data.table(params[names(params) != "b"]),
                 as.data.table(tuning[names(tuning) != "init_b"]),
                 as.data.table(curve[names(curve) != "init_values"]))
    res$mu <- mu_from_vartheta(res$vartheta, data$p, data$proportions)
    res$loc_tol <- loc_tol
    res$seed <-  seed
    res
  }

  est_power <- function(vartheta) {
    cat('.')
    power <- mean(unlist(lapply(1:curve$curve_n_sim, function(i) {
      data$vartheta <- vartheta
      count_tfp_anom(simulate_mvcapa(data, params, return_anom_only = TRUE),
                     loc_tol, data)$tp
    })))
    data.table("vartheta" = vartheta, "power" = power)
  }

  init_power_est <- function(init_values) {
    res <- do.call('rbind', Map(est_power, curve$init_values))
    while ((min(res$power) > 0.02 || max(res$power) < 0.98) && max(res$vartheta) <= 10) {
      if (min(res$power) > 0.02)
        res <- rbind(res, est_power(min(res$vartheta) / 1.2))
      if (max(res$power < 0.98))
        res <- rbind(res, est_power(max(res$vartheta) * 1.2))
    }
    res
  }

  split_inds <- function(res) {
    # exceeds_max_dist <- adjacent_dist(res[, .(vartheta, power)]) > curve$curve_max_dist
    exceeds_max_dist <- diff(res$power) > curve$curve_max_dist
    ind <- which(exceeds_max_dist) + 1
    ind <- ind[!(res[ind - 1]$power >= 0.98 & res[ind]$power >= 0.98)]
    ind <- ind[!(res[ind - 1]$power <= 0.02 & res[ind]$power <= 0.02)]
    ind
  }

  message(paste0("Estimating power curve for n=", data$n,
                 ", p=", data$p,
                 ", cost=", params$cost,
                 ", precision=", data$precision_type,
                 ", block_size=", data$block_size,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", prop=", data$proportions,
                 ", shape=", data$shape,
                 ", precision_est_struct=", params$precision_est_struct,
                 ", est_band=", params$est_band,
                 "."))
  all_params <- c(data, params, tuning, curve, list("loc_tol" = loc_tol))
  if (already_estimated(out_file, all_params, read_power_curve)) return(NULL)

  params$b <- get_tuned_penalty(data, params, tuning, seed = seed + 2)$b
  if (!is.na(seed)) set.seed(seed)
  res <- init_power_est(curve$init_values)
  ind <- split_inds(res)
  while (length(ind) > 0 && nrow(res) <= curve$curve_max_iter) {
    new_varthetas <- (res$vartheta[ind - 1] + res$vartheta[ind]) / 2
    res <- rbind(res, do.call('rbind', Map(est_power, new_varthetas)))
    res <- res[order(vartheta)]
    ind <- split_inds(res)
  }
  cat('\n')
  print(res)
  fwrite(add_setup_info(res), paste0("./results/", out_file), append = TRUE)
}

many_power_curves <- function(out_file = "power.csv",
                              variables = list("cost" = c("iid", "cor"),
                                               "rho" = c(0.7, 0.9, 0.99),
                                               "proportions" = c(0.1, 0.3, 1)),
                              data   = init_data(),
                              params = method_params(),
                              tuning = tuning_params(),
                              curve  = curve_params(max_dist = 0.1, n_sim = 300),
                              loc_tol = 10) {
  params_list <- split_lists(expand_list(c(data, params), variables),
                             list("data"   = names(data),
                                  "method" = names(params)))
  data_list <- lapply(params_list$data, function(data) {init_data_(data)})
  method_list <- lapply(params_list$method, function(method) {method_params_(method)})
  Map(power_curve,
      data   = data_list,
      params = method_list,
      seed   = get_sim_seeds(params_list, variables),
      MoreArgs = list("out_file"     = out_file,
                      "tuning"        = tuning,
                      "curve"         = curve,
                      "loc_tol"       = loc_tol))
  NULL
}

#' @export
power_runs <- function() {
  curve <- curve_params(max_dist = 0.1, n_sim = 300)
  out_file <- "power.csv"

  #### BANDED
  banded_data <- init_data(n = 100, p = 10, precision_type = "banded",
                           band = 2, locations = 50, durations = 10,
                           change_type = "adjacent")
  banded_variables <- list("cost"        = c("iid", "cor"),
                           "rho"         = rev(c(-0.3, 0.5, 0.7, 0.9, 0.99)),
                           "proportions" = rev(c(0.1, 0.3, 1)),
                           "shape"       = rev(c(0, 2, 3, 4, 5)))
  many_power_curves(out_file, banded_variables, banded_data,
                    method_params(), tuning_params(), curve)

  #### BLOCK
  block_data <- init_data(n = 100, p = 10, precision_type = "banded",
                                 band = 2, block_size = 9, locations = 50,
                                 durations = 10, change_type = "custom",
                                 changing_vars = 10)
  block_variables <- list("cost"   = c("iid", "cor"),
                          "rho"    = c(0.5, 0.7, 0.9, 0.99))
  many_power_curves(out_file, block_variables, block_data,
                    method_params(), tuning_params(), curve)

  #### LATTICE
  lattice_data <- init_data(n = 100, p = 10, precision_type = "lattice",
                            locations = 50, durations = 10,
                            change_type = "adjacent_lattice")
  lattice_variables <- list("cost"        = c("iid", "cor"),
                            "rho"         = rev(c(0.5, 0.7, 0.9, 0.99)),
                            "proportions" = rev(c(1, 4, 20)/20),
                            "shape"      = rev(c(0, 2, 3, 4, 5)),
                            "precision_est_struct" = c("correct", "banded"),
                            "est_band" = 2)
  many_power_curves(out_file, lattice_variables, lattice_data,
                    method_params(), tuning_params(), curve)

  # #### DIFFERENT BANDS
  # many_power_curves(out_file,
  #                   init_data(n = 100, p = 10, precision_type = "banded",
  #                             locations = 50, durations = 10,
  #                             change_type = "adjacent"),
  #                   rhos = 0.7, props = c(0.1, 0.3, 1), bands = c(1, 3:4),
  #                   tuning = tuning, curve = curve)
}

#' @export
plot_power_curve <- function(file_name, data = init_data(), params = method_params(),
                             tuning = tuning_params(), curve = curve_params(),
                             loc_tol = 10, vars_in_title = NA) {

  all_params <- c(data, params, tuning, curve, list("loc_tol" = loc_tol))
  res <- do.call("rbind",
                 Map(read_power_curve,
                     expand_list(all_params, list("cost" = c("iid", "cor"))),
                     MoreArgs = list(file_name = file_name)))
  title <- make_title(all_params, power_curve_title_parts(vars_in_title))
  ggplot2::ggplot(data = res, ggplot2::aes(x = vartheta, y = power, colour = cost)) +
    ggplot2::geom_line() +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous("Signal strength") +
    ggplot2::scale_y_continuous("Power")
}

#' @export
grid_plot_power <- function(variables = list("rho" = c(-0.3, 0.9),
                                             "proportions" = c(0.1, 0.3, 1)),
                            data = init_data(n = 100, p = 10,
                                             precision_type = "banded",
                                             locations = 50, durations = 10,
                                             change_type = "adjacent"),
                            params = method_params(),
                            tuning = tuning_params(),
                            curve = curve_params(max_dist = 0.1, n_sim = 300),
                            loc_tol = 10) {
  if (length(variables) > 2)
    stop("Maximum two variables at a time.")
  params_list <- split_lists(expand_list(c(data, params), variables),
                             list("data"   = names(data),
                                  "method" = names(params)))
  plots <- Map(plot_power_curve, data = params_list$data,
               params = params_list$method,
               MoreArgs = list("file_name"     = "power_test.csv",
                               "tuning"        = tuning,
                               "curve"         = curve,
                               "loc_tol"       = loc_tol,
                               "vars_in_title" = names(variables)))
  vars_in_title <- names(data)[!names(data) %in% names(variables)]
  dims <- c(length(variables[[2]]), length(variables[[1]]))
  grid_plot(plots, dims, make_title(c(data, params), vars_in_title))
}
