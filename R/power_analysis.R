#' @export
curve_params <- function(max_dist = 0.2, max_iter = 50, n_sim = 100,
                         init_values = c(0.1, 5)) {
  list("curve_max_dist" = max_dist,
       "curve_max_iter" = max_iter,
       "curve_n_sim"    = n_sim,
       "init_values"    = init_values)
}

est_power <- function(data, method, loc_tol, n_sim) {
  power <- mean(unlist(lapply(1:n_sim, function(i) {
    if (method$cost == "cor_exact" && i %% 10 == 0) cat("=")
    mvcapa_sim <- simulate_mvcapa(data, method, return_anom_only = TRUE)
    count_tfp_anom(mvcapa_sim, loc_tol, data)$tp
  })))
  data.table("vartheta" = data$vartheta, "power" = power)
}

est_power_known <- function(data, method, n_sim, parallelise = TRUE) {
  # if (parallelise) {
  #   seeds <- sample(4000:(10 * (4000 + n_sim)), n_sim)
  #   comp_cluster <- suppressMessages(setup_parallel())
  #   `%dopar%` <- foreach::`%dopar%`
  #   detection <- suppressMessages(foreach::foreach(seed = seeds,
  #                                 .packages = c("mvcapaCor", "anomaly"),
  #                                 .combine = "c") %dopar% {
  #     simulate_mvcapa_known(data, method, seed)$S_max > 0
  #   })
  #   stop_parallel(comp_cluster)
  #   power <- mean(detection)
  # } else
  power <- mean(unlist(lapply(1:n_sim, function(i) {
    simulate_mvcapa_known(data, method)$S_max > 0
  })))
  data.table("vartheta" = data$vartheta, "power" = power)
}

#' @export
power_curve <- function(out_file, data = init_data(), method = method_params(),
                        tuning = tuning_params(), curve = curve_params(),
                        loc_tol = 10, known = FALSE, seed = NA) {
  add_setup_info <- function(res) {
    which_data <- !(grepl("Sigma", names(data)) | names(data) %in% c("changing_vars", "vartheta"))
    res <- cbind(res,
                 as.data.table(data[which_data]),
                 as.data.table(method[names(method) != "b"]),
                 as.data.table(tuning[names(tuning) != "init_b"]),
                 as.data.table(curve[names(curve) != "init_values"]))
    res$mu <- mu_from_vartheta(res$vartheta, data$p, data$proportions)
    res$loc_tol <- loc_tol
    res$seed <-  seed
    res
  }

  est_power_ss <- function(vartheta) {
    data$vartheta <- vartheta
    data <- init_data_(data)
    if (method$cost == "cor_dense_optimal") {
      method$size_mu <- vartheta
      method <- method_params_(method)
    }
    method$b <- get_tuned_penalty(data, method, tuning, known, seed + 2)$b
    cat('.')
    if (known) est_power_known(data, method, curve$curve_n_sim)
    else est_power(data, method, loc_tol, curve$curve_n_sim)
  }

  init_power_est <- function(init_values) {
    min_power <- (tuning$alpha + tuning$alpha_tol) * 1.2
    max_power <- 0.98
    res <- do.call('rbind', Map(est_power_ss, curve$init_values))
    while ((min(res$power) > min_power || max(res$power) < max_power) && max(res$vartheta) <= 10) {
      if (min(res$power) > min_power)
        res <- rbind(res, est_power_ss(min(res$vartheta) / 1.2))
      if (max(res$power) < max_power)
        res <- rbind(res, est_power_ss(max(res$vartheta) * 1.2))
    }
    res
  }

  split_inds <- function(res) {
    # exceeds_max_dist <- adjacent_dist(res[, .(vartheta, power)]) > curve$curve_max_dist
    exceeds_max_dist <- diff(res$power) > curve$curve_max_dist |
                          (diff(res$power) > curve$curve_max_dist / 2 &
                             !is_in_interval(res$power[-1], c(0.1, 0.9)))
    ind <- which(exceeds_max_dist) + 1
    ind <- ind[!(res[ind - 1]$power >= 0.98 & res[ind]$power >= 0.98)]
    ind <- ind[!(res[ind - 1]$power <= 0.02 & res[ind]$power <= 0.02)]
    ind
  }

  if (known && !grepl("known", out_file))
    stop(paste0("out_file ", out_file, " must have 'known' in its name when running with known = TRUE."))

  message(paste0("Estimating power curve for n=", data$n,
                 ", p=", data$p,
                 ", cost=", method$cost,
                 ", precision=", data$precision_type,
                 ", block_size=", data$block_size,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", prop=", data$proportions,
                 ", shape=", data$shape,
                 ", precision_est_struct=", method$precision_est_struct,
                 ", est_band=", method$est_band,
                 "."))
  all_params <- c(data, method, tuning, curve, list("loc_tol" = loc_tol))
  if (known) read_func <- read_power_curve_known
  else       read_func <- read_power_curve
  # print(all_params$tuning_n_sim)
  if (already_estimated(out_file, all_params, read_func)) return(NULL)

  if (!is.na(seed)) set.seed(seed)
  res <- init_power_est(curve$init_values)
  ind <- split_inds(res)
  while (length(ind) > 0 && nrow(res) <= curve$curve_max_iter) {
    new_varthetas <- (res$vartheta[ind - 1] + res$vartheta[ind]) / 2
    res <- rbind(res, do.call('rbind', Map(est_power_ss, new_varthetas)))
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
                              method = method_params(),
                              tuning = tuning_params(),
                              curve  = curve_params(max_dist = 0.1, n_sim = 300),
                              known = FALSE, loc_tol = 10, cpus = 1) {
  params_list <- split_params(
    expand_list(c(data, method), variables),
    list("data"   = names(data),
         "method" = names(method))
  )
  if (cpus == 1)
    Map(power_curve,
        data   = params_list$data,
        method = params_list$method,
        seed   = get_sim_seeds(params_list, variables),
        MoreArgs = list("out_file"     = out_file,
                        "tuning"        = tuning,
                        "curve"         = curve,
                        "loc_tol"       = loc_tol,
                        "known"         = known))
  else if (cpus > 1) {
    parallelMap::parallelStart(mode = "multicore", cpus = cpus, show.info = FALSE)
    parallelMap::parallelLibrary("anomaly", "mvcapaCor")
    parallelMap::parallelMap(
      power_curve,
      data   = params_list$data,
      method = params_list$method,
      seed   = get_sim_seeds(params_list, variables),
      more.args = list("out_file"     = out_file,
                       "tuning"        = tuning,
                       "curve"         = curve,
                       "loc_tol"       = loc_tol,
                       "known"         = known)
    )
    parallelMap::parallelStop()
  } else stop("Number of cpus must be >= 1")
  NULL
}

#' @export
plot_power_curve <- function(file_name, variables, data = init_data(),
                             method = method_params(),
                             tuning = tuning_params(), curve = curve_params(),
                             loc_tol = 10, known = FALSE, vars_in_title = NA,
                             dodge = FALSE) {
  upper_xlim <- function(res) {
    res[power >= 0.98, .(vartheta = min(vartheta)), by = cost][, max(vartheta)]
  }

  all_params <- c(data, method, tuning, curve, list("loc_tol" = loc_tol))
  params_list <- combine_lists(split_params(
    expand_list(all_params, variables),
    list("data"    = names(data),
         "method"  = names(method),
         "tuning"  = names(tuning),
         "curve"   = names(curve),
         "loc_tol" = "loc_tol")
  ))
  if (known) read_func <- read_power_curve_known
  else       read_func <- read_power_curve
  res <- do.call("rbind",
                 Map(read_func,
                     params_list,
                     MoreArgs = list(file_name = file_name)))
  res <- rename_precision_est_struct(res)
  title <- make_title(all_params, power_curve_title_parts(vars_in_title))
  # pl <- ggplot2::ggplot(data = res, ggplot2::aes(x = vartheta, y = power,
  #                                                colour = cost,
  #                                                linetype = precision_est_struct))
  pl <- ggplot2::ggplot(data = res, ggplot2::aes(x = vartheta, y = power,
                                                 colour = precision_est_struct,
                                                 linetype = cost))
  if (dodge) {
    width <- max(res$vartheta) / 50
    pl <- pl + ggplot2::geom_line(position = ggplot2::position_dodge(width = width))
  } else
    pl <- pl + ggplot2::geom_line()
  pl +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous("Signal strength", limits = c(0, upper_xlim(res))) +
    ggplot2::scale_y_continuous("Power") +
    ggplot2::scale_colour_discrete(name = "Cost") +
    ggplot2::scale_linetype_discrete(name = "Precision estimation")
}

#' @export
grid_plot_power <- function(variables = list("rho" = c(-0.3, 0.9),
                                             "proportions" = c(1, 0.3, 0.1)),
                            data = init_data(n = 100, p = 10,
                                             precision_type = "banded",
                                             locations = 50, durations = 10,
                                             change_type = "adjacent"),
                            method = method_params(),
                            tuning = tuning_params(),
                            curve = curve_params(max_dist = 0.1, n_sim = 300),
                            file_name = "power.csv", known = FALSE, loc_tol = 10,
                            dodge = FALSE) {
  plot_var_names <- c("cost", "precision_est_struct", "est_band")
  plot_variables <- variables[names(variables) %in% plot_var_names]
  grid_variables <- variables[!names(variables) %in% plot_var_names]
  if (length(grid_variables) > 2)
    stop("Maximum two variables in addition to 'cost' and 'precision_est_struct' at a time.")
  params_list <- split_params(
    expand_list(c(data, method), grid_variables),
    list("data"   = names(data),
         "method" = names(method))
  )
  plots <- Map(plot_power_curve,
               data = params_list$data,
               method = params_list$method,
               MoreArgs = list("file_name"     = file_name,
                               "variables"     = plot_variables,
                               "tuning"        = tuning,
                               "curve"         = curve,
                               "loc_tol"       = loc_tol,
                               "known"         = known,
                               "vars_in_title" = names(variables),
                               "dodge"         = dodge))
  vars_in_title <- names(data)[!names(data) %in% names(variables)]
  dims <- c(length(grid_variables[[2]]), length(grid_variables[[1]]))
  grid_plot(plots, dims, make_title(c(data, method), vars_in_title))
}

###### Regular runs.
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
additional_power_runs <- function() {
  curve <- curve_params(max_dist = 0.1, n_sim = 300)
  out_file <- "power.csv"

  #### BANDED
  banded_data <- init_data(n = 100, p = 10, precision_type = "banded",
                           band = 2, locations = 50, durations = 10,
                           change_type = "adjacent")
  banded_variables <- list("cost"        = "cor",
                           "rho"         = c(-0.3, 0.01, 0.2, 0.5, 0.7, 0.9, 0.99),
                           "proportions" = c(0.1, 0.3, 1),
                           "shape"       = c(0, 5))
  many_power_curves(out_file, banded_variables, banded_data,
                    method_params(precision_est_struct = NA), tuning_params(), curve)

}

#' @export
low_dim_exact_power_runs <- function() {
  curve <- curve_params(max_dist = 0.2, n_sim = 100)
  out_file <- "power.csv"

  banded_data <- init_data(n = 50, p = 5, precision_type = "banded",
                           band = 2, rho = 0.9, locations = 20, durations = 5,
                           proportions = 1, change_type = "adjacent")
  method <- method_params(maxsl = 10)
  banded_variables <- list("cost"  = c("iid", "cor", "cor_exact"),
                           "precision_est_struct" = c(NA, "correct"),
                           "shape" = c(5, 0))
  many_power_curves(out_file, banded_variables, banded_data,
                    method, tuning_params(), curve, loc_tol = 5)
}

#' @export
grid_plots_power <- function() {
  curve <- curve_params(max_dist = 0.1, n_sim = 300)
  out_file <- "power.csv"
  banded_data <- init_data(n = 100, p = 10, precision_type = "banded",
                           band = 2, locations = 50, durations = 10,
                           change_type = "adjacent", shape = 0)
  grid_plot_power(list("rho" = c(0.01, 0.2, 0.5),
                       "proportions" = c(0.1, 0.3, 1)),
                  banded_data,
                  method_params(),
                  tuning_params(), curve,
                  file_name = out_file)
  grid_plot_power(list("rho" = c(0.7, 0.9, 0.99),
                       "proportions" = c(0.1, 0.3, 1)),
                  banded_data,
                  method_params(),
                  tuning_params(), curve,
                  file_name = out_file)
}

###### Known anomaly runs.
#' @export
known_anom_power_runs_MLE <- function() {
  curve <- curve_params(max_dist = 0.1, n_sim = 300)
  out_file <- "power_known_anom.csv"
  banded_data <- init_data(n = 100, p = 8, precision_type = "banded",
                           band = 2, locations = 50, durations = 10,
                           change_type = "adjacent")
  banded_variables <- list("cost"        = c("iid", "cor", "cor_exact"),
                           "rho"         = c(0.01, 0.2, 0.5, 0.7, 0.9, 0.99),
                           "proportions" = c(1, 3, 8)/8,
                           "shape"       = c(0, 5))
  many_power_curves(out_file, banded_variables, banded_data,
                    method_params(precision_est_struct = NA),
                    tuning_params(), curve, known = TRUE)
  many_power_curves(out_file, banded_variables, banded_data,
                    method_params(precision_est_struct = "correct"),
                    tuning_params(), curve, known = TRUE)
}

#' @export
grid_plot_power_MLE <- function(rho = "high", shape = 0) {
  # Shape = 0 or 5.
  curve <- curve_params(max_dist = 0.1, n_sim = 300)
  out_file <- "power_known_anom.csv"
  tuning <- tuning_params()
  banded_data <- init_data(n = 100, p = 8, precision_type = "banded",
                           band = 2, locations = 50, durations = 10,
                           change_type = "adjacent", shape = shape)
  banded_variables <- list("cost"        = c("iid", "cor", "cor_exact"),
                           "precision_est_struct" = c(NA, "correct"))
  if (rho == "low") banded_variables$rho <- c(0.01, 0.2, 0.5)
  else if (rho == "high") banded_variables$rho <- c(0.7, 0.9, 0.99)
  banded_variables$proportions <- c(1, 3, 8)/8
  grid_plot_power(banded_variables, banded_data, method_params(), tuning,
                  curve, file_name = out_file, known = TRUE, dodge = TRUE)
}

known_anom_setup <- function(p = 10, precision_type = "banded",
                             shape = c(0, 5, 6),
                             rho = c(0.99, 0.9, 0.7, 0.5, 0.3),
                             proportions = c(1/p, round(sqrt(p)) / p, 1)) {
  if (p <= 50) {
    n <- 100
    n_sim <- 1000
  } else {
    n <- 2 * p
    n_sim <- 500
  }
  durations <- 10
  locations <- round(n / 2)
  # proportions <- c(1/p, round(sqrt(p)) / p, 1)

  if (precision_type == "banded") {
    band <- 2
    precision_est_struct <- c(NA, "correct", "banded")
    est_band <- c(1, 4)
  } else if (precision_type == "lattice") {
    band <- NA
    precision_est_struct <- c(NA, "correct", "banded")
    est_band <- c(1, 2, 4)
  } else if (precision_type == "global_const") {
    band <- NA
    precision_est_struct <- "banded"
    est_band <- c(1, 2, 4)
  }

  # rho <- rev(c(0.3, 0.5, 0.7, 0.9, 0.99))
  curve <- curve_params(max_dist = 0.1, n_sim = n_sim)
  out_file <- "power_known_anom_FINAL.csv"
  data <- init_data(n = n, p = p, precision_type = precision_type,
                    band = band, locations = locations, durations = durations)
  method <- method_params()
  variables <- list("cost"        = c("iid", "cor"),
                    "precision_est_struct" = precision_est_struct,
                    "est_band"    = est_band,
                    "rho"         = rho,
                    "proportions" = proportions,
                    "shape"       = shape)
  tuning <- tuning_params(init_b = c(0.1, 1, 4), n_sim = n_sim)
  list(variables = variables, data = data, method = method,
       tuning = tuning, curve = curve, out_file = out_file)
}

#' @export
known_anom_power_runs <- function(p = 10, precision_type = "banded",
                                  shape = c(0, 5, 6),
                                  rho = c(0.99, 0.9, 0.7, 0.5, 0.3),
                                  proportions = c(1/p, round(sqrt(p)) / p, 1),
                                  cpus = 1) {
  setup <- known_anom_setup(p, precision_type, shape, rho, proportions)
  many_power_curves(setup$out_file, setup$variables, setup$data,
                    setup$method, setup$tuning, setup$curve, known = TRUE,
                    cpus = cpus)
}

#' @export
all_known_power_runs10 <- function() {
  known_anom_power_runs(10, "banded")
  known_anom_power_runs(16, "lattice")
  known_anom_power_runs(10, "global_const")
}

#' @export
all_known_power_runs100 <- function() {
  known_anom_power_runs(100, "banded")
  known_anom_power_runs(100, "global_const")
  known_anom_power_runs(100, "lattice")
}

#' @export
grid_plot_power_known_anom <- function(p = 10, precision_type = "banded",
                                       shape = 6, rho = c(0.7, 0.9, 0.99),
                                       proportions = NULL) {
  setup <- known_anom_setup(p, precision_type, shape)
  setup$data$shape <- shape
  if (is.null(proportions)) proportions <- setup$variables$proportions
  variables <- c(setup$variables[c("cost", "precision_est_struct", "est_band")],
                 list("rho" = rho, "proportions" = proportions))
  grid_plot_power(variables, setup$data, setup$method, setup$tuning,
                  setup$curve, file_name = setup$out_file, known = TRUE,
                  dodge = TRUE)
}

#' @export
known_anom_Wishart_power_runs <- function() {
  curve <- curve_params(max_dist = 0.1, n_sim = 1000)
  out_file <- "power_known_anom.csv"
  banded_variables <- list("cost"        = c("iid", "cor"),
                           "rho"         = c(0.7, 0.9, 0.99),
                           "proportions" = c(0.1, 0.3, 1),
                           "shape"       = c(0, 6))
  banded_data <- init_data(n = 100, p = 10, precision_type = "Wishart",
                           band = 2, locations = 50, durations = 10)
  method <- method_params(precision_est_struct = "banded", est_band = 2)
  tuning <- tuning_params(init_b = c(0.1, 1, 4), n_sim = 1000)
  many_power_curves(out_file, banded_variables, banded_data,
                    method, tuning, curve, known = TRUE)
}

#' @export
grid_plot_power_known_anom_Wishart <- function(shape = 6) {
  # Valid shapes: 0 and 6.
  curve <- curve_params(max_dist = 0.1, n_sim = 1000)
  out_file <- "power_known_anom.csv"
  banded_data <- init_data(n = 100, p = 10, precision_type = "Wishart",
                           band = 2, locations = 50, durations = 10,
                           shape = shape)
  banded_variables <- list("cost"        = c("iid", "cor"),
                           "rho"         = c(0.7, 0.9, 0.99),
                           "proportions" = c(0.1, 0.3, 1))
  method <- method_params(precision_est_struct = "banded", est_band = 2)
  tuning <- tuning_params(n_sim = 1000)
  grid_plot_power(banded_variables, banded_data, method, tuning,
                  curve, file_name = out_file, known = TRUE, dodge = TRUE)
}


###### Known dense anomaly runs.
#' @export
known_dense_anom_power_runs <- function() {
  curve <- curve_params(max_dist = 0.1, n_sim = 1000)
  out_file <- "power_known_anom.csv"

  banded_data <- init_data(n = 100, p = 8, precision_type = "banded",
                           band = 2, locations = 50, durations = 10,
                           proportions = 1, shape = 0,
                           change_type = "adjacent")
  banded_variables <- list("cost"        = c("iid_dense", "cor_dense",
                                             "cor_dense_optimal", "cor", "iid"),
                           "precision_est_struct" = c(NA, "correct"),
                           "rho"         = 0.99)
  tuning <- tuning_params(init_b = c(0.0001, 1, 4), n_sim = 1000)
  many_power_curves(out_file, banded_variables, banded_data,
                    method_params(cost = "cor_dense_optimal",
                                  precision_est_struct = NA, size_mu = 1),
                    tuning = tuning, curve = curve, known = TRUE)
}

#' @export
plot_power_known_dense_anom <- function() {
  curve <- curve_params(max_dist = 0.1, n_sim = 300)
  out_file <- "power_known_anom.csv"
  banded_data <- init_data(n = 100, p = 8, precision_type = "banded",
                           band = 2, locations = 50, durations = 10,
                           proportions = 1, shape = 0, rho = 0.99,
                           change_type = "adjacent")
  banded_variables <- list("cost" = c("iid", "cor", "cor_exact", "iid_dense",
                                      "cor_dense", "cor_dense_optimal"),
                           "precision_est_struct" = c(NA, "correct"))
  curve <- curve_params(max_dist = 0.1, n_sim = 1000)
  tuning <- tuning_params(n_sim = 1000)
  plot_power_curve(out_file, banded_variables, banded_data,
                   method_params(cost = "cor_dense_optimal",
                                 precision_est_struct = NA, size_mu = 1),
                   tuning = tuning, curve = curve, known = TRUE, dodge = TRUE)
}

