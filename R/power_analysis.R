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
    mvcapa_sim <- simulate_detection(data, method, standardise_output = TRUE)
    count_tfp_anom(mvcapa_sim, loc_tol, data)$tp
  })))
  data.table("vartheta" = data$vartheta, "power" = power)
}

est_power_known <- function(data, method, n_sim, cpus = 1) {
  seeds <- sample(1:10^4, n_sim)
  if (method$cost == "cor_exact" && cpus > 1) {
    comp_cluster <- setup_parallel(cpus)
    `%dopar%` <- foreach::`%dopar%`
    power <- mean(foreach::foreach(i = 1:n_sim,
                                   .packages = "capacc",
                                   .combine = "c") %dopar% {
                                     simulate_detection_known(data, method, seed = seeds[i])$S_max > 0
                                   })
    stop_parallel(comp_cluster)
  } else
    power <- mean(unlist(lapply(1:n_sim, function(i) {
      simulate_detection_known(data, method)$S_max > 0
    })))
  data.table("vartheta" = data$vartheta, "power" = power)
}

#' @export
power_curve <- function(out_file, data = init_data(), method = method_params(),
                        tuning = tuning_params(), curve = curve_params(),
                        loc_tol = 10, known = FALSE, seed = NA, cpus = 1) {
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
    if (known) est_power_known(data, method, curve$curve_n_sim, cpus)
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
  all_res <- read_results(out_file)
  if (already_estimated(all_res, all_params, read_func)) return(NULL)

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

many_power_curves <- function(out_file, variables, data = init_data(),
                              method = method_params(), tuning = tuning_params(),
                              curve  = curve_params(max_dist = 0.1, n_sim = 300),
                              known = FALSE, loc_tol = 10, cpus = 1) {
  params <- split_params(
    expand_list(c(data, method), variables),
    list("data"   = names(data),
         "method" = names(method))
  )
  seeds <- get_sim_seeds(variables)
  if (length(seeds) != length(params$data))
    stop("Bug: Length of seeds should be equal to number of parameter settings.")
  if (cpus == 1)
    Map(power_curve,
        data   = params$data,
        method = params$method,
        seed   = seeds,
        MoreArgs = list("out_file"     = out_file,
                        "tuning"        = tuning,
                        "curve"         = curve,
                        "loc_tol"       = loc_tol,
                        "known"         = known))
  else if (cpus > 1) {
    comp_cluster <- setup_parallel(cpus)
    `%dopar%` <- foreach::`%dopar%`
    res <- foreach::foreach(i = 1:length(params$data),
                            .packages = c("anomaly", "capacc")) %dopar% {
                              power_curve(data     = params$data[[i]],
                                          method   = params$method[[i]],
                                          seed     = seeds[i],
                                          out_file = out_file,
                                          tuning   = tuning,
                                          curve    = curve,
                                          loc_tol  = loc_tol,
                                          known    = known)
                            }
    stop_parallel(comp_cluster)
  } else stop("Number of cpus must be >= 1")
  NULL
}

#' @export
plot_power_curve <- function(file_name, variables, data = init_data(),
                             method = method_params(), tuning = tuning_params(),
                             curve = curve_params(), loc_tol = 10, known = FALSE,
                             vars_in_title = NA, dodge = FALSE) {
  upper_xlim <- function(res) {
    ylim <- 0.9
    if (any(res[, all(.SD$power < ylim), by = cost]$V1))
      return(res[, max(vartheta)])
    else
      return(res[power >= ylim, .(vartheta = min(vartheta)), by = cost][, max(vartheta)])
  }

  xbreaks <- function() {
    upper <- max(upper_xlim(res), 0)
    step_length <- 1
    if (upper <= 0.5) step_length <- 0.125
    else if (upper <= 1) step_length <- 0.25
    else if (upper <= 2) step_length <- 0.5
    else if (upper > 6) step_length <- 2
    seq(0, upper, step_length)
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
  all_res <- read_results(file_name)
  res <- do.call("rbind",
                 Map(read_func,
                     params_list,
                     MoreArgs = list(res = all_res)))
  res <- rename_cost(res)
  title <- make_title(all_params, power_curve_title_parts(vars_in_title))
  pl <- ggplot2::ggplot(data = res, ggplot2::aes(x = vartheta, y = power,
                                                 colour = cost))

  if (dodge) {
    width <- max(res$vartheta) / 50
    pl <- pl + ggplot2::geom_line(position = ggplot2::position_dodge(width = width))
  } else
    pl <- pl + ggplot2::geom_line()
  col_name_dt <- cost_names_colours()[name %in% res[, unique(cost)]]
  cols <- col_name_dt$colour
  names(cols) <- col_name_dt$name
  list(
    "plot" = pl +
      ggplot2::ggtitle(latex2exp::TeX(title)) +
      ggplot2::scale_x_continuous(latex2exp::TeX("$\\vartheta$"),
                                  limits = c(0, upper_xlim(res)),
                                  breaks = xbreaks()) +
      ggplot2::scale_y_continuous("Power", breaks = c(0, 0.5, 1), limits = c(0, 1)) +
      ggplot2::scale_colour_manual(name = "Method",
                                   breaks = col_name_dt$name,
                                   labels = unname(latex2exp::TeX(col_name_dt$name)),
                                   values = cols) +
      ggplot2::theme_classic() +
      ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(colour = "grey95")),
    "cost_names" = col_name_dt$name
  )
}

#' @export
grid_plot_power <- function(variables, data = init_data(),
                            method = method_params(), tuning = tuning_params(),
                            curve = curve_params(max_dist = 0.1, n_sim = 300),
                            file_name = "power.csv", known = FALSE, loc_tol = 10,
                            dodge = FALSE, out_file = NULL) {
  get_title <- function() {
    paste0("Power curves for ", make_title(c(data, method), vars_in_title))
  }

  plot_var_names <- c("cost", "precision_est_struct", "est_band")
  plot_variables <- variables[names(variables) %in% plot_var_names]
  grid_variables <- variables[!names(variables) %in% plot_var_names]
  lengths <- unlist(lapply(grid_variables, length))
  if (sum(lengths > 1) > 2)
    stop("Maximum two variables in addition to 'cost' and 'precision_est_struct' at a time.")
  params_list <- split_params(
    expand_list(c(data, method), grid_variables, prune = FALSE),
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
                               "vars_in_title" = names(grid_variables)[lengths > 1],
                               "dodge"         = dodge))
  costs <- extract_nested_element("cost_names", plots)
  plots <- extract_nested_element("plot", plots)

  vars_in_title <- power_curve_title_parts(NA)
  vars_in_title <- vars_in_title[!vars_in_title %in% names(grid_variables)[lengths > 1]]
  n_col <- length(grid_variables[lengths > 1][[1]])
  if (length(grid_variables[lengths > 1]) < 2) n_row <- 1
  else n_row <- length(grid_variables[lengths > 1][[2]])
  dims <- c(n_row, n_col)
  n_costs <- unlist(lapply(costs, length))
  # pp <- grid_plot(plots, dims, latex2exp::TeX(get_title()), which.max(n_costs))
  pp <- grid_plot(plots, dims, NULL, which.max(n_costs))

  if (!is.null(out_file))
    save_grid_plot(pp, dims, out_file, grid_variables[lengths == 1], data)
  else return(pp)
}

#################
#### Unknown single anomaly runs
#################
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

#################
#### MLE runs
#################
#' @export
single_known_anom_MLE_run <- function(p = 10, rho = 0.99, proportions = 1/p,
                                      shape = 0, cpus = 1) {
  out_file <- "power_known_anom_FINAL.csv"
  n <- 100
  data <- init_data(n = n, p = p, precision_type = "banded",
                    band = 2, locations = round(n / 2), durations = 10,
                    change_type = "adjacent", rho = rho,
                    proportions = proportions, shape = shape)
  n_sim <- 300
  tuning <- tuning_params(init_b = c(0.1, 1, 4), n_sim = n_sim)
  curve <- curve_params(max_dist = 0.1, n_sim = n_sim)
  costs <- c("cor", "cor_exact")
  seed <- sample(1:10^4, 1)
  lapply(costs, function(cost) {
    method <- method_params(cost = cost, precision_est_struct = NA)
    power_curve(out_file, data, method, tuning, curve,
                known = TRUE, seed = seed, cpus = cpus)
  })
}

#' @export
plot_single_known_anom_MLE <- function(p = 10, rho = 0.99, proportions = 1/p,
                                       shape = 0) {
  out_file <- "power_known_anom_FINAL.csv"
  n <- 100
  data <- init_data(n = n, p = p, precision_type = "banded",
                    band = 2, locations = round(n / 2), durations = 10,
                    change_type = "adjacent", rho = rho,
                    proportions = proportions, shape = shape)
  method <- method_params(precision_est_struct = NA)
  n_sim <- 300
  tuning <- tuning_params(init_b = c(0.1, 1, 4), n_sim = n_sim)
  curve <- curve_params(max_dist = 0.1, n_sim = n_sim)
  variables <- list("cost" = c("cor", "cor_exact"))
  plot_power_curve(out_file, variables, data, method, tuning, curve, known = TRUE)
}

grid_plot_MLE_sparse_highcor <- function(rho = 0.99, proportions = 1/p,
                                         shape = 0, dodge = FALSE,
                                         out_file = "power_MLE_sparse_highcor") {
  read_file <- "power_known_anom_FINAL.csv"
  n <- 100
  p <- 15
  data <- init_data(n = n, p = p, precision_type = "banded",
                    band = 2, locations = round(n / 2), durations = 10,
                    change_type = "adjacent", rho = rho,
                    proportions = proportions, shape = shape)
  method <- method_params(precision_est_struct = NA)
  n_sim <- 300
  tuning <- tuning_params(init_b = c(0.1, 1, 4), n_sim = n_sim)
  curve <- curve_params(max_dist = 0.1, n_sim = n_sim)
  variables <- list("cost" = c("cor", "cor_exact"),
                    "p"    = c(5, 10, 15),
                    "precision_type" = "banded",
                    "shape"          = shape)
  grid_plot_power(variables, data, method, tuning, curve,
                  file_name = read_file, known = TRUE, dodge = dodge,
                  out_file = out_file)
}

#' @export
known_anom_power_runs_MLE <- function(cpus = 1) {
  curve <- curve_params(max_dist = 0.1, n_sim = 1000)
  out_file <- "power_known_anom_FINAL.csv"
  data <- init_data(n = 100, p = 10, precision_type = "banded",
                    band = 2, locations = 50, durations = 10)
  variables <- list("cost"        = c("cor", "cor_exact"),
                    "precision_est_struct" = c(NA, "correct"),
                    "rho"         = c(0.5, 0.7, 0.9, 0.99),
                    "proportions" = c(1, 3, data$p)/data$p,
                    "shape"       = c(6, 5, 0))
  tuning <- tuning_params(init_b = c(0.1, 1, 2), n_sim = 1000)
  many_power_curves(out_file, variables, data, method_params(),
                    tuning, curve, known = TRUE, cpus = cpus)
}

#' @export
grid_plot_power_MLE <- function(p = 10, rho = "high", shape = 6, dodge = FALSE,
                                out_file = "power_MLE") {
  # Shape = 0 or 5.
  if (p == 8) {
    read_file <- "power_known_anom.csv"
    curve <- curve_params(max_dist = 0.1, n_sim = 300)
    tuning <- tuning_params()
  } else if (p == 10) {
    read_file <- "power_known_anom_FINAL.csv"
    curve <- curve_params(max_dist = 0.1, n_sim = 1000)
    tuning <- tuning_params(init_b = c(0.1, 1, 2), n_sim = 1000)
  }
  data <- init_data(n = 100, p = p, precision_type = "banded",
                    band = 2, locations = 50, durations = 10, shape = shape)
  variables <- list("cost"        = c("cor", "cor_exact"),
                    "precision_est_struct" = c(NA, "correct"))
  if (rho == "low") variables$rho <- c(0.01, 0.2, 0.5)
  else if (rho == "high") variables$rho <- c(0.7, 0.9, 0.99)
  variables$proportions <- c(1/p, 3 / p, 1)
  variables$precision_type <- "banded"
  variables$shape <- shape
  grid_plot_power(variables, data, method_params(), tuning, curve,
                  file_name = read_file, known = TRUE, dodge = dodge,
                  out_file = out_file)
}

#################
#### known anom runs
#################
known_anom_setup <- function(p = 10, precision_type = "banded",
                             shape = c(0, 5, 6),
                             rho = c(0.99, 0.9, 0.7, 0.5, 0.3),
                             proportions = c(1/p, round(sqrt(p)) / p, 1),
                             change_type = "adjacent") {
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
    est_band <- c(1, 2, 4)
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
  data <- init_data(n = n, p = p, precision_type = precision_type[1],
                    band = band, locations = locations, durations = durations,
                    change_type = change_type, shape = shape[1],
                    rho = rho[1], proportions = proportions[1])
  method <- method_params()
  variables <- list("cost"        = c("iid", "cor"),
                    "precision_est_struct" = precision_est_struct,
                    "est_band"    = est_band,
                    "shape"       = shape,
                    "rho"         = rho,
                    "proportions" = proportions,
                    "precision_type" = precision_type)
  tuning <- tuning_params(init_b = c(0.1, 1, 4), n_sim = n_sim)
  list(variables = variables, data = data, method = method,
       tuning = tuning, curve = curve, out_file = out_file)
}

#' @export
known_anom_power_runs <- function(p = 10, precision_type = "banded",
                                  shape = c(0, 5, 6),
                                  rho = c(0.99, 0.9, 0.7, 0.5, 0.3),
                                  proportions = c(1/p, round(sqrt(p)) / p, 1),
                                  change_type = "adjacent",
                                  cpus = 1) {
  setup <- known_anom_setup(p, precision_type, shape, rho, proportions, change_type)
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

  known_anom_power_runs(100, "banded", shape = 0, change_type = "block_scattered",
                        rho = c(0.9, 0.7, 0.5), proportions = 0.1)
  known_anom_power_runs(100, "banded", shape = 8:9, rho = c(0.9, 0.7, 0.5))
  known_anom_power_runs(100, "lattice", shape = 8:9, rho = c(0.9, 0.7, 0.5))
  known_anom_power_runs(100, "lattice", shape = c(0, 5, 6, 8), proportions = c(0.1, 1),
                        change_type = "adjacent_lattice", rho = c(0.9, 0.7, 0.5))
  known_anom_power_runs(100, "global_const", shape = 8, rho = c(0.9, 0.7, 0.5))
}

#' @export
grid_plot_power_known_anom <- function(p = 10, precision_type = "banded",
                                       shape = 6, rho = c(0.5, 0.7, 0.9),
                                       change_type = "adjacent",
                                       proportions = c(1/p, round(sqrt(p)) / p, 1),
                                       dodge = FALSE,
                                       out_file = "power_known_anom") {
  setup <- known_anom_setup(p, precision_type, shape, rho, proportions, change_type)
  grid_plot_power(setup$variables, setup$data, setup$method, setup$tuning,
                  setup$curve, file_name = setup$out_file, known = TRUE,
                  dodge = dodge, out_file = out_file)
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

#################
#### known changepoint runs
#################
known_cpt_setup <- function(p = 10, precision_type = "banded",
                            locations = NULL, shape = c(0, 5, 6),
                            rho = c(0.99, 0.9, 0.7, 0.5, 0.3),
                            proportions = c(1/p, round(sqrt(p)) / p, 1)) {
  if (p <= 50) {
    n <- 100
    n_sim <- 1000
  } else {
    n <- 2 * p
    n_sim <- 500
  }
  if (is.null(locations)) locations <- n - 30

  durations <- n - locations
  band <- 2
  precision_est_struct <- "banded"
  est_band <- c(0, 4)

  curve <- curve_params(max_dist = 0.1, n_sim = n_sim)
  out_file <- "power_known_cpt_FINAL.csv"
  data <- init_data(n = n, p = p, precision_type = precision_type[1],
                    band = band, locations = locations, durations = durations,
                    shape = shape[1], rho = rho[1], proportions = proportions[1])
  method <- method_params()
  variables <- list("cost"        = c("sinspect", "mvlrt"),
                    "precision_est_struct" = precision_est_struct,
                    "est_band"    = est_band,
                    "rho"         = rho,
                    "proportions" = proportions,
                    "shape"       = shape,
                    "precision_type" = precision_type)
  tuning <- tuning_params(init_b = c(0.001, 0.1, 1, 3, 15), n_sim = n_sim)
  list(variables = variables, data = data, method = method,
       tuning = tuning, curve = curve, out_file = out_file)
}

#' @export
known_cpt_power_runs <- function(p = 10, precision_type = "banded",
                                 locations = NULL, shape = c(0, 5, 6),
                                 rho = c(0.99, 0.9, 0.7, 0.5, 0.3),
                                 proportions = c(1/p, round(sqrt(p)) / p, 1),
                                 cpus = 1) {
  setup <- known_cpt_setup(p, precision_type, locations, shape, rho, proportions)
  many_power_curves(setup$out_file, setup$variables, setup$data,
                    setup$method, setup$tuning, setup$curve, known = TRUE,
                    cpus = cpus)
}

#' @export
all_known_cpt_runs10 <- function(locations = NULL, cpus = 1) {
  known_cpt_power_runs(10, "banded", locations, cpus = cpus)
  known_cpt_power_runs(16, "lattice", locations, cpus = cpus)
  known_cpt_power_runs(10, "global_const", locations, cpus = cpus)
}

#' @export
all_known_cpt_runs100 <- function(locations = NULL, cpus = 1) {
  known_cpt_power_runs(100, "banded", locations, cpus = cpus)
  known_cpt_power_runs(100, "lattice", locations, cpus = cpus)
  known_cpt_power_runs(100, "global_const", locations, cpus = cpus)
  known_cpt_power_runs(100, c("banded", "lattice", "global_const"),
                       locations, shape = 8, rho = c(0.5, 0.7, 0.9),
                       proportions = c(0.1, 1), cpus = cpus)
}

#' @export
grid_plot_power_known_cpt <- function(p = 10, precision_type = "banded",
                                      locations = NULL,
                                      shape = 6, rho = c(0.5, 0.7, 0.9),
                                      proportions = c(1/p, round(sqrt(p)) / p, 1),
                                      dodge = FALSE, out_file = "power_known_cpt") {
  # setup <- known_cpt_setup(p, precision_type[1], locations, shape, rho, proportions)
  # setup$variables$precision_type <- precision_type
  setup <- known_cpt_setup(p, precision_type, locations, shape, rho, proportions)
  grid_plot_power(setup$variables, setup$data, setup$method, setup$tuning,
                  setup$curve, file_name = setup$out_file, known = TRUE,
                  dodge = dodge, out_file = out_file)
}

make_all_plots <- function() {
  ### MLE vs. approx plots.
  grid_plot_MLE_sparse_highcor()
  shapes <- c(0, 5, 6)
  Map(grid_plot_power_MLE, shape = shapes)

  ### Single known anomaly plots.
  ps <- rep(c(10, 16, 10, 100, 100, 100), each = 3)
  precision_types <- rep(rep(c("banded", "lattice", "global_const"), each = 3), 2)
  shapes <- rep(c(0, 5, 6), 6)
  Map(grid_plot_power_known_anom,
      p              = ps,
      precision_type = precision_types,
      shape          = shapes)
  ps <- rep(100, 4)
  precision_types <- rep(c("banded", "lattice"), each = 2)
  shapes <- rep(8:9, 2)
  Map(grid_plot_power_known_anom,
      p              = ps,
      precision_type = precision_types,
      shape          = shapes)
  # precision_type vs. rho
  grid_plot_power_known_anom(p = 100, c("banded", "lattice", "global_const"),
                             proportions = 1/100, rho = c(0.5, 0.7, 0.9),
                             shape = 0)
  # Shape vs. precision_type plots.
  rhos <- rep(c(0.5, 0.7, 0.9), 2)
  proportions <- rep(c(0.1, 1), each = 3)
  Map(grid_plot_power_known_anom,
      rho         = rhos,
      proportions = proportions,
      MoreArgs = list(p = 100, shape = c(5, 6, 8, 0),
                      precision_type = c("banded", "lattice", "global_const")))
  # Different change_types.
  grid_plot_power_known_anom(p = 100, precision_type = "banded", shape = 0,
                             change_type = "block_scattered", proportions = 0.1,
                             out_file = "power_known_anom_block_scattered")
  grid_plot_power_known_anom(p = 100, precision_type = "banded", shape = 0,
                             change_type = "adjacent", proportions = 0.1,
                             out_file = "power_known_anom_adjacent")
  Map(grid_plot_power_known_anom, shape = c(0, 5, 6, 8),
      MoreArgs = list(p = 100, precision_type = "lattice",
                      change_type = "adjacent_lattice",
                      out_file = "power_known_anom_adjlat"))


  ### Single changepoint plots.
  ps <- rep(c(10, 16, 10, 100, 100, 100), each = 3)
  precision_types <- rep(rep(c("banded", "lattice", "global_const"), each = 3), 2)
  shapes <- rep(c(0, 5, 6), 6)
  Map(grid_plot_power_known_cpt,
      p              = ps,
      precision_type = precision_types,
      shape          = shapes)
  # precision_type vs. rho
  grid_plot_power_known_cpt(p = 100, c("banded", "lattice", "global_const"),
                            proportions = 1/100, rho = c(0.5, 0.7, 0.9),
                            shape = 0)
  # Shape vs. precision_type plots.
  rhos <- rep(c(0.5, 0.7, 0.9), 2)
  proportions <- rep(c(0.1, 1), each = 3)
  Map(grid_plot_power_known_cpt,
      rho         = rhos,
      proportions = proportions,
      MoreArgs = list(p = 100, shape = c(5, 6, 8, 0),
                      precision_type = c("banded", "lattice", "global_const")))
}
