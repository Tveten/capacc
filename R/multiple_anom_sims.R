
#' @export
classify_anom <- function(out_file, data = init_data(), method = method_params(),
                          tuning = tuning_params(), n_sim = 100, seed = NA) {
  format_data <- function(data) {
    data$locations <- data$locations[1]
    data$durations <- min(data$durations)
    data$vartheta <- data$vartheta[1]
    data$mu <- data$mu[1]
    data$proportions <- data$proportions[1]
    data$change_type <- data$change_type[1]
    data$point_locations <- data$point_locations[1]
    data$point_proportions <- data$point_proportions[1]
    data$point_mu <- data$point_mu[1]
    which_data <- !(grepl("Sigma", names(data)) | names(data) == c("changing_vars"))
    data[which_data]
  }

  add_setup_info <- function(res) {
    res <- cbind(res,
                 as.data.table(format_data(data)),
                 as.data.table(method[!names(method) %in% c("size_mu")]),
                 as.data.table(tuning[names(tuning) != "init_b"]))
    res$n_sim <- n_sim
    res$seed <-  seed
    res
  }

  sim_classify <- function() {
    sim <- simulate_detection(data, method)
    est_labs <- label_anom_est(sim, data)
    true_labs <- true_labels(data)
    out <- c(list("K_hat" = count_collective_anomalies(sim)),
             count_binary_tfp(est_labs, true_labs),
             count_rand_tfp(est_labs, true_labs))
    out <- c(out,
             list("precision"      = precision(out),
                  "recall"         = recall(out),
                  "acc"            = acc(out),
                  "bacc"           = bacc(out),
                  "rand_precision" = precision(out, "_rand"),
                  "rand_recall"    = recall(out, "_rand"),
                  "rand_acc"       = acc(out, "_rand"),
                  "rand_bacc"      = bacc(out, "_rand"),
                  "arand_acc"      = adjusted_rand_index(est_labs, true_labs)))
  }

  message(paste0("Multiple anomaly classification for n=", data$n,
                 ", p=", data$p,
                 ", precision=", data$precision_type,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", vartheta=", data$vartheta[1],
                 ", shape=", data$shape,
                 ", cost=", method$cost,
                 ", precision_est_struct=", method$precision_est_struct,
                 ", est_band=", method$est_band,
                 "."))
  all_params <- c(format_data(data), method, tuning)
  all_res <- read_results(out_file)
  if (already_estimated(all_res, all_params, read_anom_class)) return(NULL)

  if (!is.na(seed)) set.seed(seed)
  if (is.na(method$b))
    method$b <- get_tuned_penalty(data, method, tuning, FALSE, seed + 2)$b
  sim_classify_with_progress <- dot_every(5, sim_classify)
  res <- as.data.table(t(replicate(n_sim, sim_classify_with_progress())))
  cat("\n")
  print(res)
  fwrite(add_setup_info(res), paste0("./results/", out_file), append = TRUE)
}

many_classifications <- function(out_file, variables, data = init_data(),
                                 method = method_params(), tuning = tuning_params(),
                                 n_sim = 100, cpus = 1) {
  params <- split_params(
    expand_list(c(data, method), variables),
    list("data"   = names(data),
         "method" = names(method))
  )
  seeds <- get_sim_seeds(variables)
  if (length(seeds) != length(params$data))
    stop("Bug: Length of seeds should be equal to number of parameter settings.")
  NULL
  if (cpus == 1)
    Map(classify_anom,
        data   = params$data,
        method = params$method,
        seed   = seeds,
        MoreArgs = list("out_file"      = out_file,
                        "tuning"        = tuning,
                        "n_sim"         = n_sim))
  else if (cpus > 1) {
    comp_cluster <- setup_parallel(cpus)
    `%dopar%` <- foreach::`%dopar%`
    res <- foreach::foreach(i = 1:length(params$data),
                            .packages = c("anomaly", "mvcapaCor")) %dopar% {
      classify_anom(out_file, params$data[[i]], params$method[[i]],
                    tuning, n_sim, seeds[i])
    }
    stop_parallel(comp_cluster)
  } else stop("Number of cpus must be >= 1")
  NULL
}

multiple_anom_setup <- function(p = 10, precision_type = "banded",
                                location = 300, duration = 10,
                                point_anoms = FALSE, shape = c(5, 6, 8),
                                rho = c(0.9, 0.7, 0.5),
                                vartheta = c(2, 1.5, 1)) {
  data <- init_data_mc(p = p, location = location, duration = duration,
                       precision_type = precision_type,
                       point_anoms = point_anoms, band = 2)
  method <- method_params(b = NA)
  precision_est_struct <- "banded"
  est_band <- c(0, 4)
  variables <- list("cost"        = c("iid", "cor", "inspect"),
                    "precision_est_struct" = precision_est_struct,
                    "est_band"    = est_band,
                    "rho"         = rho,
                    "vartheta"    = vartheta,
                    "shape"       = shape)
  tuning <- tuning_params(init_b = c(0.1, 1, 4), n_sim = 200)
  out_file <- "multiple_anom_FINAL.csv"
  list(variables = variables, data = data, method = method,
       tuning = tuning, out_file = out_file)
}

#' @export
multiple_anom_runs <- function(p = 10, precision_type = "banded",
                               location = 300, duration = 10,
                               point_anoms = FALSE, shape = c(5, 6, 8),
                               rho = c(0.9, 0.7, 0.5),
                               vartheta = c(2, 1.5, 1), n_sim = 100, cpus = 1) {
  setup <- multiple_anom_setup(p, precision_type, location, duration,
                               point_anoms, shape, rho, vartheta)
  many_classifications(setup$out_file, setup$variables, setup$data, setup$method,
                       setup$tuning, n_sim, cpus)
}

all_multiple_anom_runs10 <- function(cpus = 1) {
  multiple_anom_runs(10, "banded", cpus = cpus)
  multiple_anom_runs(10, "banded", point_anoms = TRUE, cpus = cpus)
  multiple_anom_runs(16, "lattice", cpus = cpus)
  multiple_anom_runs(16, "lattice", point_anoms = TRUE, cpus = cpus)
  multiple_anom_runs(10, "global_const", cpus = cpus)
  multiple_anom_runs(10, "global_const", point_anoms = TRUE, cpus = cpus)
}

#' @export
all_multiple_anom_runs100 <- function(cpus = 1) {
  multiple_anom_runs(100, "lattice", cpus = cpus)
  multiple_anom_runs(100, "lattice", point_anoms = TRUE, cpus = cpus)
  multiple_anom_runs(100, "banded", cpus = cpus)
  multiple_anom_runs(100, "banded", point_anoms = TRUE, cpus = cpus)
  multiple_anom_runs(100, "global_const", cpus = cpus)
  multiple_anom_runs(100, "global_const", point_anoms = TRUE, cpus = cpus)
}

