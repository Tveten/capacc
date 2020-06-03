
#' @export
classify_anom <- function(out_file, data = init_data(), method = method_params(),
                          tuning = tuning_params(), n_sim = 100, seed = NA,
                          cpus = 1) {
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

  message(paste0("Simulating multiple anom classification for n=", data$n,
                 ", p=", data$p,
                 ", precision=", data$precision_type,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", shape=", data$shape,
                 ", cost=", method$cost,
                 ", precision_est_struct=", method$precision_est_struct,
                 ", est_band=", method$est_band,
                 "."))
  all_params <- c(format_data(data), method, tuning)
  if (already_estimated(out_file, all_params, read_anom_class)) return(NULL)

  if (!is.na(seed)) set.seed(seed)
  if (is.na(method$b))
    method$b <- get_tuned_penalty(data, method, tuning, FALSE, seed + 2)$b
  sim_classify_with_progress <- dot_every(5, sim_classify)
  res <- as.data.table(t(replicate(n_sim, sim_classify_with_progress())))
  cat("\n")
  print(res)
  # fwrite(add_setup_info(res), paste0("./results/", out_file), append = TRUE)
}


