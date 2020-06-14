sim_variable_classification <- function(data, method) {
  u <- rep(0, data$p)
  u[simulate_detection_known(data, method)$J_max] <- 1
  J_true <- get_affected_dims(data$change_type, data$proportions, data$p, data$changing_vars)
  true_labs <- rep(0, data$p)
  true_labs[J_true] <- 1
  out <- c(list(size_J = sum(u)),
           count_binary_tfp(u, true_labs),
           count_rand_tfp(u, true_labs))
  c(out,
    list("precision"      = precision(out),
         "recall"         = recall(out),
         "acc"            = acc(out),
         "bacc"           = bacc(out),
         "rand_precision" = precision(out, "_rand"),
         "rand_recall"    = recall(out, "_rand"),
         "rand_acc"       = acc(out, "_rand"),
         "rand_bacc"      = bacc(out, "_rand"),
         "arand_acc"      = adjusted_rand_index(u, true_labs)))
}

#' @export
sim_subset_est <- function(out_file, data = init_data(), method = method_params(),
                           tuning = tuning_params(), n_sim = 100, seed = NA) {
  add_setup_info <- function(res) {
    which_data <- !(grepl("Sigma", names(data)) | names(data) == "changing_vars")
    res <- cbind(res,
                 as.data.table(data[which_data]),
                 as.data.table(method),
                 as.data.table(tuning[names(tuning) != "init_b"]))
    res$sim_nr <- 1:n_sim
    res$n_sim <- n_sim
    res$seed <-  seed
    res
  }

  message(paste0("Estimating change subset for n=", data$n,
                 ", p=", data$p,
                 ", cost=", method$cost,
                 ", precision=", data$precision_type,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", prop=", data$proportions,
                 ", vartheta=", data$vartheta,
                 ", shape=", data$shape,
                 ", s=", data$locations,
                 ", dur=", data$durations,
                 ", precision_est_struct=", method$precision_est_struct,
                 ", est_band=", method$est_band,
                 "."))
  all_params <- c(data, method, tuning, list("n_sim" = n_sim))
  all_res <- read_results(out_file)
  if (already_estimated(all_res, all_params, read_subset_est)) return(NULL)

  if (!is.na(seed)) set.seed(seed)
  method$b <- get_tuned_penalty(data, method, tuning, TRUE, seed + 2)$b
  sim_with_progress <- dot_every(5, sim_variable_classification)
  res <- as.data.table(t(replicate(n_sim, sim_with_progress(data, method))))
  fwrite(add_setup_info(res), paste0("./results/", out_file), append = TRUE)
}

#' #' @export
many_subset_est <- function(out_file, variables, data = init_data(),
                            method = method_params(), tuning = tuning_params(),
                            n_sim = 100) {
  params <- split_params(
    expand_list(c(data, method), variables),
    list("data"   = names(data),
         "method" = names(method))
  )
  Map(sim_subset_est,
      data = params$data,
      method = params$method,
      seed = get_sim_seeds(variables),
      MoreArgs = list("out_file" = out_file,
                      "tuning"   = tuning,
                      "n_sim"    = n_sim))
  NULL
}

subset_est_setup <- function(p = 10, precision_type = "banded",
                             vartheta = 2, shape = 6,
                             rho = c(0.9, 0.7, 0.5),
                             proportions = c(1/p, round(sqrt(p)) / p)) {
  if (p <= 50) {
    n <- 100
    n_sim <- 1000
    costs <- c("iid", "cor", "cor_exact")
  }
  else {
    n <- 2 * p
    n_sim <- 500
    costs <- c("iid", "cor")
  }

  locations <- round(n / 10)
  durations <- 10
  band <- 2
  precision_est_struct <- "banded"
  est_band <- 4

  out_file <- "subset_est_FINAL.csv"
  data <- init_data(n = n, p = p, precision_type = precision_type, band = band,
                    vartheta = vartheta, locations = locations, durations = durations)
  method <- method_params()
  tuning <- tuning_params(init_b = c(0.1, 1, 3), n_sim = n_sim)
  variables <- list("cost"        = costs,
                    "precision_est_struct" = precision_est_struct,
                    "est_band"    = est_band,
                    "rho"         = rho,
                    "proportions" = proportions,
                    "shape"       = shape,
                    "vartheta"    = vartheta)
  list(variables = variables, data = data, method = method,
       tuning = tuning, out_file = out_file, n_sim = n_sim)
}

#' @export
subset_est_runs <- function(p = 10, precision_type = "banded",
                            vartheta = c(2, 3), shape = 6, rho = c(0.9, 0.7, 0.5),
                            proportions = c(1/p, round(sqrt(p)) / p)) {
  setup <- subset_est_setup(p, precision_type, vartheta, shape, rho, proportions)
  many_subset_est(setup$out_file, setup$variables, setup$data,
                  setup$method, setup$tuning, setup$n_sim)
}

#' @export
all_subset_est_runs <- function() {
  subset_est_runs(10, "banded")
  subset_est_runs(100, "banded")
}

subset_est_table <- function(p = 10, vartheta = 2, rho = 0.9, proportions = 1/p) {
  out_file <- "subset_est_FINAL.csv"
  all_res <- read_results(out_file)
  l <- list(p = p, vartheta = vartheta, rho = rho, prop = proportions)
  res <- all_res[p == l$p & vartheta == l$vartheta & rho == l$rho & proportions == l$prop]
  res
}

precision_recall_hist <- function(p = 10, vartheta = 2, rho = 0.9, proportions = 1/p) {
  out_file <- "subset_est_FINAL.csv"
  all_res <- read_results(out_file)
  l <- list(p = p, vartheta = vartheta, rho = rho, prop = proportions)
  res <- all_res[p == l$p & vartheta == l$vartheta & rho == l$rho & proportions == l$prop]
  print(res[, .(size_J = mean(V1), precision = mean(precision), recall = mean(recall), ari = mean(arand_acc)), by = cost])
  ggplot2::ggplot(data = res, ggplot2::aes(x = V1, fill = cost)) +
    ggplot2::geom_histogram(alpha = 0.8, position = "dodge", binwidth = 0.1)
}
