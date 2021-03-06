subset_and_check <- function(res, var, var_list, msg) {
  if (is.na(var_list[[var]]))
    i_call <- paste0("is.na(", var, ")")
  else if (is.numeric(var_list[[var]]))
    i_call <- paste0("is_equal(", var, ", ", "var_list$", var, ")")
  else
    i_call <- paste0(var, " == ", "var_list$", var)
  res <- res[eval(parse(text = i_call))]
  if (nrow(res) == 0) {
    message(paste0(msg, " for ", var, "=", var_list[[var]],
                   " does not exist for this setup."))
    res_exists <- FALSE
  } else
    res_exists <- TRUE
  list(res = res, exists = res_exists)
}

read_results <- function(file_name) {
  if (!file.exists(paste0("./results/", file_name))) {
    message(paste0("File does not exist. Will make ", file_name, " in ./results/ during run."))
    return(data.table(a = NULL))
  }
  res <- fread(paste0("./results/", file_name))
  res[precision_est_struct == "", precision_est_struct := NA]
  res
}

read_single_result <- function(res, query_params, all_params, msg) {
  if (!is.data.table(res)) stop("res must be a data.table.")
  i <- 1
  if (nrow(res) == 0) res_exists <- FALSE
  else res_exists <- TRUE
  while (res_exists && i <= length(query_params)) {
    sub_res <- subset_and_check(res, query_params[i], all_params, msg)
    res <- sub_res$res
    res_exists <- sub_res$exists
    i <- i + 1
  }
  res
}

read_anom_class <- function(res, all_params, exclude = "none") {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size",
                    "shape", "locations", "durations", "vartheta", "point_locations",
                    "n_sd_changes",
                    "cost", "minsl", "maxsl", "precision_est_struct", "est_band",
                    "alpha", "alpha_tol", "tuning_n_sim",
                    "n_sim")
  query_params <- query_params[!query_params %in% exclude]
  read_single_result(res, query_params, all_params, "Anomaly classification")
}

read_power_curve <- function(res, all_params, exclude = "none") {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size", "proportions",
                    "shape", "locations", "durations", "change_type",
                    "cost", "minsl", "maxsl", "precision_est_struct", "est_band", "size_mu",
                    "alpha", "alpha_tol", "tuning_n_sim",
                    "curve_n_sim", "curve_max_dist",
                    "loc_tol")
  query_params <- query_params[!query_params %in% exclude]
  read_single_result(res, query_params, all_params, "A power curve")
}

read_power_curve_known <- function(res, all_params, exclude = "none") {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size", "proportions",
                    "shape", "locations", "durations", "change_type",
                    "cost", "precision_est_struct", "est_band", "size_mu",
                    "alpha", "alpha_tol", "tuning_n_sim",
                    "curve_n_sim", "curve_max_dist",
                    "loc_tol")
  query_params <- query_params[!query_params %in% exclude]
  res <- read_single_result(res, query_params, all_params, "A power curve (known anom)")
}

read_subset_est <- function(res, all_params, exclude = "none") {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size", "proportions",
                    "vartheta", "shape", "locations", "durations", "change_type",
                    "cost", "minsl", "maxsl", "precision_est_struct", "est_band",
                    "alpha", "alpha_tol", "tuning_n_sim",
                    "n_sim")
  query_params <- query_params[!query_params %in% exclude]
  read_single_result(res, query_params, all_params, "Subset estimates")
}


read_cpt_est <- function(res, all_params, exclude = "none") {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size", "proportions",
                    "vartheta", "shape", "locations", "durations", "change_type",
                    "cost", "minsl", "maxsl", "precision_est_struct", "est_band",
                    "alpha", "alpha_tol", "tuning_n_sim",
                    "n_sim")
  query_params <- query_params[!query_params %in% exclude]
  read_single_result(res, query_params, all_params, "Changepoint estimates")
}

read_cpt_est_res <- function(res, all_params, exclude = "none") {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size", "proportions",
                    "vartheta", "shape", "locations", "durations", "change_type",
                    "cost", "minsl", "maxsl", "precision_est_struct", "est_band",
                    "alpha", "alpha_tol", "tuning_n_sim",
                    "n_sim")
  query_params <- query_params[!query_params %in% exclude]
  read_single_result(res, query_params, all_params, "Changepoint estimates")
}

read_penalties <- function(res, all_params, exclude = "none") {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size",
                    "n_sd_changes",
                    "cost", "minsl", "maxsl", "precision_est_struct", "est_band",
                    "alpha", "alpha_tol", "tuning_n_sim")
  query_params <- query_params[!query_params %in% exclude]
  read_single_result(res, query_params, all_params, "A penalty")
}

read_penalties_known <- function(res, all_params, exclude = "none") {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size",
                    "locations", "durations", "n_sd_changes",
                    "cost", "precision_est_struct", "est_band", "size_mu",
                    "alpha", "alpha_tol", "tuning_n_sim")
  query_params <- query_params[!query_params %in% exclude]
  read_single_result(res, query_params, all_params, "A penalty (known anom)")
}

already_estimated <- function(res, all_params, read_func, out = TRUE) {
  if (!is.data.table(res)) stop("res must be a data.table.")
  res <- read_func(res, all_params)
  if (nrow(res) > 0) {
    if (out) message("Results for this setup already exists. Exiting.")
    return(TRUE)
  } else return(FALSE)
}

