subset_and_check <- function(res, var, var_list, approx, msg) {
  if (is.na(var_list[[var]]))
    i_call <- paste0("is.na(", var, ")")
  else {
    if (approx)
      i_call <- paste0("is_equal(", var, ", ", "var_list$", var, ")")
    else
      i_call <- paste0(var, " == ", "var_list$", var)
  }
  res <- res[eval(parse(text = i_call))]
  if (nrow(res) == 0) {
    message(paste0(msg, " for ", var, "=", var_list[[var]],
                   " does not exist for this setup."))
    res_exists <- FALSE
  } else
    res_exists <- TRUE
  list(res = res, exists = res_exists)
}

read_single_result <- function(file_name, query_params, all_params, msg) {
  res <- fread(paste0("./results/", file_name), na.strings = "")
  approx <- vapply(query_params, function(var) is.numeric(all_params[[var]]), logical(1))
  i <- 1
  res_exists <- TRUE
  while (res_exists && i <= length(query_params)) {
    sub_res <- subset_and_check(res, query_params[i], all_params, approx[i], msg)
    res <- sub_res$res
    res_exists <- sub_res$exists
    i <- i + 1
  }
  res
}

read_power_curve <- function(file_name, all_params) {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size", "proportions",
                    "shape", "locations", "durations", "change_type",
                    "cost", "minsl", "maxsl", "precision_est_struct", "est_band", "size_mu",
                    "alpha", "alpha_tol", "tuning_n_sim",
                    "curve_n_sim", "curve_max_dist",
                    "loc_tol")
  read_single_result(file_name, query_params, all_params, "A power curve")
}

read_power_curve_known <- function(file_name, all_params) {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size", "proportions",
                    "shape", "locations", "durations", "change_type",
                    "cost", "precision_est_struct", "est_band", "size_mu",
                    "alpha", "alpha_tol", "tuning_n_sim",
                    "curve_n_sim", "curve_max_dist",
                    "loc_tol")
  read_single_result(file_name, query_params, all_params, "A power curve (known anom)")
}


read_cpt_distr <- function(file_name, all_params) {
  query_params <- c("n", "p", "rho", "precision_type", "band", "block_size", "proportions",
                    "vartheta", "shape", "locations", "durations", "change_type",
                    "cost", "minsl", "maxsl", "b", "precision_est_struct", "est_band", "size_mu",
                    "n_sim")
  read_single_result(file_name, query_params, all_params, "A changepoint distribution")
}

read_penalties <- function(file_name, all_params) {
    query_params <- c("n", "p", "rho", "precision_type", "band", "block_size",
                      "cost", "minsl", "maxsl", "precision_est_struct", "est_band", "size_mu",
                      "alpha", "alpha_tol", "tuning_n_sim")
    read_single_result(file_name, query_params, all_params, "A penalty")
}

read_penalties_known <- function(file_name, all_params) {
    query_params <- c("n", "p", "rho", "precision_type", "band", "block_size",
                      "cost", "precision_est_struct", "est_band", "size_mu",
                      "alpha", "alpha_tol", "tuning_n_sim")
    read_single_result(file_name, query_params, all_params, "A penalty (known anom)")
}

already_estimated <- function(file_name, all_params, read_func, out = TRUE) {
  if (!file.exists(paste0("./results/", file_name))) {
    message(paste0("File does not exist. Making file ", file_name, " in ./results/"))
    return(FALSE)
  }
  res <- read_func(file_name, all_params)
  if (nrow(res) > 0) {
    if (out) message("Results for this setup already exists. Exiting.")
    return(TRUE)
  } else return(FALSE)
}

