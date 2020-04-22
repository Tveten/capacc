subset_and_check <- function(res, var, var_list, approx = FALSE) {
  if (approx)
    i_call <- paste0("is_equal(", var, ", ", "var_list$", var, ")")
  else
    i_call <- paste0(var, " == ", "var_list$", var)
  res <- res[eval(parse(text = i_call))]
  if (nrow(res) == 0) {
    message(paste0("Results for ", var, "=", eval(parse(text = paste0("var_list$", var))),
                   " does not exist for the this setup."))
    res_exists <- FALSE
  } else
    res_exists <- TRUE
  list(res = res, exists = res_exists)
}

add_iid_costs <- function(res) {
  res <- res[est_band == 0, "cost" := paste0(cost, ".iid")]
  res
}

add_precision_est_struct_to_cost <- function(res) {
  res <- res[cost == "cor" & is.na(est_band), "cost" := paste0(cost, ".", precision_est_struct)]
  res <- res[cost == "cor" & !is.na(est_band), "cost" := paste0(cost, ".", est_band, precision_est_struct)]
  res
}

remove_iid_duplicates <- function(res) {
  iid_est_structs <- unique(res[cost == "iid"]$precision_est_struct)
  if (length(iid_est_structs) > 1)
    res <- res[!(cost == "iid" & precision_est_struct %in%
                   iid_est_structs[2:length(iid_est_structs)] )]
  res
}

read_single_result <- function(file_name, query_params, all_params) {
  res <- fread(paste0("./results/", file_name))
  approx <- vapply(query_params, function(var) is.numeric(all_params[[var]]), logical(1))
  i <- 1
  res_exists <- TRUE
  while (res_exists && i <= length(query_params)) {
    sub_res <- subset_and_check(res, query_params[i], all_params, approx[i])
    res <- sub_res$res
    res_exists <- sub_res$exists
    i <- i + 1
  }
  res <- add_precision_est_struct_to_cost(res)
  res <- remove_iid_duplicates(res)
  res <- add_iid_costs(res)
  res
}

read_power_curve <- function(file_name, all_params) {
  query_params <- c("n", "p", "rho", "precision_type", "block_size", "proportions",
                    "shape", "locations", "durations", "change_type", "band",
                    "cost", "minsl", "maxsl",
                    "alpha", "alpha_tol", "tuning_n_sim",
                    "curve_n_sim", "curve_max_dist",
                    "loc_tol")
  read_single_result(file_name, query_params, all_params)
}

read_cpt_distr <- function(file_name, all_params) {
  query_params <- c("n", "p", "rho", "precision_type", "block_size", "proportions",
                    "vartheta", "shape", "locations", "durations", "change_type", "band",
                    "cost", "minsl", "maxsl", "b",
                    "n_sim")
  read_single_result(file_name, query_params, all_params)
}

already_estimated <- function(file_name, all_params, read_func) {
  if (!file.exists(paste0("./results/", file_name))) {
    message(paste0("File does not exist. Making file ", file_name, " in ./results/"))
    return(FALSE)
  }
  res <- read_func(file_name, all_params)
  if (nrow(res) > 0) {
    message("Results for this setup already exists. Exiting.")
    return(TRUE)
  } else return(FALSE)
}

