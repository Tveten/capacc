
small_classify_setup <- function(precision_type = "banded",
                                 point_anoms = FALSE,
                                 shape = c(5, 6, 8),
                                 rho = c(0.9, 0.5),
                                 vartheta = c(2, 1)) {
  p <- 10
  n <- 200
  location <- 50
  duration <- 5
  data <- init_data_mc(p = p, n = n, location = location, duration = duration,
                       precision_type = precision_type,
                       point_anoms = point_anoms, band = 2)
  method <- method_params(b = NA)
  precision_est_struct <- "banded"
  est_band <- c(0, 4)
  variables <- list("cost"        = c("iid", "decor", "cor", "inspect", "gflars", "var_pgl"),
                    "precision_est_struct" = precision_est_struct,
                    "est_band"    = est_band,
                    "rho"         = rho,
                    "vartheta"    = vartheta,
                    "shape"       = shape)
  tuning <- tuning_params(tol = 0.01, init_b = c(0.001, 0.1, 1, 3, 30), n_sim = 500)
  out_file <- "classify_anom_extra.csv"
  list(variables = variables, data = data, method = method,
       tuning = tuning, out_file = out_file)
}

#' @export
small_classify_runs <- function(precision_type = "banded",
                               point_anoms = FALSE, shape = c(5, 6, 8),
                               rho = c(0.9, 0.7, 0.5),
                               vartheta = c(2, 1.5, 1), n_sim = 100, cpus = 1) {
  setup <- small_classify_setup(precision_type, point_anoms, shape, rho, vartheta)
  many_classifications(setup$out_file, setup$variables, setup$data, setup$method,
                       setup$tuning, n_sim, cpus)
}

#' @export
all_small_classify_runs <- function() {
  # small_classify_runs("banded", shape = 6, rho = 0.9, vartheta = 2, n_sim = 10)
  # small_classify_runs("banded", shape = 6, rho = 0.9, vartheta = 2, n_sim = 100)
  small_classify_runs("banded", shape = c(5, 6), rho = 0.9, vartheta = c(1, 2), n_sim = 100)
}

small_classify_results <- function() {
  multi_anom_table(p = 10, vartheta = 2, precision_type = "banded", rho = 0.9,
                   shape = 5, point_anom = FALSE, latex = FALSE, file_name = "classify_anom_extra.csv")
  multi_anom_table(p = 10, vartheta = 1, precision_type = "banded", rho = 0.9,
                   shape = c(5, 6), point_anom = FALSE, latex = FALSE, file_name = "classify_anom_extra.csv")
  multi_anom_table(p = 10, vartheta = 2, precision_type = "banded", rho = 0.9,
                   shape = c(5, 6), point_anom = FALSE, latex = FALSE, file_name = "classify_anom_extra.csv")
  rbind(cbind(data.table(vartheta = rep(1, 2)),
              multi_anom_table(p = 10, vartheta = 1, precision_type = "banded", rho = 0.9,
                               shape = c(5, 6), point_anom = FALSE, latex = FALSE, file_name = "classify_anom_extra.csv")),
        cbind(data.table(vartheta = rep(2, 2)),
              multi_anom_table(p = 10, vartheta = 2, precision_type = "banded", rho = 0.9,
                               shape = c(5, 6), point_anom = FALSE, latex = FALSE, file_name = "classify_anom_extra.csv")))
}





## -----------------------------------------------------------------------------
# CHANGING VARIANCE
## -----------------------------------------------------------------------------
changing_var_setup <- function(precision_type = "banded",
                               point_anoms = FALSE,
                               shape = c(5, 6, 8),
                               rho = c(0.9, 0.5),
                               vartheta = c(2, 1)) {
  p <- 10
  n <- 200
  location <- 50
  duration <- 5
  data <- init_data_mc(p = p, n = n, location = location, duration = duration,
                       precision_type = precision_type,
                       point_anoms = point_anoms, band = 2, n_sd_changes = 10)
  method <- method_params(b = NA)
  precision_est_struct <- "banded"
  est_band <- c(0, 4)
  variables <- list("cost"        = c("iid", "decor", "cor", "inspect", "gflars"),
                    "precision_est_struct" = precision_est_struct,
                    "est_band"    = est_band,
                    "rho"         = rho,
                    "vartheta"    = vartheta,
                    "shape"       = shape)
  tuning <- tuning_params(init_b = c(0.001, 0.1, 1, 3, 30), n_sim = 200)
  out_file <- "classify_anom_sd_change.csv"
  list(variables = variables, data = data, method = method,
       tuning = tuning, out_file = out_file)
}

#' @export
changing_var_runs <- function(precision_type = "banded",
                               point_anoms = FALSE, shape = c(5, 6, 8),
                               rho = c(0.9, 0.7, 0.5),
                               vartheta = c(2, 1.5, 1), n_sim = 100, cpus = 1) {
  setup <- changing_var_setup(precision_type, point_anoms, shape, rho, vartheta)
  many_classifications(setup$out_file, setup$variables, setup$data, setup$method,
                       setup$tuning, n_sim, cpus)
}

#' @export
all_small_classify_runs <- function() {
  changing_var_runs("banded", shape = c(5, 6), rho = 0.9, vartheta = c(1, 2), n_sim = 100)
}

small_classify_results <- function() {
  rbind(cbind(data.table(vartheta = rep(1, 2)),
              multi_anom_table(p = 10, vartheta = 1, precision_type = "banded", rho = 0.9,
                               shape = c(5, 6), point_anom = FALSE, latex = FALSE, file_name = "classify_anom_sd_change.csv")),
        cbind(data.table(vartheta = rep(2, 2)),
              multi_anom_table(p = 10, vartheta = 2, precision_type = "banded", rho = 0.9,
                               shape = c(5, 6), point_anom = FALSE, latex = FALSE, file_name = "classify_anom_sd_change.csv")))
}


adjust_files_sd_changes <- function(file_name) {
  res <- fread(paste0("results/", file_name))
  res[, n_sd_changes := 0]
  cost_ind <- which(names(res) == "cost")
  new_col_ordering <- c(1:(cost_ind - 1), ncol(res), cost_ind:(ncol(res) - 1))
  new_res <- res[, ..new_col_ordering]
  fwrite(new_res, file = paste0("results/", file_name))
}


