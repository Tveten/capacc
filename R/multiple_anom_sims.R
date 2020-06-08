
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

correct_res <- function() {
  file_name <- "./results/multiple_anom_FINAL.csv"
  res <- fread(file_name)
  print(dim(res))
  print(res[, length(unique(seed)), by = .(p, precision_type, shape, rho, vartheta, point_locations)]$V1)
  print(res[, .N, by = .(p, precision_type, shape, rho, vartheta, point_locations)]$N)
  res <- res[, seed := .SD[cost == "iid", seed][1],
             by = .(p, precision_type, shape, rho, vartheta, point_locations)]
  print(dim(res))
  print(res[, length(unique(seed)), by = .(p, precision_type, shape, rho, vartheta, point_locations)]$V1)
  print(res[, .N, by = .(p, precision_type, shape, rho, vartheta, point_locations)]$N)
}

#' @export
classify_anom <- function(out_file, data = init_data(), method = method_params(),
                          tuning = tuning_params(), n_sim = 100, seed = NA) {
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

  get_seed <- function() {
    all_params <- c(format_data(data), method, tuning, list(n_sim = n_sim))
    all_res <- read_results(out_file)
    exclude <- c("cost", "precision_est_struct", "est_band")
    prev_res <- read_anom_class(all_res, all_params, exclude = exclude)
    prev_seed <- unique(prev_res$seed)
    if (length(prev_seed) > 0) {
      message(paste0("Overriding input seed with the same seed used for other costs: ",
                     paste(prev_seed, collapse = ", ")))
      return(prev_seed[1])
    } else return(seed)
  }

  message(paste0("Multiple anomaly classification for n=", data$n,
                 ", p=", data$p,
                 ", precision=", data$precision_type,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", vartheta=", data$vartheta[1],
                 ", shape=", data$shape,
                 ", point_anom=", data$point_locations[1],
                 ", cost=", method$cost,
                 ", precision_est_struct=", method$precision_est_struct,
                 ", est_band=", method$est_band,
                 "."))
  all_params <- c(format_data(data), method, tuning, list(n_sim = n_sim))
  all_res <- read_results(out_file)
  if (already_estimated(all_res, all_params, read_anom_class)) return(NULL)

  seed <- get_seed()
  if (!is.na(seed)) set.seed(seed)
  if (is.na(method$b) && !is.null(method$b))
    method$b <- get_tuned_penalty(data, method, tuning, FALSE, seed + 2)$b
  sim_classify_with_progress <- dot_every(5, sim_classify)
  res <- as.data.table(t(replicate(n_sim, sim_classify_with_progress())))
  cat("\n")
  print(paste0("Results for cost=", method$cost,
               ", precision_est_struct=", method$precision_est_struct,
               ", est_band=", method$est_band))
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
  tuning <- tuning_params(init_b = c(0.001, 0.1, 1, 3, 30), n_sim = 200)
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

#' @export
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

multi_anom_row <- function(pm_dt, perf_metric, pt, rh, sh, pa) {
  pm_dt <- pm_dt[precision_type == pt & rho == rh & shape == sh & point_anom == pa]
  seeds <- unique(pm_dt$seed)
  if (length(seeds) > 1) {
    warning("Number of seeds is greater than 1. Picking the first one.")
    pm_dt <- pm_dt[seed == seeds[1]]
  }
  pm_dt <- pm_dt[, .(pt = pt,
                     rh = rh,
                     sh = sh,
                     pa = pa,
                     pm = mean(eval(parse(text = perf_metric)))),
                 by = cost]
  pm_dt <- pm_dt[order(cost)]

  pm_vec <- pm_dt$pm
  which_max <- which.max(pm_vec)
  which_notmax <- (1:length(pm_vec))[-which_max]
  pm_vec <- round(pm_vec, 2)
  pm_vec <- vapply(pm_vec, function(x) sprintf("%.2f", x), character(1))
  pm_vec[which_max] <- paste0("$\\mathbf{", pm_vec[which_max] , "}$")
  pm_vec[which_notmax] <- paste0("$", pm_vec[which_notmax], "$")
  wide_pm <- matrix(c(pt, rh, sh, pa, pm_vec), ncol = length(pm_vec) + 4)
  colnames(wide_pm) <- c("$\\mathbf{Q}$", "$\\rho$", "shape",
                          "pt_anoms", pm_dt$cost)

  if (wide_pm[1, 1] == "banded") wide_pm[1, 1] <- "$\\mathbf{Q}(2)$"
  else if (wide_pm[1, 1] == "lattice") wide_pm[1, 1] <- "$\\mathbf{Q}_\\text{lat}$"
  else if (wide_pm[1, 1] == "global_const") wide_pm[1, 1] <- "$\\mathbf{Q}_\\text{con}$"

  wide_pm[, c(1:4, 6, 5, 8, 7)]
}

# Varying cost, rho, vartheta shape and precision_type in sims.
# Table:
# given p, vartheta or not.
# precision_type rho, shape, point_anom, cost1 ... cost4
multi_anom_table <- function(p = 10, vartheta = 2, perf_metric = "arand_acc",
                             precision_type = c("banded", "lattice", "global_const"),
                             rho = c(0.5, 0.7, 0.9),
                             shape = c(5, 6, 8),
                             point_anom = c(FALSE, TRUE),
                             latex = FALSE) {
  v <- list(p = p, vartheta = vartheta)
  if (p == 10) v$p <- c(10, 16)
  pm_dt <- fread("./results/multiple_anom_FINAL.csv")
  pm_dt[, "point_anom" := !is.na(point_locations)]
  pm_dt <- pm_dt[p %in% v$p & vartheta == v$vartheta]
  pm_dt <- rename_cost(pm_dt)
  layout <- expand.grid(rh = rho,
                        sh = shape,
                        pa = point_anom,
                        pt = precision_type,
                        stringsAsFactors = FALSE)
  table <- do.call("rbind", Map(multi_anom_row,
                                pt = layout$pt,
                                rh = layout$rh,
                                sh = layout$sh,
                                pa = layout$pa,
                                MoreArgs = list(pm_dt = pm_dt,
                                                perf_metric = perf_metric)))
  rownames(table) <- NULL
  table <- as.data.table(table)
  if (latex) return(latex_ari_table(table, p, vartheta))
  else return(table)
}


latex_ari_table <- function(x, p, vartheta) {
  rename_shape <- function(shape) {
    # if (shape == 0) shape_text <- "$\\mu_{(1)}$"
    # else if (shape == 5) shape_text <- "$\\mu_{(0)}$"
    # else if (shape == 6) shape_text <- "$\\mu_{(\\Sigma)}$"
    # else if (shape == 8) shape_text <- "$\\mu_{(0.8)}$"
    # else if (shape == 9) shape_text <- "$\\mu_{(0.9)}$"
    # else shape_text <- paste0("sh=", shape)
    if (shape == 0) shape_text <- "$1$"
    else if (shape == 5) shape_text <- "$0$"
    else if (shape == 6) shape_text <- "$\\Sigma$"
    else if (shape == 8) shape_text <- "$0.8$"
    else if (shape == 9) shape_text <- "$0.9$"
    else shape_text <- shape
    return(shape_text)
  }

  x[, shape := unlist(lapply(shape, rename_shape))]
  x[pt_anoms == TRUE, pt_anoms := "\\checkmark"]
  x[pt_anoms == FALSE, pt_anoms := "--"]

  caption <- paste0("ARI for ",
                    "$p = ", p,
                    "$, $\\vartheta = ", vartheta,
                    "$. The largest value is given in bold.")
  label <- paste0("tab:ari_p", p, "_vartheta", vartheta)

  date_line <- paste0('% ', date())
  begin_table <- paste('\\begin{table}[htb]',
                       paste0('\\caption{', caption, '}'),
                       paste0('\\label{', label ,'}'),
                       '\\centering',
                       paste0('\\begin{tabular}{',
                              paste(rep('c', ncol(x)), collapse = ""),
                              '}'),
                       '\\toprule', sep = ' \n')
  end_table <- paste('\\bottomrule',
                     '\\end{tabular}',
                     '\\end{table}', sep = ' \n')

  mid_sep <- ' \n\\midrule'
  colnames(x) <- c("$\\bQ$", "$\\rho$", "$\\mu_{(\\cdot)}$", "Pt. anoms",
                   "MVCAPA($\\hat{\\bQ}(4)$)",
                   "MVCAPA($\\bI$)",
                   "inspect($\\hat{\\bQ}$)",
                   "inspect($\\bI$)")
  heading <- paste(colnames(x), collapse = " & ")

  latex_table <- paste(date_line, begin_table, heading, sep = ' \n')
  latex_table <- paste0(latex_table, " \\\\", mid_sep)
  for (i in 1:nrow(x)) {
    table_line <- paste(x[i, ], collapse = " & ")
    table_line <- paste0(table_line, " \\\\")
    latex_table <- paste(latex_table, table_line, sep = " \n")
  }
  latex_table <- paste0(latex_table, ' \n', end_table)
  cat(latex_table)
}

tables_to_show <- function() {
  multi_anom_table(p = 10, vartheta = 1.5)
  multi_anom_table(p = 10, vartheta = 1.5, shape = c(5, 8), rho = c(0.5, 0.9), latex = TRUE)
  multi_anom_table(p = 100, vartheta = 1.5)
  multi_anom_table(p = 100, vartheta = 1.5, shape = c(5, 8), rho = c(0.5, 0.9))
}




