
#' @export
sim_cpt_est <- function(out_file,
                        data = init_data(n = 200, p = 10, rho = 0.8, band = 2,
                                         locations = 80, durations = 120),
                        method = method_params(cost = "mvlrt"),
                        tuning = tuning_params(),
                        n_sim = 100, seed = NA) {
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

  message(paste0("Estimating changepoints for n=", data$n,
                 ", p=", data$p,
                 ", cost=", method$cost,
                 ", precision=", data$precision_type,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", prop=", data$proportions,
                 ", vartheta=", data$vartheta,
                 ", shape=", data$shape,
                 ", cpt=", data$locations,
                 ", precision_est_struct=", method$precision_est_struct,
                 ", est_band=", method$est_band,
                 "."))
  all_params <- c(data, method, tuning, list("n_sim" = n_sim))
  all_res <- read_results(out_file)
  if (already_estimated(all_res, all_params, read_cpt_est)) return(NULL)

  if (!is.na(seed)) set.seed(seed)
  method$b <- get_tuned_penalty(data, method, tuning, FALSE, seed + 2)$b
  res <- as.data.table(t(replicate(n_sim, simulate_detection(data, method, NULL, FALSE)[c("cpt", "value")])))
  fwrite(add_setup_info(res), paste0("./results/", out_file), append = TRUE)
}

#' #' @export
many_cpt_est <- function(out_file, variables, data = init_data(),
                         method = method_params(), tuning = tuning_params(),
                         n_sim = 100) {
  params <- split_params(
    expand_list(c(data, method), variables),
    list("data"   = names(data),
         "method" = names(method))
  )
  Map(sim_cpt_est,
      data = params$data,
      method = params$method,
      seed = get_sim_seeds(variables),
      MoreArgs = list("out_file" = out_file,
                      "tuning"   = tuning,
                      "n_sim"    = n_sim))
  NULL
}

#' @export
plot_cpt_distr <- function(file_name, variables,
                           data = init_data(),
                           method = method_params(),
                           tuning = tuning_params(),
                           n_sim = 1000, vars_in_title = NA) {
  set_breaks <- function(res) {
    res[, unique(cost)]
  }

  all_params <- c(data, method, tuning, list("n_sim" = n_sim))
  params <- combine_lists(split_params(
    expand_list(all_params, variables),
    list("data"    = names(data),
         "method"  = names(method),
         "tuning"  = names(tuning),
         "n_sim" = "n_sim")
  ))
  all_res <- read_results(out_file)
  res <- do.call("rbind",
                 Map(read_cpt_est,
                     params,
                     MoreArgs = list(res = all_res)))
  res <- rename_cost(res)
  print(mse_cpt_est(res))

  # res <- res[ "all_detect" := .SD[all(value > 0)], by = seed][all_detect]
  title <- make_title(c(data, method), cpt_distr_title_parts(vars_in_title), "cpt")
  ggplot2::ggplot(res, ggplot2::aes(cpt, colour = cost)) +
    ggplot2::stat_bin(alpha = 0.8, ggplot2::aes(y=..density..), geom = "step",
                      position = "identity", binwidth = 1) +
    ggplot2::ggtitle(latex2exp::TeX(title)) +
    ggplot2::scale_x_continuous(name = "Estimated cpt",
                                limits = c(0, res$n[1]),
                                breaks = res$location[1]) +
    ggplot2::scale_colour_discrete(name = "Cost",
                                   breaks = set_breaks(res),
                                   labels = unname(latex2exp::TeX(set_breaks(res))))
}

#' @export
boxplot_cpt_est <- function(file_name, variables,
                            data = init_data(n = 100, p = 10, rho = 0.9, band = 2,
                                             vartheta = 1, proportions = 0.1,
                                             locations = 50, durations = 50),
                            method = method_params(),
                            tuning = tuning_params(),
                            n_sim = 1000, vars_in_title = NA) {

  set_breaks <- function(res) {
    res[, unique(cost)]
  }

  all_params <- c(data, method, tuning, list("n_sim" = n_sim))
  params <- combine_lists(split_params(
    expand_list(all_params, variables),
    list("data"    = names(data),
         "method"  = names(method),
         "tuning"  = names(tuning),
         "n_sim" = "n_sim")
  ))
  all_res <- read_results(out_file)
  res <- do.call("rbind",
                 Map(read_cpt_est,
                     params,
                     MoreArgs = list(res = all_res)))
  res <- rename_cost(res)
  print(mse_cpt_est(res))
  loc_tol <- 10

  title <- make_title(c(data, method), cpt_distr_title_parts(vars_in_title), "cpt")
  ggplot2::ggplot(res, ggplot2::aes(x = cost, y = cpt)) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_x_discrete("Cost",
                              breaks = set_breaks(res),
                              labels = unname(latex2exp::TeX(set_breaks(res)))) +
    ggplot2::scale_y_continuous("Cpt estimate") +
    ggplot2::ggtitle(latex2exp::TeX(title))
}

#' @export
grid_plot_cpt_est <- function(plot_func,
                              variables = list("rho" = c(0.7, 0.9, 0.99),
                                               "proportions" = c(0.1, 1)),
                              data = init_data(n = 100, p = 10, band = 2,
                                               vartheta = 1, locations = 50,
                                               durations = 50),
                              method = method_params(),
                              tuning = tuning_params(),
                              n_sim = 1000,
                              file_name = "cpt_est_FINAL.csv") {
  plot_var_names <- c("cost", "precision_est_struct", "est_band")
  plot_variables <- variables[names(variables) %in% plot_var_names]
  grid_variables <- variables[!names(variables) %in% plot_var_names]
  if (length(grid_variables) > 2)
    stop("Maximum two variables in addition to 'cost' and 'precision_est_struct' at a time.")
  params <- split_params(
    expand_list(c(data, method), grid_variables),
    list("data"   = names(data),
         "method" = names(method))
  )
  plots <- Map(plot_func,
               data = params$data,
               method = params$method,
               MoreArgs = list("file_name"     = file_name,
                               "variables"     = plot_variables,
                               "tuning"        = tuning,
                               "n_sim"         = n_sim,
                               "vars_in_title" = names(variables)))

  vars_in_title <- cpt_distr_title_parts(NA)
  vars_in_title <- vars_in_title[!vars_in_title %in% names(grid_variables)]
  dims <- c(length(grid_variables[[2]]), length(grid_variables[[1]]))
  grid_plot(plots, dims, latex2exp::TeX(make_title(c(data, method), vars_in_title, "cpt")))
}


cpt_est_setup <- function(p = 10, precision_type = "banded",
                          vartheta = 1.5, shape = c(0, 5, 6),
                          rho = c(0.99, 0.9, 0.7, 0.5, 0.3),
                          proportions = c(1/p, round(sqrt(p)) / p, 1)) {
  if (p <= 50) {
    n <- 100
    p_type <- "lowp"
  } else {
    n <- 2 * p
    p_type <- "highp"
  }
  n_sim <- 1000

  locations <- round(7 / 10 * n)
  durations <- n - locations
  band <- 2
  precision_est_struct <- "banded"
  est_band <- c(0, 4)

  out_file <- paste0("cpt_est_", p_type, "_FINAL.csv")
  data <- init_data(n = n, p = p, precision_type = precision_type, band = band,
                    vartheta = vartheta, locations = locations, durations = durations)
  method <- method_params()
  tuning <- tuning_params(init_b = c(0.001, 0.1, 1, 3, 20), n_sim = n_sim)
  variables <- list("cost"        = c("sinspect", "mvlrt"),
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
cpt_est_runs <- function(p = 10, precision_type = "banded",
                         vartheta = c(1, 2, 3), shape = c(0, 5, 6),
                         rho = c(0.99, 0.9, 0.7, 0.5, 0.3),
                         proportions = c(1/p, round(sqrt(p)) / p, 1)) {
  setup <- cpt_est_setup(p, precision_type, vartheta, shape, rho, proportions)
  many_cpt_est(setup$out_file, setup$variables, setup$data,
               setup$method, setup$tuning, setup$n_sim)
}

#' @export
all_cpt_est_runs <- function() {
  cpt_est_runs(10, "banded")
  cpt_est_runs(16, "lattice")
  cpt_est_runs(10, "global_const")
  cpt_est_runs(100, "banded")
  cpt_est_runs(100, "lattice")
  cpt_est_runs(100, "global_const")
}

#' @export
plot_cpt_ests <- function(plot_func = plot_cpt_distr, p = 10,
                          precision_type = "banded", vartheta = 1.5,
                          shape = 6, rho = c(0.5, 0.7, 0.9),
                          proportions = NULL, dodge = FALSE) {
  setup <- cpt_est_setup(p, precision_type, vartheta, shape)
  setup$data$shape <- shape
  if (is.null(proportions)) proportions <- setup$variables$proportions
  variables <- c(setup$variables[c("cost", "precision_est_struct", "est_band")],
                 list("rho" = rho, "proportions" = proportions))
  grid_plot_cpt_est(plot_func, variables, setup$data, setup$method, setup$tuning,
                    setup$n_sim, setup$out_file)
}



####################################################
#### MSE table
####################################################
mse_cpt_est <- function(res, loc_tol = 10) {
  res[,
      "all_tp" := all(value >= 0 &
                        is_in_interval(cpt, c(locations[1] - loc_tol,
                                              locations[1] + loc_tol))),
      by = .(seed, sim_nr)]
  # print(res[, sum(all_tp) / length(unique(cost))])
  out <- cbind(res[, .(n = n[1],
                       p = p[1],
                       shape = shape[1],
                       locations = locations[1],
                       vartheta = vartheta[1],
                       precision_type = precision_type[1],
                       rho  = rho[1],
                       size_J = p[1] * proportions[1],
                       # p_detect = sum(value >= 0) / n_sim[1],
                       power = sum(value >= 0 & is_in_interval(cpt, c(locations[1] - loc_tol, locations[1] + loc_tol))) / n_sim[1],
                       # mad = mean(abs(cpt - locations)),
                       rmse = sqrt(mean((cpt - locations)^2))),
                   by = cost],
               # res[value >= 0 & is_in_interval(cpt, c(locations[1] - loc_tol, locations[1] + loc_tol)),
               #     .(rmse_given_tp = sqrt(mean((cpt - locations)^2)),
               #       n_tp          = .N / length(unique(cost))),
               #     by = cost][, .(rmse_given_tp, n_tp)],
               res[(all_tp),
                   .(rmse_given_alltp = sqrt(mean((cpt - locations)^2)),
                     n_alltp          = .N / length(unique(cost))),
                   by = cost][, .(rmse_given_alltp, n_alltp)])
  out
}

#' @export
write_cpt_mse <- function(p = 10,
                          precision_type = "banded", vartheta = 1.5,
                          shape = 6, rho = c(0.3, 0.5, 0.7, 0.9, 0.99),
                          proportions = NULL, loc_tol = 10) {

  write_mse <- function(res, variables, data = init_data(), method = method_params(),
                        tuning = tuning_params(), n_sim = 1000, loc_tol = 10) {
    write_single_mse <- function(res, variables, data = init_data(),
                                 method = method_params(), tuning = tuning_params(),
                                 n_sim = 1000, loc_tol = 10) {
      all_params <- c(data, method, tuning, list("n_sim" = n_sim))
      params <- combine_lists(split_params(
        expand_list(all_params, variables, prune = FALSE),
        list("data"    = names(data),
             "method"  = names(method),
             "tuning"  = names(tuning),
             "n_sim" = "n_sim")
      ))
      res <- do.call("rbind",
                     Map(read_cpt_est_res,
                         params,
                         MoreArgs = list(res = res)))
      res <- rename_cost(res)
      fwrite(mse_cpt_est(res, loc_tol), "./results/cpt_mse_FINAL.csv", append = TRUE)
    }

    node_var_names <- c("cost", "precision_est_struct", "est_band")
    node_variables <- variables[names(variables) %in% node_var_names]
    grid_variables <- variables[!names(variables) %in% node_var_names]
    params <- split_params(
      expand_list(c(data, method), grid_variables, prune = FALSE),
      list("data"   = names(data),
           "method" = names(method))
    )
    Map(write_single_mse,
        data = params$data,
        method = params$method,
        MoreArgs = list("res"           = res,
                        "variables"     = node_variables,
                        "tuning"        = tuning,
                        "n_sim"         = n_sim,
                        "loc_tol"       = loc_tol))
  }

  setup <- cpt_est_setup(p, precision_type, vartheta, shape)
  setup$data$shape <- shape
  if (is.null(proportions)) proportions <- setup$variables$proportions
  variables <- c(setup$variables[c("cost", "precision_est_struct", "est_band")],
                 list("rho"            = rho,
                      "proportions"    = proportions))
  all_res <- read_results(setup$out_file)
  write_mse(all_res, variables, setup$data, setup$method, setup$tuning,
            setup$n_sim, loc_tol)
}

#' @export
write_all_mses <- function() {
  for (vartheta in 1:3) {
    print(paste0("vartheta = ", vartheta))
    ps <- rep(c(10, 16, 10, 100, 100, 100), each = 3)
    precision_types <- rep(rep(c("banded", "lattice", "global_const"), each = 3), 2)
    shapes <- rep(c(0, 5, 6), 6)
    for (i in 1:length(ps)) {
      print(paste0("p = ", ps[i],
                   ", precision = ", precision_types[i],
                   ", shape = ", shapes[i]))
      write_cpt_mse(p = ps[i], precision_type = precision_types[i],
                    shape = shapes[i], vartheta = vartheta)
      print(gc())
    }
  }
}

mse_row <- function(mse_dt, pt, rh, pr, sh) {
  p <- mse_dt[precision_type == pt]$p[1]
  if (is_equal(pr * p, 1)) search_sh <- 0
  else search_sh <- sh
  mse_dt <- mse_dt[precision_type == pt & rho == rh & size_J == pr * p & shape == search_sh]
  mse_dt[, "shape" := sh]
  unique_costs <- unique(mse_dt$cost)
  if (nrow(mse_dt) > length(unique_costs))
    mse_dt <- mse_dt[1:length(unique_costs)]
  mse_dt <- mse_dt[order(cost)]

  mse_vec <- mse_dt$rmse
  which_min <- which.min(mse_vec)
  which_notmin <- (1:length(mse_vec))[-which_min]
  mse_vec <- round(mse_vec, 2)
  mse_vec <- vapply(mse_vec, function(x) sprintf("%.2f", x), character(1))
  mse_vec[which_min] <- paste0("$\\mathbf{", mse_vec[which_min] , "}$")
  mse_vec[which_notmin] <- paste0("$", mse_vec[which_notmin], "$")
  wide_mse <- matrix(c(mse_dt$precision_type[1],
                       mse_dt$rho[1],
                       mse_dt$size_J[1],
                       mse_vec),
                     ncol = length(mse_vec) + 3)
  if (wide_mse[1, 1] == "banded") wide_mse[1, 1] <- "$\\mathbf{Q}(2)$"
  else if (wide_mse[1, 1] == "lattice") wide_mse[1, 1] <- "$\\mathbf{Q}_\\text{lat}$"
  else if (wide_mse[1, 1] == "global_const") wide_mse[1, 1] <- "$\\mathbf{Q}_\\text{con}$"
  colnames(wide_mse) <- c("$\\mathbf{Q}$", "$\\rho$", "$J$", mse_dt$cost)
  wide_mse[, c(1:3, 5, 4, 7, 6)]
}

rename_precision_type <- function(res) {
  res[precision_type == "banded", "precision_type" := "$Q(2)$"]
  res[precision_type == "lattice", "precision_type" := "$Q_lat$"]
  res[precision_type == "global_const", "precision_type" := "$Q_con$"]
}

latex_mse_table <- function(x, p, vartheta, shape) {
  if (shape == 0) shape_text <- "equal changes"
  else if (shape == 5) shape_text <- "iid changes"
  else if (shape == 6) shape_text <- "cor changes"
  caption <- paste0("RMSE for ",
                    "$p = ", p,
                    "$, $\\vartheta = ", vartheta,
                    "$ and ", shape_text,
                    ". The smallest value is given in bold.")
  label <- paste0("tab:mse_p", p, "_vartheta", vartheta, "_shape", shape)

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
  colnames(x) <- c("$\\bQ$", "$\\rho$", "$J$",
                   "CPT-CC($\\hat{\\bQ}(4)$)",
                   "CPT-CC($\\bI$)",
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

#' @export
cpt_mse_table <- function(p = 10, vartheta = 2, shape = 6,
                          rho = c(0.5, 0.7, 0.9), loc_tol = 10,
                          latex = FALSE) {
  v <- list(p = p, vartheta = vartheta, shape = shape)
  if (p == 10) v$p <- c(10, 16)
  mse_dt <- fread("./results/cpt_mse_FINAL.csv")
  mse_dt <- mse_dt[p %in% v$p & vartheta == v$vartheta]
  proportions <- c(1/p, round(sqrt(p)) / p, 1)
  precision_type <- c("banded", "lattice", "global_const")
  layout <- expand.grid(rh = rho,
                        pt = precision_type,
                        pr = proportions,
                        stringsAsFactors = FALSE)
  if (p == 10) {
    sparse_ind <- layout$pt == "lattice" & layout$pr == 1/10
    medparse_ind <- layout$pt == "lattice" & layout$pr == round(sqrt(10)) / 10
    layout[sparse_ind, "pr"] <- 1 / 16
    layout[medparse_ind, "pr"] <- round(sqrt(16)) / 16
  }
  table <- do.call("rbind", Map(mse_row,
                                pt = layout$pt,
                                rh = layout$rh,
                                pr = layout$pr,
                                MoreArgs = list(mse_dt = mse_dt,
                                                sh = shape)))
  rownames(table) <- NULL
  if (latex) return(latex_mse_table(table, p, vartheta, shape))
  else return(as.data.table(table))
}

most_relevant_tables <- function() {
  cpt_mse_table(100, 2, 5, latex = TRUE)
  cpt_mse_table(100, 2, 6, latex = TRUE)
  cpt_mse_table(100, 2, 0, latex = TRUE)
  cpt_mse_table(100, 3, 5, latex = TRUE)
  cpt_mse_table(100, 3, 6, latex = TRUE)
  cpt_mse_table(100, 3, 0, latex = TRUE)

  # Table in main text:
  cpt_mse_table(100, 3, 6, rho = c(0.5, 0.9), latex = TRUE)
}

