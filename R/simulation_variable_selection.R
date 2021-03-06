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
                             vartheta = 2, shape = c(5, 6),
                             rho = c(0.9, 0.7, 0.5),
                             proportions = c(1/p, round(sqrt(p)) / p),
                             long_anom = FALSE, alpha = 0.05) {
  if (long_anom) {
    n <- 1000
    locations <- round(n / 2)
    durations <- 200
    if (p <= 50) costs <- c("iid", "cor", "cor_exact")
    else costs <- c("iid", "cor")
  } else {
    if (p <= 50) {
      n <- 100
      costs <- c("iid", "cor", "cor_exact")
    }
    else {
      n <- 2 * p
      costs <- c("iid", "cor")
    }
    locations <- round(n / 10)
    durations <- 10
  }
  n_sim <- 1000

  band <- 2
  precision_est_struct <- c(NA, "banded")
  est_band <- 4

  out_file <- "subset_est_FINAL.csv"
  data <- init_data(n = n, p = p, precision_type = precision_type, band = band,
                    vartheta = vartheta, locations = locations, durations = durations)
  method <- method_params()
  tuning <- tuning_params(init_b = c(0.1, 1, 3), n_sim = n_sim, alpha = alpha,
                          tol = alpha / 2)
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
                            vartheta = c(2, 3, 5), shape = c(5, 6),
                            rho = c(0.9, 0.7, 0.5),
                            proportions = c(1/p, round(sqrt(p)) / p),
                            long_anom = FALSE, alpha = 0.05) {
  setup <- subset_est_setup(p, precision_type, vartheta, shape, rho,
                            proportions, long_anom, alpha)
  many_subset_est(setup$out_file, setup$variables, setup$data,
                  setup$method, setup$tuning, setup$n_sim)
}

#' @export
all_subset_est_runs <- function() {
  # subset_est_runs(10, "banded")
  # subset_est_runs(100, "banded")
  subset_est_runs(100, "banded", proportions = c(2/100, 5/100))
  subset_est_runs(10, "banded", alpha = 0.005)
  subset_est_runs(100, "banded", proportions = c(1/100, 2/100, 5/100, 10/100),
                  alpha = 0.005)
  # subset_est_runs(100, "banded", long_anom = TRUE)
  # subset_est_runs(10, "banded", long_anom = TRUE)
  # subset_est_runs(10, "banded", vartheta = 5)
  # subset_est_runs(100, "banded", vartheta = 5)
}

subset_est_table <- function(p = 10, vartheta = c(2, 5), rho = c(0.5, 0.9),
                             proportions = c(1/p, round(sqrt(p)) / p),
                             shape = 6, durations = 10, alpha = 0.005,
                             latex = FALSE) {
  out_file <- "subset_est_FINAL.csv"
  all_res <- read_results(out_file)
  # if (isTRUE(all.equal(proportions, 1/p))) shape <- 0
  l <- list(p = p, vartheta = vartheta, rho = rho, prop = proportions,
            shape = c(0, shape), durations = durations, alpha = alpha)
  res <- all_res[p %in% l$p & vartheta %in% l$vartheta & rho %in% l$rho & proportions %in% l$prop & shape %in% l$shape & durations == l$durations & alpha == l$alpha]
  # l <- list(p = p, vartheta = vartheta, rho = rho, prop = proportions,
  #           durations = durations)
  # res <- all_res[p %in% l$p & vartheta %in% l$vartheta & rho %in% l$rho & proportions %in% l$prop & durations == l$durations]
  res <- res[!precision_est_struct %in% "correct"]
  res <- rename_cost(res)
  res[, J := round(p * proportions)]
  res <- res[, .SD[, .(size_J    = sprintf("%.2f", mean(size_J)),
                       precision = sprintf("%.2f", mean(precision)),
                       recall    = sprintf("%.2f", mean(recall))),
                 by = cost],
             by = .(J, vartheta, rho)]
  res <- res[order(J), .SD[order(rho)], by = .(J, vartheta)]
  # res <- rbind(res[p == 10 & J == 1 & shape == 0],
  #              res[p == 10 & J == 3 & shape == 6],
  #              res[p == 10 & J == 3 & shape == 5],
  #              res[p == 100 & J == 1 & shape == 0],
  #              res[p == 100 & J == 10 & shape == 6],
  #              res[p == 100 & J == 10 & shape == 5])
  # res[, "shape" := unlist(lapply(shape, rename_shape))]
  # res[shape == "$1$", "shape" := "--"]

  if (latex) return(latex_subset_est_table(res, p, shape, alpha))
  else return(res)
}

latex_subset_est_table <- function(x, p, shape, alpha) {
  if (p == 10) n <- 100
  else if (p == 100) n <- 200
  caption <- paste0("Average precision, recall and $\\hat{J}$ over 1000 repetitions for $p = ", p,
                    "$ and  $n = ", n, "$. ",
                    "Other parameters: $\\bQ = \\bQ(2)$, $s = n / 10$ and $e = s + 10$, ",
                    "$\\bmu_{(", rename_shape(shape), ")}$, $\\alpha = ", alpha, "$.")
  label <- paste0("tab:subset_est_p", p)

  date_line <- paste0('% ', date())
  begin_table <- paste('\\begin{table}[!htb]',
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
  colnames(x) <- c("$J$", "$\\vartheta$", "$\\rho$", "Method",
                   "$\\hat{J}$",
                   "Precision",
                   "Recall")
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

size_J_hist <- function(p = 10, vartheta = 2, shape = 6, rho = 0.9,
                        proportions = 1/p, durations = 10, alpha = 0.05,
                        bw = FALSE) {
  out_file <- "subset_est_FINAL.csv"
  all_res <- read_results(out_file)
  if (isTRUE(all.equal(proportions, 1/p))) shape <- 0
  l <- list(p = p, vartheta = vartheta, rho = rho, prop = proportions,
            shape = shape, durations = durations, alpha = alpha)
  res <- all_res[p %in% l$p & vartheta %in% l$vartheta & rho %in% l$rho & proportions %in% l$prop & shape %in% l$shape & durations == l$durations & alpha == l$alpha]
  res <- res[!precision_est_struct %in% "correct"]
  res <- rename_cost(res)
  # res <- res[size_J > 0]

  graphics_name_dt <- cost_names_graphics()[name %in% res[, unique(cost)]]
  res[, "cost" := factor(cost, levels = graphics_name_dt$name)]
  if (bw) {
    n_costs <- length(unique(res$cost))
    cols <- RColorBrewer::brewer.pal(n_costs + 1, "Greys")[2:(n_costs + 1)]
  } else {
    cols <- graphics_name_dt$colour
    names(cols) <- graphics_name_dt$name
  }
  print(unique(res$size_J))
  if (shape == 0)
    title_text <- paste0("$J = ", l$prop * p, "$")
  else if (shape == 5)
    title_text <- paste0("$J = ", l$prop * p, "$, $\\mu_{(0)}$")
  else if (shape == 6)
    title_text <- paste0("$J = ", l$prop * p, "$, $\\mu_{(\\Sigma)}$")

  ggplot2::ggplot(data = res, ggplot2::aes(x = as.factor(size_J), fill = cost)) +
    ggplot2::geom_bar(position = "dodge") +
    ggplot2::scale_x_discrete(latex2exp::TeX("$\\hat{J}$")) +
    ggplot2::scale_y_continuous("Count per 1000") +
    ggplot2::scale_fill_manual(name = "Method",
                               breaks = graphics_name_dt$name,
                               labels = unname(latex2exp::TeX(graphics_name_dt$name)),
                               values = cols) +
    ggplot2::ggtitle(latex2exp::TeX(title_text)) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey90"),
                   panel.grid.minor = ggplot2::element_line(colour = "grey95"))
}



grid_plot_size_J <- function(p = 10, vartheta = 2, rho = 0.7, durations = 10,
                             proportions = c(1/p, round(sqrt(p)) / p), alpha = 0.05,
                             shape = 6, out_file = "sizeJ", bw = FALSE) {
  # plots <- c(lapply(proportions, size_J_hist, p = p, vartheta = vartheta,
  #                   rho = rho, shape = 6, durations = durations, alpha = alpha),
  #            list(size_J_hist(p, vartheta, 5, rho, proportions[2], durations, alpha)))
  plots <- lapply(proportions, size_J_hist, p = p, vartheta = vartheta,
                  rho = rho, shape = shape, durations = durations, alpha = alpha,
                  bw = bw)
  # print(list(p, vartheta, 5, rho, proportions[2]))
  # plots <- size_J_hist(p, vartheta, 5, rho, proportions[2])
  # print(plots)
  dims <- c(1, length(proportions))
  pp <- grid_plot(plots, dims, "")
  if (!is.null(out_file)) {
    file_name <- paste0("./images/", out_file,
                        "_p", p,
                        "_rho", rho,
                        "_vartheta", vartheta,
                        "_shape", shape,
                        "_dur", durations,
                        "_alpha", alpha,
                        ".png")
    width <- 8
    height <- 3
    show(pp)
    ggplot2::ggsave(file_name, width = width, height = height,
                    units = "in", dpi = 800)
  }
  else return(pp)
}

all_grid_J_plots <- function() {
  for (rho in c(0.7, 0.9)) {
    for (vartheta in c(2, 3, 5)) {
      grid_plot_size_J(p = 10, vartheta = vartheta, rho = rho, alpha = 0.05)
      grid_plot_size_J(p = 10, vartheta = vartheta, rho = rho, alpha = 0.005)
      grid_plot_size_J(p = 100, vartheta = vartheta, proportions = c(0.01, 0.05, 0.1),
                       rho = rho, alpha = 0.05)
      grid_plot_size_J(p = 100, vartheta = vartheta, proportions = c(0.01, 0.05, 0.1),
                       rho = rho, alpha = 0.005)
    }
  }
  grid_plot_size_J(p = 10, vartheta = 2, rho = 0.9, alpha = 0.005,
                   out_file = "sizeJ_bw", bw = TRUE)
}

all_J_tables <- function() {
  subset_est_table(p = 10, vartheta = c(2, 5), proportions = c(0.1, 0.3))
  subset_est_table(p = 10, vartheta = c(2, 5), proportions = c(0.1, 0.3), latex = TRUE)
  subset_est_table(p = 100, vartheta = c(2, 5), proportions = c(0.01, 0.05, 0.1))
  subset_est_table(p = 100, vartheta = c(2, 5), proportions = c(0.01, 0.05, 0.1), latex = TRUE)
}
