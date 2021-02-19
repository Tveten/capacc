
grid_plot <- function(plots, dims, title = NULL, legend_ind = 1) {
  figure <- ggpubr::ggarrange(plotlist = plots, nrow = dims[1], ncol = dims[2],
                              common.legend = TRUE,
                              legend.grob = ggpubr::get_legend(plots[[legend_ind]]),
                              legend = "right")
  if (!is.null(title)) {
    title <- ggpubr::text_grob(title, face = "bold", size = 14)
    ggpubr::annotate_figure(figure, top = title)
  }
  figure
}

save_grid_plot <- function(pp, dims, prefix, const_vars, data) {
  potential_const_vars <- c("precision_type", "rho", "proportions", "shape")
  ids <- c()
  if (any(names(const_vars) == "precision_type")) {
    if (const_vars$precision_type == "banded")
      precision_text <- paste0(data$band, "-banded")
    else
      precision_text <- const_vars$precision_type
    ids[[length(ids) + 1]] <- precision_text
  }
  if (any(names(const_vars) == "shape"))
    ids[[length(ids) + 1]] <- paste0("shape", const_vars$shape)
  if (any(names(const_vars) == "rho")) {
    rho_parts <- strsplit(as.character(const_vars$rho), "[.]")[[1]]
    ids[[length(ids) + 1]] <- paste0("rho", rho_parts[1], rho_parts[2])
  }
  if (any(names(const_vars) == "proportions"))
    ids[[length(ids) + 1]] <- paste0("J", const_vars$proportions * data$p)

  file_name <- paste0("./images/",
                      prefix,
                      "_n", data$n,
                      "_p", data$p,
                      "_", ids[[1]],
                      "_", ids[[2]],
                      "_loc", data$locations,
                      "_dur", data$durations,
                      ".png")
  width <- min(7, 1 + 2 * dims[2])
  height <- min(8, 0.5 + 1.5 * dims[1])
  show(pp)
  ggplot2::ggsave(file_name, width = width, height = height,
                  units = "in", dpi = 800)
}

# test_quote <- function() {
#   rho_expr <- quote(rho)
#   rho <- 0.9
#   rho_part <- bquote(.(rho_expr) == .(rho))
#   precision_part <- ", 2-banded"
#   title <- bquote(.(rho_part) ~ .(precision_part))
#   parts <- list(rho_part, precision_part)
#   title <- lapply(1:length(parts), function(i) {
#     paste0(".(parts[[", i, "]])")
#   })
#   title <- paste0("bquote(", paste(title, collapse = "*"), ")")
#   plot(1, main = eval(parse(text = title)))
# }
#
# title_expr <- function(parts) {
#   title <- lapply(1:length(parts), function(i) {
#     paste0(".(parts[[", i, "]])")
#   })
#   parse(text = paste0("bquote(", paste(title, collapse = "*"), ")"))
# }

make_title <- function(params,
                       which_parts = c("precision_type", "rho", "p", "n",
                                       "locations", "durations", "proportions",
                                       "shape"),
                       type = "anom") {
  if (params$precision_type == "banded")
    precision_text <- paste0("$Q(", params$band, ")$")
  else if (params$precision_type == "global_const")
    precision_text <- "$Q_{con}$"
  else if (params$precision_type == "lattice")
    precision_text <- "$Q_{lat}$"
  else precision_text <- params$precision_type
  if (params$block_size < params$p)
    precision_text <- paste0(precision_text, ", m=", params$block_size)

  if (type == "anom") location_text <- paste0("s=", params$locations + 1)
  else if (type == "cpt") location_text <- paste0("cpt=", params$locations)

  if (params$shape == 0) shape_text <- "$\\mu_{(1)}$"
  else if (params$shape == 5) shape_text <- "$\\mu_{(0)}$"
  else if (params$shape == 6) shape_text <- "$\\mu_{(\\Sigma)}$"
  else if (params$shape == 8) shape_text <- "$\\mu_{(0.8)}$"
  else if (params$shape == 9) shape_text <- "$\\mu_{(0.9)}$"
  else shape_text <- paste0("sh=", params$shape)

  alpha_text <- paste0("$\\alpha =", params$alpha, "\\pm ", params$alpha_tol, "$")

  title_parts <- list("precision_type" = precision_text,
                      "rho"            = paste0("$\\rho =", params$rho, "$"),
                      "p"              = paste0("p=", params$p),
                      "n"              = paste0("n=", params$n),
                      "locations"      = location_text,
                      "durations"      = paste0("e=", params$locations + params$durations),
                      "proportions"    = paste0("$J=", params$proportions * params$p, "$"),
                      "vartheta"       = paste0("$\\vartheta =", params$vartheta, "$"),
                      "shape"          = shape_text,
                      "b"              = paste0("b=", params$b),
                      "alpha"          = alpha_text,
                      "tuning_n_sim"   = paste0("n_{sim} =", params$tuning_n_sim))
  if (any(is.na(which_parts)))
    which_parts <- names(title_parts)
  if (is_equal(params$proportions * params$p, 1))
    which_parts <- which_parts[which_parts != "shape"]
  paste0(title_parts[names(title_parts) %in% which_parts], collapse = ", ")
}

power_curve_title_parts <- function(vars_in_title) {
  if (any(is.na(vars_in_title)))
    return(c("precision_type", "rho", "p", "n",
             "locations", "durations", "proportions", "shape"))
  else return(vars_in_title)
}

cpt_distr_title_parts <- function(vars_in_title) {
  if (any(is.na(vars_in_title)))
    return(c("precision_type", "rho", "p", "n",
             "locations", "vartheta", "proportions", "shape"))
  else return(vars_in_title)
}

penalties_title_parts <- function(vars_in_title) {
  if (any(is.na(vars_in_title)))
    return(c("precision_type", "p", "n", "alpha", "tuning_n_sim"))
  else return(vars_in_title)
}


add_iid_costs <- function(res) {
  res <- res[est_band == 0, "cost" := paste0(cost, ".iid")]
  res
}

add_precision_est_struct_to_cost <- function(res) {
  res[grepl("cor", cost) & precision_est_struct == "correct",
      "precision_est_struct" := "correct_adj"]
  res[grepl("cor", cost) & is.na(precision_est_struct),
      "cost" := paste0(cost, ".", "true")]
  res[grepl("cor", cost) & is.na(est_band) & !is.na(precision_est_struct),
      "cost" := paste0(cost, ".", precision_est_struct)]
  res[grepl("cor", cost) & !is.na(est_band),
      "cost" := paste0(cost, ".", est_band, precision_est_struct)]
  res
}

cost_names_colours <- function() {
  mvcor_cols <- RColorBrewer::brewer.pal(9, "Reds")
  mviid_cols <- RColorBrewer::brewer.pal(9, "Blues")
  mvdecor_col <- "purple"
  mvcor_sparse_col <- "aquamarine3"
  # ml_cols <- c("cyan3", "dodgerblue2")
  ml_cols <- c("chocolate2", "orange2", "gold3")
  inspect_cols <- RColorBrewer::brewer.pal(9, "Greens")
  gflars_col <- "cyan3"
  var_pgl_col <- "magenta1"
  rbind(
    data.table(cost = "cor", precision_est_struct = NA, est_band = NA,
               name = "CAPA-CC($Q$)", colour = mvcor_cols[9]),
    data.table(cost = "cor", precision_est_struct = "correct", est_band = NA,
               name = "CAPA-CC($\\hat{Q}(W^*)$)", colour = mvcor_cols[7]),
    data.table(cost = "cor", precision_est_struct = "banded", est_band = 4,
               name = "CAPA-CC($\\hat{Q}(4)$)", colour = mvcor_cols[6]),
    data.table(cost = "cor", precision_est_struct = "banded", est_band = 2,
               name = "CAPA-CC($\\hat{Q}(2)$)", colour = mvcor_cols[5]),
    data.table(cost = "cor", precision_est_struct = "banded", est_band = 1,
               name = "CAPA-CC($\\hat{Q}(1)$)", colour = mvcor_cols[4]),
    data.table(cost = "cor_sparse", precision_est_struct = "banded", est_band = 4,
               name = "CAPA-CCs($\\hat{Q}(4)$)", colour = mvcor_sparse_col),
    data.table(cost = "iid", precision_est_struct = "banded", est_band = 0,
               name = "MVCAPA", colour = mviid_cols[6]),
    data.table(cost = "decor", precision_est_struct = NA, est_band = NA,
               name = "Whiten + MVCAPA", colour = mvdecor_col),
    data.table(cost = "cor_exact", precision_est_struct = NA, est_band = NA,
               name = "ML($Q$)", colour = ml_cols[3]),
    data.table(cost = "cor_exact", precision_est_struct = "correct", est_band = NA,
               name = "ML($\\hat{Q}(W^*))$)", colour = ml_cols[2]),
    data.table(cost = "cor_exact", precision_est_struct = "banded", est_band = 4,
               name = "ML($\\hat{Q}(4)$)", colour = ml_cols[1]),
    data.table(cost = "inspect", precision_est_struct = NA, est_band = NA,
               name = "inspect($Q$)", colour = inspect_cols[8]),
    data.table(cost = "inspect", precision_est_struct = "correct", est_band = NA,
               name = "inspect($\\hat{Q}$)", colour = inspect_cols[6]),
    data.table(cost = "inspect", precision_est_struct = "banded", est_band = 0,
               name = "inspect($I$)", colour = inspect_cols[4]),
    data.table(cost = "mvlrt", precision_est_struct = NA, est_band = NA,
               name = "CPT-CC($Q$)", colour = mvcor_cols[9]),
    data.table(cost = "mvlrt", precision_est_struct = "correct", est_band = NA,
               name = "CPT-CC($\\hat{Q}(W^*)$)", colour = mvcor_cols[8]),
    data.table(cost = "mvlrt", precision_est_struct = "banded", est_band = 0,
               name = "CPT-CC($I$)", colour = mviid_cols[6]),
    data.table(cost = "mvlrt", precision_est_struct = "banded", est_band = 1,
               name = "CPT-CC($\\hat{Q}(1)$)", colour = mvcor_cols[4]),
    data.table(cost = "mvlrt", precision_est_struct = "banded", est_band = 2,
               name = "CPT-CC($\\hat{Q}(2)$)", colour = mvcor_cols[5]),
    data.table(cost = "mvlrt", precision_est_struct = "banded", est_band = 4,
               name = "CPT-CC($\\hat{Q}(4)$)", colour = mvcor_cols[6]),
    data.table(cost = "gflars", precision_est_struct = NA, est_band = NA,
               name = "Group Fused LARS", colour = gflars_col),
    data.table(cost = "var_pgl", precision_est_struct = NA, est_band = NA,
               name = "VAR DP", colour = var_pgl_col)
  )
}

rename_cost <- function(res) {
  # Anomaly costs
  cost_names <- cost_names_colours()
  calls <- c(
    'cost == "cor" & is.na(precision_est_struct)',
    'cost == "cor" & precision_est_struct == "correct"',
    'cost == "cor" & precision_est_struct == "banded" & est_band == 1',
    'cost == "cor" & precision_est_struct == "banded" & est_band == 2',
    'cost == "cor" & precision_est_struct == "banded" & est_band == 4',
    'cost == "cor_sparse" & precision_est_struct == "banded" & est_band == 4',
    'cost == "iid"',
    'cost == "cor_exact" & is.na(precision_est_struct)',
    'cost == "cor_exact" & precision_est_struct == "correct"',
    'cost == "cor_exact" & precision_est_struct == "banded" & est_band == 4',
    'grepl("inspect", cost) & is.na(precision_est_struct)',
    'grepl("inspect", cost) & precision_est_struct == "correct"',
    'grepl("inspect", cost) & precision_est_struct == "banded"',
    'cost == "mvlrt" & is.na(precision_est_struct)',
    'cost == "mvlrt" & precision_est_struct == "correct"',
    'cost == "mvlrt" & precision_est_struct == "banded" & est_band == 0',
    'cost == "mvlrt" & precision_est_struct == "banded" & est_band == 1',
    'cost == "mvlrt" & precision_est_struct == "banded" & est_band == 2',
    'cost == "mvlrt" & precision_est_struct == "banded" & est_band == 4',
    'cost == "decor"',
    'cost == "gflars"',
    'cost == "var_pgl"'
  )
  for (call in calls) {
    res[eval(parse(text = call)), "cost" := cost_names[eval(parse(text = call)), name]]
  }
  res

  # res[cost == "cor" & precision_est_struct == "correct",
  #     "cost" := "MVCAPA($\\hat{Q}(W^*)$)"]
  # res[cost == "cor" & precision_est_struct == "banded",
  #     "cost" := paste0("MVCAPA($\\hat{Q}(W(", est_band, "))$)")]
  # res[cost == "iid",
  #     "cost" := paste0("MVCAPA($I$)")]
  # res[cost == "cor_exact" & is.na(precision_est_struct), "cost" :=
  #       "ML($Q$)"]
  # res[cost == "cor_exact" & precision_est_struct == "correct",
  #     "cost" := "ML(\\hat{Q}(W^*))$)"]

  # Changepoint costs
  # res[grepl("inspect", cost) & is.na(precision_est_struct),
  #     "cost" := "inspect($Q$)"]
  # res[grepl("inspect", cost) & precision_est_struct == "correct",
  #     "cost" := "inspect($\\hat{Q}$)"]
  # res[grepl("inspect", cost) & precision_est_struct == "banded",
  #     "cost" := "inspect($I$)"]
  # res[cost == "mvlrt" & is.na(precision_est_struct), "cost" :=
  #       "MVCPT($Q$)"]
  # res[cost == "mvlrt" & precision_est_struct == "correct",
  #     "cost" := "MVCPT($\\hat{Q}(W^*)$)"]
  # res[cost == "mvlrt" & precision_est_struct == "banded" & est_band != 0,
  #     "cost" := paste0("MVCPT($\\hat{Q}(W(", est_band, "))$)")]
  # res[cost == "mvlrt" & precision_est_struct == "banded" & est_band == 0,
  #     "cost" := "MVCPT($I$)"]
  # res
}
rename_precision_est_struct <- function(res) {
  res <- res[is.na(precision_est_struct), precision_est_struct := "true"]
  res <- res[precision_est_struct == "correct", precision_est_struct := "true adj mat"]
  res <- res[precision_est_struct == "banded",
             precision_est_struct := paste0(est_band, "-", precision_est_struct)]
  res <- res[grepl("iid", cost), precision_est_struct := "true"]
}

rename_shape <- function(shape) {
  if (shape == 0) shape_text <- "1"
  else if (shape == 5) shape_text <- "0"
  else if (shape == 6) shape_text <- "\\Sigma"
  else if (shape == 8) shape_text <- "0.8"
  else if (shape == 9) shape_text <- "0.9"
  else shape_text <- shape
  return(shape_text)
}
