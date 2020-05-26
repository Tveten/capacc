
grid_plot <- function(plots, dims, title) {
  figure <- ggpubr::ggarrange(plotlist = plots, nrow = dims[1], ncol = dims[2],
                              common.legend = TRUE, legend = "right")
  title <- ggpubr::text_grob(title, face = "bold", size = 14)
  ggpubr::annotate_figure(figure, top = title)
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
    precision_text <- paste0(params$band, "-banded")
  else if (params$precision_type == "global_const")
    precision_text <- "global cor"
  else precision_text <- params$precision_type
  if (params$block_size < params$p)
    precision_text <- paste0(precision_text, ", m=", params$block_size)

  if (type == "anom") location_text <- paste0("s=", params$locations + 1)
  else if (type == "cpt") location_text <- paste0("cpt=", params$locations)

  alpha_text <- paste0("$\\alpha =", params$alpha, "\\pm ", params$alpha_tol, "$")

  title_parts <- list("precision_type" = precision_text,
                      "rho"            = paste0("$\\rho =", params$rho, "$"),
                      "p"              = paste0("p=", params$p),
                      "n"              = paste0("n=", params$n),
                      "locations"      = location_text,
                      "durations"      = paste0("e=", params$locations + params$durations),
                      "proportions"    = paste0("$|J|=", params$proportions * params$p, "$"),
                      "shape"          = paste0("sh=", params$shape),
                      "b"              = paste0("b=", params$b),
                      "alpha"          = alpha_text,
                      "tuning_n_sim"   = paste0("n_{sim}=", params$tuning_n_sim))
  if (any(is.na(which_parts)))
    which_parts <- names(title_parts)
  paste0(title_parts[names(title_parts) %in% which_parts], collapse = ", ")
}

power_curve_title_parts <- function(vars_in_title) {
  if (any(is.na(vars_in_title)))
    vars_in_title <- c("precision_type", "rho", "p", "n",
                       "locations", "durations", "proportions")
  else vars_in_title
}

cpt_distr_title_parts <- function(vars_in_title) {
  if (any(is.na(vars_in_title)))
    vars_in_title <- c("precision_type", "rho", "p", "n",
                       "locations", "proportions", "shape")
  else vars_in_title
}

penalties_title_parts <- function(vars_in_title) {
  if (any(is.na(vars_in_title)))
    vars_in_title <- c("precision_type", "p", "n", "alpha", "tuning_n_sim")
  else vars_in_title
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

rename_cost <- function(res) {
  res[cost == "cor" & is.na(precision_est_struct), "cost" :=
        "MVCAPA($Q$)"]
  res[cost == "cor" & precision_est_struct == "correct",
      "cost" := "MVCAPA($\\hat{Q}(W)$)"]
  res[cost == "cor" & precision_est_struct == "banded",
      "cost" := paste0("MVCAPA($\\hat{Q}(W(", est_band, "))$)")]
  res[cost == "iid",
      "cost" := paste0("MVCAPA($I$)")]
  res[cost == "inspect" & is.na(precision_est_struct),
      "cost" := "inspect($Q$)"]
  res[cost == "inspect" & precision_est_struct == "correct",
      "cost" := "inspect($\\hat{Q}$)"]
  res[cost == "inspect" & precision_est_struct == "banded",
      "cost" := "inspect($I$)"]
  res[cost == "mvlrt" & is.na(precision_est_struct), "cost" :=
        "MVCPT($Q$)"]
  res[cost == "mvlrt" & precision_est_struct == "correct",
      "cost" := "MVCPT($\\hat{Q}(W)$)"]
  res[cost == "mvlrt" & precision_est_struct == "banded" & est_band != 0,
      "cost" := paste0("MVCPT($\\hat{Q}(W(", est_band, "))$)")]
  res[cost == "mvlrt" & precision_est_struct == "banded" & est_band == 0,
      "cost" := "MVCPT($I$)"]
  res
}
rename_precision_est_struct <- function(res) {
  res <- res[is.na(precision_est_struct), precision_est_struct := "true"]
  res <- res[precision_est_struct == "correct", precision_est_struct := "true adj mat"]
  res <- res[precision_est_struct == "banded",
             precision_est_struct := paste0(est_band, "-", precision_est_struct)]
  res <- res[grepl("iid", cost), precision_est_struct := "true"]
}
