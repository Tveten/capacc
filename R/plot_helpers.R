
grid_plot <- function(plots, dims, title) {
  figure <- ggpubr::ggarrange(plotlist = plots, nrow = dims[1], ncol = dims[2],
                              common.legend = TRUE, legend = "right")
  title <- ggpubr::text_grob(title, face = "bold", size = 14)
  ggpubr::annotate_figure(figure, top = title)
}

make_title <- function(params,
                       which_parts = c("cor_mat_type", "rho", "p", "n",
                                       "locations", "durations", "proportions",
                                       "shape"),
                       type = "anom") {
  if (params$precision_type == "lattice")
    precision_text <- "lattice"
  else if (params$precision_type == "banded")
    precision_text <- paste0(params$band, "-banded")
  if (params$block_size < params$p)
    precision_text <- paste0(precision_text, ", m=", params$block_size)

  if (type == "anom") location_text <- paste0("s=", params$locations + 1)
  else if (type == "cpt") location_text <- paste0("cpt=", params$locations)

  title_parts <- list("precision_type" = precision_text,
                      "rho"            = paste0("rho=", params$rho),
                      "p"              = paste0("p=", params$p),
                      "n"              = paste0("n=", params$n),
                      "locations"      = location_text,
                      "durations"      = paste0("e=", params$locations + params$durations),
                      "proportions"    = paste0("pr=", round(params$proportions, 2)),
                      "shape"          = paste0("sh=", params$shape),
                      "b"              = paste0("b=", params$b))
  if (any(is.na(which_parts)))
    which_parts <- names(title_parts)
  paste0(title_parts[names(title_parts) %in% which_parts], collapse = ", ")
}

power_curve_title_parts <- function(vars_in_title) {
  if (any(is.na(vars_in_title)))
    vars_in_title <- c("cor_mat_type", "rho", "p", "n",
                       "locations", "durations", "proportions")
  else vars_in_title
}

cpt_distr_title_parts <- function(vars_in_title) {
  if (any(is.na(vars_in_title)))
    vars_in_title <- c("cor_mat_type", "rho", "p", "n",
                       "locations", "proportions", "b")
  else vars_in_title
}
