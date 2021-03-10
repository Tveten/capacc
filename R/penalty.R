get_penalty <- function(penalty_regime, n, p, b = 1, a = 2) {
  if (penalty_regime == 'sparse') {
    penalty_obj <- linear_penalty(n, p, b = b, a = a)
    beta <- penalty_obj$beta
    alpha_lin <- penalty_obj$alpha
    alpha_const <- Inf
  } else if (penalty_regime == 'dense') {
    alpha_lin <- 0
    beta <- 0
    alpha_const <- const_penalty(n, p, b = b, a = a)
  } else if (penalty_regime == "combined") {
    alpha_lin <- get_penalty("sparse", n, p, b)[["alpha_lin"]]
    beta <- get_penalty("sparse", n, p, b)[["beta"]]
    alpha_const <- get_penalty("dense", n, p, b)[["alpha_const"]]
  } else if (penalty_regime == 'point') {
    penalty_obj <- linear_penalty(n, p, b = b, a = a)
    beta <- penalty_obj$beta + penalty_obj$alpha
    alpha_lin <- 0
    alpha_const <- Inf
  }
  list("alpha_lin"   = alpha_lin,
       "alpha_const" = alpha_const,
       "beta"        = beta,
       "k_star"      = k_star(alpha_lin, alpha_const, beta))
}

get_penalty_vec <- function(penalty_regime, n, p, b) {
  penalty <- get_penalty(penalty_regime, n, p, b)
  penalty_vec <- 0:p * penalty$beta + penalty$alpha_lin
  if (penalty$k_star <= p)
    penalty_vec[ceiling(penalty$k_star):p + 1] <- penalty$alpha_const
  list("vec" = penalty_vec, "k_star" = penalty$k_star)
}

linear_penalty <- function(n, p, b = 1, a = 2, C = 2) {
  a <- as.numeric(a)
  psi <- a * log(n)
  # k_star <- (p + C * sqrt(p * psi)) / (2 * log(p))
  # list('alpha' = b * C * psi, 'beta' = b * C * log(p), 'k_star' = k_star)
  list('alpha' = b * C * psi, 'beta' = b * C * log(p))
}

const_penalty <- function(n, p, b = 1, a = 2, C = 2) {
  a <- as.numeric(a)
  psi <- a * log(n)
  b * (p + C * psi + C * sqrt(p * psi))
}

k_star <- function(alpha_lin, alpha_const, beta) {
  if (beta > 0) return((alpha_const - alpha_lin) / beta)
  else return(0);
}

iid_penalty <- function(n, p, b = 1) {
    s = 2 * log(n)
    a_vector = qchisq(seq(p-1,0)/p, 1)

    penalty_1 = 2*s + 2*1:p*log(p)
    penalty_2 = rep(p + 2*s + 2*sqrt(p*s),p)

    penalty_3    = 2*(s+log(p)) + 1:p + 2*p*a_vector*dchisq(a_vector, 1) + 2*sqrt((1:p+2*p*a_vector*dchisq(a_vector, 1))*(s+log(p)))
    penalty_3[p] = 2*(s+log(p)) + p + 2*sqrt(p*(s+log(p)))

    beta = b * diff( c(0,pmin(penalty_1,penalty_2,penalty_3)))
    beta
}

iid_point_penalty <- function(n, p, b = 1) {
  point_pen <- get_penalty("point", n, p, b)
  point_pen$beta
}

plot_penalty_func <- function(n = 200, p = 100) {
  pen <- get_penalty_vec("combined", n, p, 1)$vec
  pen_df <- data.frame("J" = 0:p, "pen" = pen)
  ggplot2::ggplot(data = pen_df, ggplot2::aes(x = J, y = pen)) +
    ggplot2::geom_line(color = "red") +
    ggplot2::scale_y_continuous(latex2exp::TeX("P(|\\mathbf{J}|)"),
                                limits = c(0, max(pen_df$pen))) +
    ggplot2::scale_x_continuous(latex2exp::TeX("|\\mathbf{J}|"),
                                breaks = round(seq(0, p, round(p / 5)))) +
    # ggplot2::ggtitle(paste0("n = ", n, ", p = ", p)) +
    ggplot2::theme_classic() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "grey90"),
                   panel.grid.minor = ggplot2::element_line(colour = "grey95"))
}

save_plot_penalty_func <- function(n = 200, p = 100) {
  plot_penalty_func(n, p)
  file_name <- paste0("penalty_func_n", n, "_p", p, ".png")
  ggsave(paste0("./images/", file_name), width = 3, height = 2, units = "in", dpi = 800)
}







