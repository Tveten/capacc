
compare_signal_strength <- function(mu, Q) {
  data.frame("type"            = c("cor", "iid"),
             "signal_strength" = c(sqrt(as.numeric(mu %*% Q %*% mu)),
                                   sqrt(as.numeric(mu %*% mu))))
}

#' @export
plot_signal_strength <- function(p, mu, rho, band = 2) {
  Q <- car_precision_mat(banded_neighbours(band, p), rho = rho)
  ss <- do.call("rbind", lapply(1:p, function(i) {
    mu_vec <- rep(0, p)
    mu_vec[1:i] <- mu
    cbind(compare_signal_strength(mu_vec, Q), data.frame("prop" = round(i/p, 4)))
  }))
  ggplot2::ggplot(data = ss, ggplot2::aes(prop, signal_strength, colour = type)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous("Signal strength", limits = c(0, max(ss$signal_strength))) +
    ggplot2::scale_x_continuous("Proportion") +
    ggplot2::ggtitle(paste0("2-banded, p = ", p, ", rho = ", rho, ", mu = ", mu))
}
