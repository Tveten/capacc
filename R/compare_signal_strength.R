compare_signal_strength <- function(vartheta, k, shape, Q) {
  p <- ncol(Q)
  theta <- generate_change(vartheta, k, shape)
  theta <- c(theta, rep(0, ncol(Q) - length(theta)))
  data.table("cost"            = c("cor", "iid"),
             "signal_strength" = c(sqrt(as.numeric(theta %*% Q %*% theta)),
                                   sqrt(as.numeric(theta %*% theta))),
             "vartheta"        = vartheta,
             "k"               = k,
             "shape"           = shape)
}
#' @export
plot_signal_strength <- function(p, vartheta = 1, shape = 0, rho = 0.9, band = 2) {
  Q <- car_precision_mat(banded_neighbours(band, p), rho = rho)
  ss <- do.call("rbind", lapply(0:5, function(shape) {
    do.call("rbind", lapply(1:p, function(k) {
      compare_signal_strength(vartheta, k, shape, Q)
    }))
  }))
  ss <- ss[(cost == "iid" & shape == 0) | cost == "cor"]
  ggplot2::ggplot(data = ss, ggplot2::aes(k, signal_strength, colour = interaction(cost, shape))) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous("Signal strength", limits = c(0, max(ss$signal_strength))) +
    ggplot2::scale_x_continuous("Proportion") +
    ggplot2::ggtitle(paste0("2-banded, p=", p, ", rho=", rho,
                            ", vartheta=", vartheta, ", shape=", shape))
}
