compare_signal_strength <- function(vartheta, J, shape, Q) {
  p <- ncol(Q)
  theta <- generate_change(vartheta, J, shape, solve(Q))
  # theta <- c(theta, rep(0, ncol(Q) - length(theta)))
  data.table("cost"            = c("cor", "iid"),
             "signal_strength" = c(sqrt(as.numeric(theta %*% Q %*% theta)),
                                   sqrt(as.numeric(theta %*% theta))),
             "vartheta"        = vartheta,
             "k"               = length(J),
             "shape"           = shape)
}

expected_savings <- function(vartheta, J, shape, Q, n, seed = NA, n_sim = 100) {
  expected_saving <- function(Q) {
    theta <- generate_change(vartheta, J, shape, solve(Q), seed)
    p + n * as.numeric(theta %*% Q %*% theta)
  }

  p <- ncol(Q)
  if (shape %in% 4:10)
    expected_saving_cor <- mean(replicate(n_sim, expected_saving(Q)))
  else
    expected_saving_cor <- expected_saving(Q)
  data.table("cost"               = c("cor", "iid"),
             "expected_saving_H1" = c(expected_saving_cor,
                                      expected_saving(diag(p))),
             "expected_saving_H0" = rep(p, 2),
             "vartheta"           = vartheta,
             "k"                  = length(J),
             "shape"              = shape)
}

rename_shape <- function(dt) {
  dt <- dt[(cost == "iid" & shape == 0) | cost == "cor"]
  dt[cost == "iid", c("shape_str", "shape") := list("all", NA)]
  dt[shape == 0, shape_str := "equal"]
  dt[shape == 1, shape_str := "linear"]
  dt[shape == 2, shape_str := "quadratic"]
  dt[shape == 3, shape_str := "sqrt_inv"]
  dt[shape == 4, shape_str := "chisq_2"]
  dt[shape == 5, shape_str := "norm_iid"]
  dt[shape == 6, shape_str := "norm_cor"]
  dt[shape == 7, shape_str := "norm_cor07"]
  dt[shape == 8, shape_str := "norm_cor08"]
  dt[shape == 9, shape_str := "norm_cor09"]
  dt[shape == 10, shape_str := "norm_cor099"]
  dt$shape <- dt$shape_str
  dt_shape_str <- NULL
  dt
}

#' @export
plot_expected_saving <- function(p, vartheta = 1, precision_type = "banded",
                                 rho = 0.9, band = 2, change_type = "adjacent",
                                 seed = NA, n_sim = 100) {
  if (precision_type == "global_const")
    Q <- solve(constant_cor_mat(p, rho))
  else
    Q <- block_precision_mat(p                 = p,
                             m                 = p,
                             rho               = rho,
                             within_block_type = precision_type,
                             band              = band)
  n <- 1
  expect_dt <- do.call("rbind", lapply(0:10, function(shape) {
    do.call("rbind", lapply(1:p, function(k) {
      J <- get_affected_dims(change_type, k/p, p)
      expected_savings(vartheta, J, shape, Q, n, seed, n_sim)
    }))
  }))
  expect_dt <- rename_shape(expect_dt)
  expect_dt[, "expected_diff" := expected_saving_H1 - expected_saving_H0]

  ggplot2::ggplot(data = expect_dt, ggplot2::aes(k, expected_diff,
                                                 colour = interaction(cost, shape))) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous(expression(paste(E, "[", "S | ", H[1], "] - ", E, "[", "S | ", H[0], "]")),
                                limits = c(0, max(expect_dt$expected_diff))) +
    ggplot2::scale_x_continuous("|J|", breaks = round(seq.int(0, p, length.out = 11))) +
    ggplot2::ggtitle(paste0("2-banded, rho=", rho,  ", n=", n, ", p=", p,
                            ", vartheta=", vartheta, ", ", change_type, " changes")) +
    ggplot2::labs(colour = "cost.shape")
}
