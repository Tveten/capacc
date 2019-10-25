ar_cor_mat <- function(p, phi) {
  stopifnot(phi < 1 && phi > -1)

  Sigma <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      Sigma[i, j] <- phi^(abs(i - j))
    }
  }
  Sigma
}

constant_cor_mat <- function(p, rho) {
  stopifnot(rho > 0 && rho < 1)

  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  Sigma
}

ar_precision_mat <- function(p, phi) {
  Sigma_inv <- matrix(0, nrow = p, ncol = p)
  diag(Sigma_inv) <- 1 + phi^2
  diag(Sigma_inv)[c(1, p)] <- 1

  delta <- row(Sigma_inv) - col(Sigma_inv)
  Sigma_inv[off_diag_ind(c(-1, 1), p)] <- -phi

  Sigma_inv / (1 - phi^2)
}

car_precision_mat <- function(p, band = 2, rho = 0.5, sigma = 1) {
  W <- matrix(0, nrow = p, ncol = p)
  W[off_diag_ind(c(1:band, -(1:band)), p)] <- 1
  eigen_values_W <- eigen(W)$values
  if (rho <= 1 / eigen_values_W[p])
    stop(paste0('rho must be greater than 1 / smallest_eigen_value(W) = ', 1 / eigen_values_W[p]))
  1 / sigma^2 * (diag(as.vector(W %*% rep(1, p))) - rho * W)
}
