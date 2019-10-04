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
