#' @param X A data matrix with p columns and n = (e - s) rows.
#' @param A A precision matrix.
#' @param l The maximum size of the subset.
optimise_subset_AR1 <- function(X, A, l = p) {
  p <- ncol(X)
  n <- nrow(X)

  assert_cov_mat(A)
  assert_equal_ncol(X, A)
  assert_integer_in_interval(l, c(2, p))
  # Restrict n to be greater than p?

  # Initialise B1 and B0.
  S <- get_S(X, A)
  B1 <- matrix(NA, nrow = p, ncol = p)
  B1[, 1] <- S$single
  B0 <- matrix(NA, nrow = p, ncol = p)
  B0[2:p, 1] <- cummax(B1[-p, 1])

  # Initialise J1 and J.
  J1 <- array(0, dim = c(p, p, p)) # One {0, 1}^p-vector at [d, k].
  J1[, , 1] <- diag(1, p)
  J <- array(0, dim = c(p, p, p)) # One {0, 1}^p-vector at [d, k].
  for (d in 1:p) {
    update_J_indicator <- which.max(c(B0[d, 1], B1[d, 1]))
    if (update_J_indicator == 1) J[, d, 1] <- J[, d - 1, 1]
    if (update_J_indicator == 2) J[d, d, 1] <- 1
  }

  # Update B1, B0, J1 and J for the remaining subset sizes k.
  for (k in 2:l) {
    for (d in k:p) {
      current_potential_savings <- c(B0[d - 1, k - 1] + S$single[d],
                                     B1[d - 1, k - 1] + S$interaction[d])
      B1[d, k] <- max(current_potential_savings, na.rm = TRUE)
      if(k < p && d >= k + 1)
        B0[d, k] <- max(B0[d - 1, k], B1[d - 1, k], na.rm = TRUE)

      update_J1_indicator <- which.max(current_potential_savings)
      if (update_J1_indicator == 1) {
        J1[, d, k] <- J[, d - 2, k - 1]
      } else if(update_J1_indicator == 2) {
        J1[, d, k] <- J1[, d - 1, k - 1]
      }
      J1[d, d, k] <- 1

      update_J_indicator <- which.max(c(B0[d, k], B1[d, k]))
      if (update_J_indicator == 1) {
        J[, d, k] <- J[, d - 1, k]
      } else if (update_J_indicator == 2) {
        J[, d, k] <- J1[, d, k]
      }
    }
  }
  B <- pmax(B0, B1, na.rm = TRUE)

  list('B0' = B0, 'B1' = B1, 'B' = B, 'J1' = J1, 'J' = J)
}

get_S <- function(X, A) {
  p <- ncol(X)
  n <- nrow(X)
  mean_x <- colMeans(X, na.rm = TRUE)
  S_single <- mean_x^2 * diag(A)
  S_interaction <- S_single[-1] + 2 * off_diag(1, A) * mean_x[-1] * mean_x[-p]
  list('single' = S_single, 'interaction' = c(0, S_interaction))
}

normalised_savings <- function(J, X, A) {
  mean_x <- colMeans(X, na.rm = TRUE)
  as.numeric(mean_x[J] %*% A[J, J] %*% mean_x[J])
}

penalised_savings <- function(X_subset, A) {
  n <- nrow(X_subset)
  p <- ncol(X_subset)
  beta <- adjusted_penalty(2, n, p)
  optimal_subset_res <- optimise_subset_AR1(X_subset, A)
  optimal_savings <- optimal_subset_res$B[p, ]
  savings_k <- n * optimal_savings  - cumsum(beta)
  optimal_k <- which.max(savings_k)
  optimal_J <- which(as.logical(optimal_subset_res$J[, p, optimal_k]))
  optimal_savings <- max(savings_k)
  return(list('k_max' = optimal_k, 'J_max' = optimal_J, 'S_max' = optimal_savings))
}
