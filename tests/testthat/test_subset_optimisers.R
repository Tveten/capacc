test_get_S <- function(n = 10, p = 4, phi = 0.5) {
  X <- matrix(2, nrow = n, ncol = p)
  A <- ar_precision_mat(p, phi)
  list(colMeans(X), A, get_S(X, A))
}

test_optimise_AR1 <- function(n = 10, p = 5, phi = 0.5, prop = 0.2, mu = c(0, 1),
                              return_all = FALSE) {
  X <- matrix(mu[1], nrow = n, ncol = round((1 - prop) * p))
  X <- cbind(X, matrix(mu[2], nrow = n, ncol = round(prop * p)))
  A <- ar_precision_mat(p, phi)
  optim_subset_res <- optimise_subset_AR1_SN(X, A)
  if (return_all) return(list('X' = X, 'A' = A, 'optim_subset_res' = optim_subset_res))
  else return(optim_subset_res)
}

expect_all_correct_subset_sizes <- function(J) {
  for (k in 1:dim(J)[3]) {
    for (d in k:dim(J)[2]) {
      expect_true(sum(J[, d, k]) == k)
    }
  }
}

test_that('Subsets are of correct size', {
  expect_all_correct_subset_sizes(test_optimise_AR1(10, 5, 0.5, 0.2, c(0, 1))$J)
  expect_all_correct_subset_sizes(test_optimise_AR1(10, 5, 0.5, 0.2, c(0.1, 1))$J)
  expect_all_correct_subset_sizes(test_optimise_AR1(100, 50, 0.9, 0.05, c(0.04, 1))$J)
})

expect_identical_vectors <- function(x, y) {
  expect_true(all(x == y))
}

test_that('Optimal subsets are reasonable', {
  n <- 20
  p <- 10
  J1 <- test_optimise_AR1(n, p, phi = 0.2)$J
  expect_identical_vectors(which(as.logical(J1[, p, 2])), c(9, 10))
  expect_identical_vectors(which(as.logical(J1[, p, 3])), c(1, 9, 10))

  J2 <- test_optimise_AR1(n, p, phi = 0.99)$J
  expect_identical_vectors(which(as.logical(J2[, p, 2])), c(1, 9))
  expect_identical_vectors(which(as.logical(J2[, p, 3])), c(1, 2, 9))

  J3 <- test_optimise_AR1(n, p, phi = 0.4, mu = c(0.1, 1))$J
  expect_identical_vectors(which(as.logical(J3[, p, 2])), c(9, 10))
  expect_identical_vectors(which(as.logical(J3[, p, 3])), c(2, 9, 10))
  expect_identical_vectors(which(as.logical(J3[, p, 4])), c(2, 4, 9, 10))

  J4 <- test_optimise_AR1(n, p, phi = 0.4, prop = 0.8, mu = c(1, 0.1))$J
  expect_identical_vectors(which(as.logical(J4[, p, 2])), c(1, 2))
  expect_identical_vectors(which(as.logical(J4[, p, 3])), c(1, 2, 4))
  expect_identical_vectors(which(as.logical(J4[, p, 4])), c(1, 2, 4, 6))

  J5 <- test_optimise_AR1(n, p, phi = 0.6, mu = c(0.1, 1))$J
  expect_identical_vectors(which(as.logical(J5[, p, 2])), c(2, 9))
  expect_identical_vectors(which(as.logical(J5[, p, 3])), c(2, 4, 9))
  expect_identical_vectors(which(as.logical(J5[, p, 4])), c(2, 4, 6, 9))

  J6 <- test_optimise_AR1(n, p, phi = 0.6, prop = 0.8, mu = c(1, 0.1))$J
  expect_identical_vectors(which(as.logical(J6[, p, 2])), c(2, 4))
  expect_identical_vectors(which(as.logical(J6[, p, 3])), c(2, 4, 6))
  expect_identical_vectors(which(as.logical(J6[, p, 4])), c(2, 4, 6, 8))
})

test_that('B matrix is of correct format', {
  p <- 5
  B <- test_optimise_AR1(p = p)$B
  expect_true(all(is.na(off_diag(-(1:4), B))))
  expect_true(all(!is.na(off_diag(0:4, B))))
})

test_that('J is argmax corresponding to max in B', {
  n <- 10
  p <- 5
  res <- test_optimise_AR1(n = n, p = p, mu = c(0.1, 1), return_all = TRUE)
  B <- res$optim_subset_res$B
  J <- res$optim_subset_res$J
  for (k in 1:p) {
    for (d in k:p) {
      compare_savings <- c(normalised_savings(res$X, res$A, which(as.logical(J[, d, k]))), B[d, k])
      expect_equal(compare_savings[1], compare_savings[2])
    }
  }
})

compare_OP_BF <- function(method = 'ar', n = 100, p = 10, r = 2, rho = 0.9,
                          prop = 0.8, regime = 'sparse', mu = 1, set_zero = 0,
                          return_all = FALSE) {
  A <- car_precision_mat(p, r, rho)
  x <- simulate_cor(n, p, mu, solve(A), 1, n - 2, prop)
  if (set_zero != 0) {
    A[set_zero, set_zero - 1] <- 0
    A[set_zero - 1, set_zero] <- 0
  }
  optim_BF <- penalised_savings_BF(x, n, A, regime)
  if (method == 'ar') optim_CAR_OP <- penalised_savings_ar(x, n, A, regime)
  if (method == 'car') optim_CAR_OP <- penalised_savings_car(x, n, A, regime)
  if (return_all) return(list('x' = x, 'A' = A, 'BF' = optim_BF, 'OP' = optim_CAR_OP))
  else return(list('BF' = optim_BF, 'OP' = optim_CAR_OP))
}

test_that('AR and CAR OP returns equally as brute force', {
  n <- 100
  p <- 10
  r <- 1:4
  set.seed(5)
  rho <- c(-0.4, -0.2, 0.2, 0.5, 0.9)
  prop <- c(0.2, 0.5, 0.8, 1)
  regime <- 'sparse'
  for (k in seq_along(r)) {
    for (j in seq_along(prop)) {
      for (i in seq_along(rho)) {
        res_ar <- compare_OP_BF('ar', n, p, r[k], rho[i], prop[j], regime)
        expect_equal(res_ar$OP$B_max, res_ar$BF$B_max)
        expect_equal(res_ar$OP$J_max, res_ar$BF$J_max)

        res_car <- compare_OP_BF('car', n, p, r[k], rho[i], prop[j], regime)
        expect_equal(res_car$OP$B_max, res_car$BF$B_max)
        expect_equal(res_car$OP$J_max, res_car$BF$J_max)

        # res2 <- compare_penalised_savings_AR1(n, p, rho[i], prop[j], set_zero = 4)
        # expect_equal(res2$OP$B_max, res2$SN$B_max)
        # expect_equal(res2$OP$J_max, res2$SN$J_max)
        #
        # res3 <- compare_penalised_savings_AR1(n, p, rho[i], prop[j], set_zero = 2)
        # expect_equal(res3$OP$B_max, res3$SN$B_max)
        # expect_equal(res3$OP$J_max, res3$SN$J_max)
      }
    }
  }
})
