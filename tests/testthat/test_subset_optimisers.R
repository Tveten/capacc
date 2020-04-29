test_savings_MLE <- function(mu_est = mu_MLE(A)) {
  p <- 6
  A <- car_precision_mat(banded_neighbours(1, p))
  x <- simulate_cor(p = p, Sigma = solve(A))
  J <- 1:3
  unlist(lapply(lapply(1:6, function(i) 1:i), function(J) {
    as.numeric(savings(J, x, A, mu_est))
  }))
}

test_that('MLE savings always increases as more components are added in affected subset', {
  for (i in 1:100) expect_true(all(diff(test_savings_MLE()) >= 0))
})

test_that('OP C++ implementation returns equally as brute force in R', {
  data <- init_data(n = 50, p = 6, locations = 20, durations = 10, vartheta = 2)
  bands <- 1:4
  props <- c(0, 1/6, 3/6, 1)
  rhos <- c(-0.35, -0.2, 0.2, 0.5, 0.9)
  OP_method <- method_params("cor", precision_est_struct = NA, b = 0.9)
  BF_method <- method_params("cor_BF", precision_est_struct = NA, b = 0.9)
  BF_exact_method <- method_params("cor_exact", precision_est_struct = NA, b = 0.9)
  ineq_tol <- 1e-6
  set.seed(5)
  for (k in seq_along(bands)) {
    for (j in seq_along(props)) {
      for (i in seq_along(rhos)) {
        data$band <- bands[k]
        data$proportions <- props[j]
        data$rho <- rhos[i]
        data <- init_data_(data)
        seed <- sample(1:10^3, 1)
        OP_res <- simulate_mvcapa_known(data, OP_method, seed)
        BF_res <- simulate_mvcapa_known(data, BF_method, seed)
        BF_exact_res <- simulate_mvcapa_known(data, BF_exact_method, seed)
        expect_equal(OP_res$S_max, BF_res$S_max)
        expect_equal(OP_res$J_max, BF_res$J_max)
        expect_true(BF_exact_res$S_max >= OP_res$S_max - ineq_tol)
      }
    }
  }

  data$precision_type <- "random"
  data <- init_data_(data)
  for (j in seq_along(props)) {
    for (i in seq_along(rhos)) {
      for (seed in 1:10) {
        data$proportions <- props[j]
        data$rho <- rhos[i]
        set.seed(seed + 3)
        data <- init_data_(data)
        OP_res <- simulate_mvcapa_known(data, OP_method, seed)
        BF_res <- simulate_mvcapa_known(data, BF_method, seed)
        BF_exact_res <- simulate_mvcapa_known(data, BF_exact_method, seed)
        expect_equal(OP_res$S_max, BF_res$S_max)
        expect_equal(OP_res$J_max, BF_res$J_max)
        expect_true(BF_exact_res$S_max >= OP_res$S_max - ineq_tol)
      }
    }
  }
})
