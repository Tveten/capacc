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

compare_OP_BF <- function(n = 50, p = 5, r = 2, rho = 0.9,
                          prop = 0.1, cor_mat_type = 'random', mu = 1) {
  if (cor_mat_type == 'random')
    Q <- cuthill_mckee(car_precision_mat(random_neighbours(p), rho))
  else if (cor_mat_type == 'banded')
    Q <- car_precision_mat(banded_neighbours(r, p), rho)
  Q_sparse <- Matrix::Matrix(Q, sparse = TRUE)

  lower_nbs <- lower_nbs(Q_sparse)
  extended_nbs <- extended_lower_nbs(lower_nbs)
  x <- simulate_cor(n, p, mu, solve(Q), 1, n - 2, prop)
  n_full <- 1000
  sparse_penalty <- get_penalty('sparse', n_full, p)
  dense_penalty <- get_penalty('dense', n_full, p)
  OP_res <- optimise_mvnormal_saving(x, Q_sparse, lower_nbs, extended_nbs,
                                     dense_penalty$alpha, sparse_penalty$beta, sparse_penalty$alpha)
  BF_res <- optim_penalised_savings_BF(n_full, Q)(x)
  list('BF' = list('B_max' = BF_res$B_max, 'J_max' = BF_res$J_max),
       'OP' = list('B_max' = OP_res$B_max, 'J_max' = OP_res$J_max))
}

test_that('OP C++ implementation returns equally as brute force in R', {
  n <- 100
  p <- 10
  r <- 1:4
  prop <- c(0, 0.2, 0.5, 0.8, 1)
  rho <- c(-0.35, -0.2, 0.2, 0.5, 0.9)
  set.seed(5)
  for (k in seq_along(r)) {
    for (j in seq_along(prop)) {
      for (i in seq_along(rho)) {
        # print(c(k, j, i))
        res <- compare_OP_BF(n, p, r[k], rho[i], prop[j], 'banded')
        # print(res)
        expect_equal(res$OP$B_max, res$BF$B_max)
        expect_equal(res$OP$J_max, res$BF$J_max)
      }
    }
  }
  for (j in seq_along(prop)) {
    for (i in seq_along(rho)) {
      for (k in 1:100) {
        set.seed(k)
        res <- compare_OP_BF(n, p, rho = rho[i], prop = prop[j],
                             cor_mat_type = 'random')
        expect_equal(res$OP$B_max, res$BF$B_max)
        expect_equal(res$OP$J_max, res$BF$J_max)
      }
    }
  }
})
