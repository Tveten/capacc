## usethis namespace: start
#' @useDynLib mvcapaCor, .registration = TRUE
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

####
#### Helper function. ----------------------------------------------------------
####
get_Q <- function(x, A, sparse = FALSE) {
  mean_x <- matrix(colMeans(x, na.rm = TRUE), ncol = ncol(x), nrow = 1)
  Q <- t(mean_x) %*% mean_x * A
  if (sparse) return(Matrix::Matrix(Q, sparse = TRUE))
  else return(Q)
}

get_w <- function(x, A) {
  mean_x <- colMeans(x, na.rm = TRUE)  # Row vector
  as.numeric(mean_x * (A %*% mean_x))
}

get_neighbours_less_than <- function(i, sparse_mat) {
  nbs <- which(!is_zero(sparse_mat[i, ]))
  rev(nbs[nbs < i])
}

remaining_neighbours_below <- function(d, nb_list, p) {
  if (d <= 1) integer(0)
  else if (d > p) stop('d must be between (or equal to) 1 and p.')
  else return(rev(intersect(1:(d - 1), unlist(nb_list[d:p]))))
}

band <- function(x) {
  distance_from_diagonal <- rep(0, nrow(x) - 1)
  for (i in 2:nrow(x)) {
    lower_nonzeros <- which(!is_zero(x[i, 1:(i - 1)]))
    if (length(lower_nonzeros) > 0)
      distance_from_diagonal[i] <- i - min(lower_nonzeros)
    else
      distance_from_diagonal[i] <- 0
  }
  max(distance_from_diagonal)
}

indicator_power_set <- function(size, include_empty_set = TRUE, as_list = FALSE) {
  if (size >= 1) {
    set <- 1:size
    ind_power_set <- vapply(2^(1:size - 1), function(k) {
      rep(c(0, 1), each = k, times = 2^size / (2 * k))
    }, numeric(2^size))
    if (!include_empty_set) ind_power_set <- ind_power_set[-1, ]
    if (as_list) ind_power_set <- split(ind_power_set, seq(nrow(ind_power_set)))
    return(ind_power_set)
  } else return(integer(0))
}

power_set <- function(size, include_empty_set = FALSE) {
  set <- 1:size
  p_set <- lapply(indicator_power_set(size, as_list = TRUE), function(u) {
    set[as.logical(u)]
  })
  names(p_set) <- NULL
  if (include_empty_set) return(p_set)
  else return(p_set[-1])
}

####
#### Penalized savings optimisers. ---------------------------------------------
####

mu_MLE <- function(Q) {
  Q <- as.matrix(Q)
  function(mean_x, J) {
    W <- solve(Q[J, J, drop = FALSE]) %*% Q[J, -J, drop = FALSE]
    mean_x[J] + W %*% mean_x[-J]
  }
}

mu_aMLE <- function() {
  function(mean_x, J) {
    mean_x[J]
  }
}

savings <- function(J, mean_x, n, Q, mu_est = mu_MLE(Q)) {
  if (length(J) == 0) return(0)
  else {
    mu_hat <- rep(0, length(mean_x))
    mu_hat[J] <- mu_est(mean_x, J)
    return(as.numeric(n * t(2 * mean_x - mu_hat) %*% Q %*% mu_hat))
  }
}

optimal_mvnormal_dense_saving <- function(x, Q, mu, alpha) {
  n <- nrow(x)
  x_bar <- colMeans(x, na.rm = TRUE)
  saving <- as.numeric(n * t(2 * x_bar - mu) %*% Q %*% mu - alpha)
  expected_saving <- - as.numeric(n * t(mu) %*% Q %*% mu)
  saving - expected_saving
}

savings_difference <- function(J, x, Q) {
  if (length(J) > 0) return(savings(J, x, Q) - savings(J, x, Q, mu_aMLE()))
  else return(0)
}

optimise_mvnormal_iid_saving <- function(x, penalty) {
  mean_x2 <- colMeans(x, na.rm = TRUE)^2
  order_mean_x2 <- order(mean_x2, decreasing = TRUE)
  sorted_mean_x2 <- mean_x2[order_mean_x2]
  savings_k <- cumsum(nrow(x) * sorted_mean_x2) - penalty
  optimal_savings <- max(savings_k)
  optimal_k <- which.max(savings_k)
  optimal_J <- order_mean_x2[1:optimal_k]
  optimal_u <- rep(0, ncol(x))
  optimal_u[optimal_J] <- 1
  list("S_max" = optimal_savings, 'J_max' = optimal_J, 'u_max' = optimal_u)
}

optimise_mvnormal_saving_BF <- function(x, Q, penalty, mu_est = mu_aMLE()) {
  p <- ncol(x)
  P_J <- power_set(p, include_empty_set = TRUE)
  lengths <- vapply(P_J, length, numeric(1))
  P_J <- P_J[lengths <= penalty$k_star | lengths == p]
  mean_x <- colMeans(x, na.rm = TRUE)
  S <- vapply(P_J, function(J) {
    savings(J, mean_x, nrow(x), Q, mu_est) - penalty$vec[length(J) + 1]
  }, numeric(1))
  S_max <- max(S)
  S_which_max <- which.max(S)
  J_max <- P_J[[S_which_max]]
  u_max <- rep(0, p)
  u_max[J_max] <- 1
  list('S_max' = S_max, 'J_max' = J_max, 'u_max' = u_max, 'S' = S)
}

optim_penalised_savings_c <- function(x, n_full, precision_mat_obj, penalty_regime = 'combined',
                                      b = 1, sparse = FALSE, adjusted = FALSE) {
  optim_func <- function(penalty, sparse) {
    if (sparse) return(BQP_optim_sparse(Q, w, penalty, precision_mat_obj$extended_nbs))
    else return(BQP_optim(Q, w, penalty, precision_mat_obj$extended_nbs))
  }

  n <- nrow(x)
  p <- ncol(x)
  Q <- n * get_Q(x, precision_mat_obj$A, sparse = sparse)
  w <- n * get_w(x, precision_mat_obj$A)

  if (penalty_regime == 'combined') {
    sparse_penalty <- get_penalty('sparse', n_full, p, b)
    sparse_res <- optim_func(sparse_penalty, sparse)
    dense_penalty <- get_penalty('dense', n_full, p, b)
    dense_B_max <- savings(1:p, x, precision_mat_obj$A) - dense_penalty$alpha
    dense_res <- list('B_max' = dense_B_max, 'J_max' = p:1)
    if (adjusted) {
      if (sparse_res$k_min > sparse_penalty$k_star) return(dense_res)
      else {
        savings_diff <- savings_difference(sparse_res$J_max, x, precision_mat_obj$A)
        sparse_res$B_max <- sparse_res$B_max + savings_diff
        if (sparse_res$B_max >= dense_B_max) return(sparse_res)
        else return(dense_res)
      }
    } else {
      maximiser <- which.max(c(sparse_res$B_max, dense_res$B_max))
      return(list(sparse_res, dense_res)[[maximiser]])
    }
  } else if (penalty_regime == 'sparse_point') {
    point_penalty <- get_penalty(penalty_regime, n_full, p, b)
    return(optim_func(point_penalty, sparse))
  }
}

optim_mvnormal_savings_test <- function(n = 50, p = 10, r = 2, rho = 0.9,
                                        prop = 1/p, mu = 1, compare = FALSE,
                                        benchmark = TRUE) {
  n_full <- 1000
  Q <- car_precision_mat(banded_neighbours(r, p), rho)
  Q_sparse <- car_precision_mat(banded_neighbours(r, p), rho, sparse = TRUE)
  Q_obj <- init_precision_mat(Q_sparse)
  sparse_penalty <- get_penalty('sparse', n_full, p)
  dense_penalty <- get_penalty('dense', n_full, p)
  x <- simulate_cor(n, p, mu, solve(Q), 1, n - 2, prop)
  if (compare) {
    C_res <- optimise_mvnormal_saving(x, Q_sparse, Q_obj$nbs,
                                      Q_obj$extended_nbs,
                                      dense_penalty$alpha,
                                      sparse_penalty$beta,
                                      sparse_penalty$alpha)
    print('C++')
    print(C_res)
    R_res <- optim_penalised_savings_BF(n_full, Q)(x)
    print('R')
    print(list(R_res$S_max, R_res$J_max))
  }

  if (benchmark) {
    microbenchmark::microbenchmark(
      optimise_mvnormal_saving(x, Q_sparse, Q_obj$nbs, Q_obj$extended_nbs, dense_penalty$alpha, sparse_penalty$beta, sparse_penalty$alpha),
    # penalised_savings_car_aMLE(x, n_full, precision_mat_obj),
      optim_penalised_savings_iid(n_full)(x),
    times = 1000)
  }
}

optim_penalised_savings_c_test <- function(n = 50, p = 5, r = 2, rho = 0.9,
                                           prop = 1/p, mu = 1) {
  n_full <- 1000
  A <- car_precision_mat(banded_neighbours(r, p), rho)
  A_sparse <- car_precision_mat(banded_neighbours(r, p), rho, sparse = TRUE)
  # A <- car_precision_mat(lattice_neighbours(p), rho)
  precision_mat_obj <- init_precision_mat(A)
  x <- simulate_cor(n, p, mu, solve(A), 1, n - 2, prop)
  C_res <- optim_penalised_savings_c(x, n_full, precision_mat_obj)

  print('C++')
  print(C_res)
  R_res <- penalised_savings_car_aMLE(x, n_full, precision_mat_obj)
  R_res <- optim_penalised_savings_BF(n_full, precision_mat_obj$A)(x)
  print('R')
  print(list(R_res$S_max, R_res$J_max))
  Q <- get_Q(x, precision_mat_obj$A)
  Q_sparse <- get_Q(x, precision_mat_obj$A, sparse = TRUE)
  w <- get_w(x, precision_mat_obj$A)
  sparse_penalty <- get_penalty('sparse', n_full, p)
  dense_penalty <- get_penalty('dense', n_full, p)
  BQP_optim(Q, w, sparse_penalty, precision_mat_obj$extended_nbs)
  print(list("dense" = dense_penalty))
  # optimise_mvnormal_savings(x, A_sparse, precision_mat_obj$extended_nbs, list("dense" = dense_penalty))

  # microbenchmark::microbenchmark(
    # BQP_optim(Q, w, penalty, precision_mat_obj$extended_nbs),
    # optim_penalised_savings_c(x, n_full, precision_mat_obj),
    # {get_Q(x, precision_mat_obj$A); get_w(x, precision_mat_obj$A)},
    # savings(1:p, x, A),
    # BQP_optim_sparse(Q_sparse, w, sparse_penalty, precision_mat_obj$extended_nbs),
    # optim_penalised_savings_c(x, n_full, precision_mat_obj, sparse = TRUE),
    # {get_Q(x, precision_mat_obj$A, sparse = TRUE); get_w(x, precision_mat_obj$A)},
    # # penalised_savings_car_aMLE(x, n_full, precision_mat_obj),
    # optim_penalised_savings_iid(n_full)(x)
  # )
}

optim_test <- function(n = 50, p = 5, r = 2, rho = 0.9, prop = 0.1, mu = 1) {
  n_full <- 1000
  A <- car_precision_mat(banded_neighbours(r, p), rho)
  # A <- car_precision_mat(lattice_neighbours(p), rho)
  precision_mat_obj <- init_precision_mat(A)
  x <- simulate_cor(n, p, mu, solve(A), 1, n - 2, prop)
  n <- nrow(x)
  p <- ncol(x)
  Q <- n * get_Q(x, precision_mat_obj$A)
  w <- n * get_w(x, precision_mat_obj$A)
  print(Q)
  print(test_subsetting(Q, 1, 2))
}

optim_savings_memory_test <- function(n = 500, p = 100, r = 2, rho = 0.9,
                                      prop = 0.1, mu = 1, n_sim = 10e6) {
  n_full <- 1000
  A <- car_precision_mat(banded_neighbours(r, p), rho)
  lower_nbs <- lapply(1:p, get_neighbours_less_than, sparse_mat = A)
  extended_nbs <- lapply(1:p, function(d) remaining_neighbours_below(d, lower_nbs, p))
  x <- simulate_cor(n, p, mu, solve(A), 1, n - 2, prop)
  penalty <- get_penalty('sparse', n_full, p)
  print(A)
  Q <- n * get_Q(x, A)
  w <- n * get_w(x, A)
  # for (i in 1:n_sim) {
  #   print(i)
  #   optimise_savings(Q, w, penalty, extended_nbs)
  # }
}

####
#### Testing functions ------------------------------------------
####
init_setup_MLE <- function(proportions = 0.3, mu = 0.7, n = 100, p = 8,
                           locations = n - durations - 1, durations = 10,
                           change_type = 'adjacent') {
  return(list('n'           = n,
              'p'           = p,
              'proportions' = proportions,
              'mu'          = mu,
              'locations'   = locations,
              'durations'   = durations,
              'change_type' = change_type))
}

MLE_compare_dt <- function(n = 12, p = 10, r = 2, rho = 0.9, prop = 0.3, mu = 1) {
  n_full <- 1000
  A <- car_precision_mat(banded_neighbours(r, p), rho)
  x <- simulate_cor(n, p, mu, solve(A), 1, n - 2, prop)
  print(round(A, 2))
  print(round(solve(A), 2))
  print(eigen(A, symmetric = TRUE)$values)
  mean_x <- colMeans(x)
  methods <- c('MLE', 'aMLE')
  res_list <- lapply(methods, function(method) {
    if (method == 'MLE') mu_func <- mu_MLE(A)
    if (method == 'aMLE') mu_func <- mu_aMLE()
    res <- optim_penalised_savings_BF(n_full, A, mu_est = mu_func)(x)
    res_mat <- cbind(res$B, 1:length(res$B),
                     indicator_power_set(p, include_empty_set = FALSE))
    res_dt <- as.data.table(res_mat)
    colnames(res_dt) <- c('savings', 'set_ind', paste0('u_', 1:p))
    res_dt$method <- method
    res_dt
  })
  res_dt <- do.call('rbind', res_list)
  u_max <- unlist(res_dt[which.max(savings)][, 3:(3 + p - 1)])
  J_max <- (1:p)[as.logical(u_max)]
  print(round(mean_x, 4))
  print(as.vector(mu_MLE(A)(mean_x, J_max)))
  print(solve(A[J_max, J_max]))
  u_max_aMLE <- unlist(res_dt[method == 'aMLE'][which.max(savings)][, 3:(3 + p - 1)])
  J_max_aMLE <- (1:p)[as.logical(u_max_aMLE)]
  W <- solve(A[J_max, J_max]) %*% A[J_max, -J_max]
  aMLE_bound(mean_x, J_max, A, W, n)
  # print(rowSums(W))
  res_dt
}

aMLE_bound <- function(mean_x, Q,  J_hat, n) {
  p <- length(mean_x)
  W <- solve(Q[J_hat, J_hat]) %*% Q[J_hat, -J_hat, drop = FALSE]
  W_extended <- matrix(0, nrow = p, ncol = p)
  W_extended[J_hat, -J_hat] <- W
  # sigma_max <- svd(A %*% W_extended)$d[1]
  # lambda_max <- eigen(Q %*% W_extended + t(Q %*% W_extended))$values[1]
  lambda_max <- eigen(Q %*% W_extended)$values[1]
  squared_norm_nonchange_mean_x <- sum(mean_x[-J_hat]^2)
  # print(paste0('sigma_max = ', sigma_max))
  # print(paste0('lambda_max = ', lambda_max))
  # print(paste0('Squared euclidean length of x(-J_hat) = ', squared_norm_nonchange_mean_x))
  savings_bound <- 2 * n * lambda_max * squared_norm_nonchange_mean_x
  # print(paste0('Savings bound = ', savings_bound))
  list("lambda_max" = lambda_max, "savings_bound" = savings_bound)
}

sim_aMLE_bound <- function(n, p, rho, size_J_hat, band = 2) {
  mean_x <- rnorm(p) / sqrt(n)
  Q <- car_precision_mat(banded_neighbours(band, p), rho = rho, sparse = FALSE)
  J_hat <- 1:size_J_hat
  aMLE_bound(mean_x, Q, J_hat, n)
}

MLE_compare <- function(n = 12, p = 10, rho = 0.9, prop = 0.3, mu = 1, r = 2) {
  compare_dt <- MLE_compare_dt(n = n, p = p, r = r, rho = rho, prop = prop, mu = mu)
  compare_dt[, 'diff_MLE_aMLE' := .SD[method == 'MLE']$savings - .SD[method == 'aMLE']$savings, set_ind]
  compare_dt <- compare_dt[order(set_ind), , ]
  show(ggplot2::qplot(diff_MLE_aMLE, data = compare_dt, geom = 'histogram'))
  # print(compare_dt)
  max_dt <- compare_dt[, .SD[which.max(savings)], method]
  max_diff <- max_dt[method == 'MLE']$savings - max_dt[method == 'aMLE']$savings
  print(paste0('MLE_max - aMLE_max: ', max_diff))
  # sum_aMLE_greater <- sum(compare_dt$diff_MLE_aMLE < 0)
  # print(paste0('Times aMLE >= MLE: ', sum_aMLE_greater))
  return(max_dt)
}

test_MLE_greater <- function(n_sim = 50) {
  p <- c(4, 7, 10)
  rho <- c(-0.3, 0.99)
  prop <- c(0.1, 0.3, 0.5, 0.8)
  mu <- c(0.2, 0.5, 1)
  r <- c(1:5)
  param_combs <- expand.grid(p, rho, prop, mu, r)
  colnames(param_combs) <- c('p', 'rho', 'prop', 'mu', 'r')
  res <- rep(0, nrow(param_combs))
  for (i in 1:nrow(param_combs)) {
    print(c(i, nrow(param_combs)))
    curr_params <- param_combs[i, ]
    print(curr_params)
    sub_res <- rep(0, n_sim)
    for (j in 1:n_sim) {
      sub_res[j] <- MLE_compare(4, curr_params$p, curr_params$rho, curr_params$prop, curr_params$mu, curr_params$r)
    }
    res[i] <- sum(sub_res)
  }
  # Result: 0.
  sum(res)
}

compare_BF_MLE <- function(n = 100, p = 6, r = 2, rho = 0.9, prop = 0.5, mu = 1) {
  A <- car_precision_mat(banded_neighbours(r, p), rho)
  x <- simulate_cor(n, p, mu, solve(A), 1, n - 2, prop)
  print(colMeans(x))
  res <- lapply(c(mu_MLE(A), mu_aMLE()), function(mu_func) {
    res <- optim_penalised_savings_BF(n, A, mu_est = mu_func)(x)
    c(res$B_max, res$u_max)
  })
  res <- as.data.frame(do.call('rbind', res))
  colnames(res) <- c('B_max', paste0('u_', 1:p))
  res <- cbind(data.frame('method' = c('MLE', 'aMLE')), res)
  res
}

compare_OP_BF2 <- function(method = 'ar', n = 100, p = 10, r = 2, rho = 0.9,
                          prop = 0.8, regime = 'sparse', mu = 1, set_zero = 0,
                          return_all = FALSE) {
  A <- car_precision_mat(banded_neighbours(p, r), rho)
  x <- simulate_cor(n, p, mu, solve(A), 1, n - 2, prop)
  if (set_zero != 0) {
    A[set_zero, set_zero - 1] <- 0
    A[set_zero - 1, set_zero] <- 0
  }
  optim_BF <- penalised_savings_BF(x, n, A, regime)
  if (method == 'ar') optim_CAR_OP <- penalised_savings_ar(x, n, A, regime)
  if (method == 'car') optim_CAR_OP <- penalised_savings_car(x, n, A, regime)
  if (return_all) return(list('x' = x, 'A' = A, 'BF' = optim_BF, 'OP' = optim_CAR_OP))
  else return(list('BF' = optim_BF$u_max, 'OP' = optim_CAR_OP$u_max))
}

compare_many_times <- function(method = 'ar', p = 10, r = 1, rho = 0.9, prop = 0.8) {
  compare_res <- lapply(1:100, function(i) {
    res <- compare_OP_BF2(method, p = p, r = r, rho = rho, prop = prop)
    c(res$BF, res$OP)
  })
  compare_res <- do.call('rbind', compare_res)
  colnames(compare_res) <- c('BF', 'OP')
  diff <- apply(compare_res, 1, function(res) res[2] - res[1])
  n_equal <- sum(abs(diff) <= sqrt(.Machine$double.eps))
  list('method' = method, 'max_diff' = max(abs(diff)), 'n_equal' = n_equal)
}

B_compare <- function() {
  n <- 100
  p <- 10
  r <- 1:4
  set.seed(5)
  rho <- c(-0.35, -0.2, 0.2, 0.5, 0.9)
  prop <- c(0.2, 0.5, 0.8, 1)
  regime <- 'dense'
  j <- 4
  i <- 5
  k <- 26
  set.seed(k)
  A <- cuthill_mckee(car_precision_mat(list('random', p), rho[i]))
  x <- simulate_cor(n, p, 1, solve(A), 1, n - 2, prop[j])
  res_bf <- penalised_savings_BF(x, n, A, regime)
  res_car <- penalised_savings_car(x, n, A, regime)
  print(res_bf$B_max - res_car$B_max)
}

time_test <- function(n = 50, p = 7, r = 2, rho = 0.9,
                      prop = 0.1, mu = 1) {
  n_full <- 1000
  set.seed(306)
  A <- car_precision_mat(lattice_neighbours(p), rho)
  x <- simulate_cor(n, p, mu, solve(A), 1, n - 2, prop)
  Q <- n * get_Q(x, A)
  w <- n * get_w(x, A)
  lower_nbs <- lapply(1:p, get_neighbours_less_than, sparse_mat = A)
  all_next_nbs <- lapply(1:p, function(d) remaining_neighbours_below(d, lower_nbs, p))
  penalty <- get_penalty('sparse', n_full, p)
  microbenchmark::microbenchmark(test_grow_tree(Q, w, penalty, all_next_nbs),
                                 times = 1000)
}

approx_savings_prunable <- function(x, A) {
    p <- ncol(x)
    n <- nrow(x)
    P_J <- power_set(p)
    B <- vapply(P_J, function(J) {
      savings(J, x, A, mu_aMLE())
    }, numeric(1))
    max_B <- max(B)
    B_split_max <- matrix(0, ncol = n - 1, nrow = 2)
    for (t in 1:(n - 1)) {
      B_split <- vapply(P_J, function(J) {
        c(savings(J, x[1:t, , drop = FALSE], A, mu_aMLE()),
          savings(J, x[(t + 1):n, , drop = FALSE], A, mu_aMLE()))
      }, numeric(2))
      B_split_max[, t] <- apply(B_split, 1, max)
    }
    sum_B <- colSums(B_split_max)
    # print(max_B)
    # print(sum_B)
    all(sum_B >= max_B)
}

run_prunable <- function(n = 100, p = 5, rho = 0.8, mu = 1, prop = 0) {
  A <- car_precision_mat(list('random', p), rho)
  x <- simulate_cor(n, p, mu, solve(A), 1, n - 2, prop)
  approx_savings_prunable(x, A)
}

####
#### Old functions -------------------------------------------------------------
####

normalised_savings <- function(x, A, J = 1:ncol(x)) {
  mean_x <- colMeans(x, na.rm = TRUE)
  as.numeric(mean_x[J] %*% A[J, J] %*% mean_x[J])
}

dense_savings2 <- function(x, n_full, A, b = 1) {
  n <- nrow(x)
  p <- ncol(x)
  beta <- const_penalty(b, 2, n_full, p)
  n * normalised_savings(x, A) - beta
}

penalised_savings_ar1 <- function(x_subset, n, A, b = 1) {
  n_subset <- nrow(x_subset)
  p <- ncol(x_subset)
  penalty <- penalty_func(b, 2, n, p)
  sum_beta <- penalty$sum_beta
  k_star <- penalty$k_star

  optimal_subset_res <- optimise_subset_AR1(x_subset, A, l = k_star)
  optimal_savings <- optimal_subset_res$B[p, ]

  savings_k <- rep(0, p)
  savings_k[1:k_star] <- n_subset * optimal_savings[1:k_star] - sum_beta[1:k_star]
  # savings_k[p] <- n_subset * normalised_savings(x_subset, A) - sum_beta[length(beta)]
  savings_k[p] <- dense_savings(x_subset, n, A, b)
  optimal_savings <- max(savings_k)
  optimal_k <- which.max(savings_k)

  if (optimal_k <= k_star) optimal_J <- optimal_subset_res$J[, p, optimal_k]
  else optimal_J <- rep(1, p)

  list('k_max' = optimal_k, 'u_max' = optimal_J, 'B_max' = optimal_savings,
       'S_obj' = optimal_subset_res$S, 'B0' = optimal_subset_res$B0,
       'B1' = optimal_subset_res$B1)
}

optimise_subset_AR1_SN <- function(x, A, l = p) {
  p <- ncol(x)
  n <- nrow(x)

  assert_cov_mat(A)
  assert_equal_ncol(x, A)
  assert_integer_in_interval(l, c(2, p))
  # Restrict n to be greater than p?

  # Initialise B1 and B0.
  S <- get_S(x, A)
  B1 <- matrix(NA, nrow = p, ncol = l)
  B1[, 1] <- S$single
  B0 <- matrix(NA, nrow = p, ncol = l)
  B0[2:p, 1] <- cummax(B1[-p, 1])

  # Initialise J1 and J.
  J1 <- array(0, dim = c(p, p, l)) # One {0, 1}^p-vector at [d, k].
  J1[, , 1] <- diag(1, p)
  J <- array(0, dim = c(p, p, l)) # One {0, 1}^p-vector at [d, k].
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

  list('B0' = B0, 'B1' = B1, 'B' = B, 'J1' = J1, 'J' = J, 'S' = S)
}

get_S <- function(x, A) {
  p <- ncol(x)
  n <- nrow(x)
  mean_x <- colMeans(x, na.rm = TRUE)
  S_single <- mean_x^2 * diag(A)
  S_interaction <- S_single[-1] + 2 * off_diag(1, A) * mean_x[-1] * mean_x[-p]
  list('single' = S_single, 'interaction' = c(0, S_interaction))
}


penalised_savings_ar1_OP <- function(x, n_full, A, b = 1) {
  p <- ncol(x)
  n <- nrow(x)

  assert_cov_mat(A)
  assert_equal_ncol(x, A)
  # Restrict n to be greater than p?

  # Initialise penalty.
  penalty <- penalty_per_comp(b, 2, n_full, p)
  beta <- penalty$beta
  k_star <- penalty$k_star

  # Initialise B1 and B/B0.
  S <- get_S(x, A)
  B1 <- rep(NA, p + 1)
  B1[2] <- n * S$single[1] - beta
  B0 <- rep(NA, p + 1)
  B0[2] <- 0

  # Initialise J1 and J/J0.
  J1 <- rbind(rep(NA, p), diag(1, p)) # (d, d) is always 1 by def of J1
  J <- matrix(0, nrow = p + 1, ncol = p) # One {0, 1}^p-vector at [d, k].
  if (which.max(c(B0[2], B1[2])) == 2) J[2, ] <- J1[2, ]

  # Update B1, B0, J1 and J for the remaining subset sizes k.
  d <- 3  # d is d - 1 along rows, but d along columns.
  while (d <= p + 1 && sum(J[d - 1, ]) <= k_star) {
    B0[d] <- max(B0[d - 1], B1[d - 1], na.rm = TRUE)
    alternative_savings_1 <- c(B0[d - 1] + n * S$single[d - 1],
                               B1[d - 1] + n * S$interaction[d - 1]) - beta
    B1[d] <- max(alternative_savings_1, na.rm = TRUE)

    update_J1_indicator <- which.max(alternative_savings_1)
    if (update_J1_indicator == 1) J1[d, ] <- J[d - 2, ]
    else if(update_J1_indicator == 2) J1[d, ] <- J1[d - 1, ]
    J1[d, d - 1] <- 1  # Rows start at 0, columns at 1.

    update_J_indicator <- which.max(c(B0[d], B1[d]))
    if (update_J_indicator == 1) J[d, ] <- J[d - 1, ]
    else if (update_J_indicator == 2) J[d, ] <- J1[d, ]

    d <- d + 1
  }

  if (d < p + 2 || sum(J[p + 1, ]) > k_star) {
    B_max <- dense_savings(x, n_full, A, b)
    u_max <- rep(1, p)
  } else {
    B_max <- max(B0[p + 1], B1[p + 1], na.rm = TRUE)
    u_max <- J[p + 1, ]
  }

  list('B_max' = B_max, 'u_max' = u_max,
       'B0' = B0, 'B1' = B1, 'J1' = J1, 'J' = J, 'S' = S)
}

penalised_savings_ar1_NB1 <- function(x, n_full, A, b = 1) {
  p <- ncol(x)
  n <- nrow(x)

  assert_cov_mat(A)
  assert_equal_ncol(x, A)
  # Restrict n to be greater than p?

  # Initialise penalty.
  penalty <- penalty_per_comp(b, 2, n_full, p)
  beta <- penalty$beta
  k_star <- penalty$k_star

  # Initialise B1 and B/B0.
  Q <- get_Q(x, A)
  B1 <- rep(NA, p + 1)
  B1[2] <- n * Q[1, 1] - beta
  B0 <- rep(NA, p + 1)
  B0[2] <- 0

  # Initialise J1 and J/J0.
  J1 <- rbind(rep(NA, p), diag(1, p)) # (d, d) is always 1 by def of J1
  J <- matrix(0, nrow = p + 1, ncol = p) # One {0, 1}^p-vector at [d, k].
  if (which.max(c(B0[2], B1[2])) == 2) J[2, ] <- J1[2, ]

  # Update B1, B0, J1 and J for the remaining subset sizes k.
  d <- 2  # d is d - 1 along rows, but d along columns.
  while (d <= p && sum(J[d, ]) <= k_star) {
    B0[d + 1] <- max(B0[d], B1[d], na.rm = TRUE)

    N_d <- get_neighbours_less_than(d, Q)
    alternative_savings_1 <- B0[d] + n * Q[d, d] - beta
    if (length(N_d) > 0)
      alternative_savings_1 <- c(alternative_savings_1,
                                 B1[d] + n * (Q[d, d] + 2 * sum(Q[d, N_d])) - beta)
    B1[d + 1] <- max(alternative_savings_1, na.rm = TRUE)

    update_J1_indicator <- which.max(alternative_savings_1)
    if (update_J1_indicator == 1) J1[d + 1, ] <- J[d - 1, ]
    else if(update_J1_indicator == 2) J1[d + 1, ] <- J1[d, ]
    J1[d + 1, d] <- 1  # Rows start at 0, columns at 1.

    update_J_indicator <- which.max(c(B0[d + 1], B1[d + 1]))
    if (update_J_indicator == 1) J[d + 1, ] <- J[d, ]
    else if (update_J_indicator == 2) J[d + 1, ] <- J1[d + 1, ]

    d <- d + 1
  }

  if (d < p + 1 || sum(J[p + 1, ]) > k_star) {
    B_max <- dense_savings(x, n_full, A, b)
    u_max <- rep(1, p)
  } else {
    B_max <- max(B0[p + 1], B1[p + 1], na.rm = TRUE)
    u_max <- J[p + 1, ]
  }

  list('B_max' = B_max, 'u_max' = u_max,
       'B0' = B0, 'B1' = B1, 'J1' = J1, 'J' = J, 'Q' = Q)
}

penalised_savings_ar1_NB2 <- function(x, n_full, A, b = 1) {
  ind <- function(u, d) {
    if (any(u > 1 || u < 0)) stop('u must be 0-1-vector.')
    if (length(u) < band_Q) u <- c(u, rep(0, band_Q - length(u)))

    if (length(d) <= 2) return(matrix(c(d, u + 1), ncol = band_Q + length(d)))
    else {
      # ind_mat <- matrix(0, ncol = band_Q + 2, nrow = length(d) - 1)
      # for (i in 2:length(d)) {
      #   ind_mat[i - 1, ] <- c(d[1], d[i], u + 1)
      # }
      l_d <- length(d)
      ind_mat <- matrix(rep(c(d[1], 0, u + 1), l_d - 1),
                        ncol = band_Q + 2, nrow = l_d - 1, byrow = TRUE)
      ind_mat[, 2] <- d[2:l_d]
      return(ind_mat)
    }
  }

  p <- ncol(x)
  n <- nrow(x)

  assert_cov_mat(A)
  assert_equal_ncol(x, A)
  # Restrict n to be greater than p?

  # Initialise penalty.
  penalty <- penalty_per_comp(b, 2, n_full, p)
  beta <- penalty$beta
  k_star <- penalty$k_star

  # Initialise B1 and B/B0.
  Q <- get_Q(x, A)
  B <- array(NA, dim = c(p + 1, rep(2, band(Q))))
  band_Q <- band(Q)
  B[ind(0, 2)] <- 0
  B[ind(1, 2)] <- n * Q[1, 1] - beta

  # Initialise J1 and J/J0.
  J <- array(0, dim = c(p + 1, p, rep(2, band(Q))))
  J[ind(1, c(2, 1))] <- 1
  if (which.max(c(B[ind(0, 2)], B[ind(1, 2)])) == 2)
    J[ind(0, c(2, 1:p))] <- J[ind(1, c(2, 1:p))]

  # Update B1, B0, J1 and J for the remaining subset sizes k.
  d <- 2  # d is d - 1 along rows, but d along columns.
  while (d <= p && sum(J[ind(0, c(d, 1:p))]) <= k_star) {
    B[ind(0, d + 1)] <- max(B[ind(0, d)], B[ind(1, d)], na.rm = TRUE)

    N_d <- get_neighbours_less_than(d, Q)
    alternative_savings_1 <- rep(0, 2^length(N_d))
    alternative_savings_1[1] <- B[ind(0, d)] + n * Q[d, d] - beta
    if (length(N_d) > 0) {
      for (i in 2:2^length(N_d))
        alternative_savings_1[i] <- B[ind(1, d)] + n * (Q[d, d] + 2 * sum(Q[d, N_d])) - beta
    }
    B[ind(1, d + 1)] <- max(alternative_savings_1, na.rm = TRUE)

    update_J1_indicator <- which.max(alternative_savings_1)
    if (update_J1_indicator == 1) J[ind(1, c(d + 1, 1:p))] <- J[ind(0, c(d - 1, 1:p))]
    else if(update_J1_indicator == 2) J[ind(1, c(d + 1, 1:p))] <- J[ind(1, c(d, 1:p))]
    J[ind(1, c(d + 1, d))] <- 1  # Rows start at 0, columns at 1.

    update_J_indicator <- which.max(c(B[ind(0, d + 1)], B[ind(1, d + 1)]))
    if (update_J_indicator == 1) J[ind(0, c(d + 1, 1:p))] <- J[ind(0, c(d, 1:p))]
    else if (update_J_indicator == 2) J[ind(0, c(d + 1, 1:p))] <- J[ind(1, c(d + 1, 1:p))]

    d <- d + 1
  }

  if (d < p + 1 || sum(J[ind(0, c(p + 1, 1:p))]) > k_star) {
    B_max <- dense_savings(x, n_full, A, b)
    u_max <- rep(1, p)
  } else {
    B_max <- max(B[ind(0, p + 1)], B[ind(1, p + 1)], na.rm = TRUE)
    u_max <- J[ind(0, c(p + 1, 1:p))]
  }

  list('B_max' = B_max, 'u_max' = u_max, 'B' = B, 'J' = J, 'Q' = Q)
}

penalised_savings_car_old <- function(x, n_full, A, penalty_regime = 'sparse', b = 1) {
  init_B <- function() {
    B <- array(NA, dim = c(p + 1, rep(3, band_nbs + 1)))
    B[ind(0, 2)] <- 0
    B[ind(1, 2)] <- n * Q[1, 1] - beta
    B[ind(2, 2)] <- max(B[ind(0, 2)], B[ind(1, 2)])
    # B[ind(c(0, 0), 2)] <- 0
    # B[ind(c(0, 1), 2)] <- 0
    # B[ind(c(1, 0), 2)] <- 0
    # B[ind(c(1, 1), 2)] <- n * Q[1, 1] - beta
    B
  }

  potential_savings <- function(u, nbs) {
    nbs_powerset <- subsets_from_indicators(u, nbs)
    unlist(lapply(nbs_powerset, function(nb_subset) {
      n * (Q[d, d] + 2 * sum(Q[d, nb_subset]))
    }))
  }

  p <- ncol(x)
  n <- nrow(x)

  assert_cov_mat(A)
  assert_equal_ncol(x, A)
  # Restrict n to be greater than p?

  # Initialise penalty.
  penalty <- get_penalty(penalty_regime, n_full, p, b)
  beta <- penalty$beta
  alpha <- penalty$alpha

  # Initialise Q, B and J and neibhbourshood structure.
  Q <- get_Q(x, A)
  # band_A <- band(A)
  lower_nbs <- lapply(1:p, get_neighbours_less_than, sparse_mat = A)
  all_next_nbs <- lapply(1:p, function(d) remaining_neighbours_below(d, lower_nbs, p))
  band_nbs <- band_from_list(all_next_nbs)
  ind <- init_indexing_func(band_nbs + 1)
  B <- init_B()

  maxing_steps <- c(TRUE, rep(FALSE, p - 2), TRUE)

  for (d in 2:p) {
    # print(paste0('d = ', d))
    nbs <- all_next_nbs[[d]]
    # print(nbs)
    # print(band_nbs)
    u <- extended_indicator_set(d, nbs, band_nbs)
    zero_u <- cbind(rep(0, nrow(u)), u)
    no_change_inds <- ind(zero_u, d + 1)
    one_u <- cbind(rep(1, nrow(u)), u)
    change_inds <- ind(one_u, d + 1)

    B[change_inds] <- B[ind(u, d)] + potential_savings(u, nbs) - beta
    B[no_change_inds] <- B[ind(u, d)]
    B[ind(1, d + 1)] <- max(B[change_inds])
    B[ind(0, d + 1)] <- max(B[ind(0, d)], B[ind(1, d)])

    if (d < p) {
      next_nbs <- all_next_nbs[[d + 1]]
      if (length(next_nbs) <= length(nbs)) {
        u_next <- extended_indicator_set(d + 1, next_nbs, band_nbs + 1)
        # print(u_next)
        u_extended <- sort_indicators(rbind(zero_u, one_u))
        # print(rbind(zero_u, one_u))
        # print(u_extended)
        corresponding_rows <- which_rows_correspond(u_extended, u_next)
        # print(corresponding_rows)
        for (i in 1:nrow(u_next)) {
          u_i <- u_next[i, , drop = FALSE]
          u_corresponds <- u_extended[corresponding_rows[i, ], ]
          B[ind(u_i, d + 1)] <- max(B[ind(u_corresponds, d + 1)])
        }
        maxing_steps[d] <- TRUE
      }
    }
  }

  B_max <- max(B[ind(0, p + 1)], B[ind(1, p + 1)])
  B_max <- B_max - alpha

  u_max <- backtrack_u_max(maxing_steps, B, lower_nbs)

  list('B_max' = B_max, 'u_max' = u_max, 'B' = B, 'Q' = Q)
  # list('B_max' = B_max, 'B' = B, 'Q' = Q)
}

penalised_savings_car <- function(x, n_full, A, b = 1) {
  init_B <- function(beta) {
    B <- array(NA, dim = c(p + 1, rep(3, band_nbs + 1)))
    B[ind(0, 2)] <- 0
    B[ind(1, 2)] <- n * Q[1, 1] - beta
    B[ind(2, 2)] <- max(B[ind(0, 2)], B[ind(1, 2)])
    # B[ind(c(0, 0), 2)] <- 0
    # B[ind(c(0, 1), 2)] <- 0
    # B[ind(c(1, 0), 2)] <- 0
    # B[ind(c(1, 1), 2)] <- n * Q[1, 1] - beta
    B
  }

  get_potential_savings <- function(u, nbs) {
    nbs_powerset <- subsets_from_indicators(u, nbs)
    unlist(lapply(nbs_powerset, function(nb_subset) {
      n * (Q[d, d] + 2 * sum(Q[d, nb_subset]))
    }))
  }

  p <- ncol(x)
  n <- nrow(x)

  assert_cov_mat(A)
  assert_equal_ncol(x, A)
  # Restrict n to be greater than p?

  # Initialise penalty.
  sparse_penalty <- get_penalty('sparse', n_full, p, b)
  dense_penalty <- get_penalty('dense', n_full, p, b)
  beta <- list('s' = sparse_penalty$beta, 'd' = dense_penalty$beta)
  alpha <- list('s' = sparse_penalty$alpha, 'd' = dense_penalty$alpha)

  # Initialise Q, B and J and neibhbourshood structure.
  Q <- get_Q(x, A)
  # band_A <- band(A)
  lower_nbs <- lapply(1:p, get_neighbours_less_than, sparse_mat = A)
  all_next_nbs <- lapply(1:p, function(d) remaining_neighbours_below(d, lower_nbs, p))
  band_nbs <- band_from_list(all_next_nbs)
  ind <- init_indexing_func(band_nbs + 1)
  B <- list('s' = init_B(beta$s), 'd' = init_B(beta$d))

  maxing_steps <- c(TRUE, rep(FALSE, p - 2), TRUE)
  for (d in 2:p) {
    # print(paste0('d = ', d))
    nbs <- all_next_nbs[[d]]
    # print(nbs)
    # print(band_nbs)
    u <- extended_indicator_set(d, nbs, band_nbs)
    zero_u <- cbind(rep(0, nrow(u)), u)
    no_change_inds <- ind(zero_u, d + 1)
    one_u <- cbind(rep(1, nrow(u)), u)
    change_inds <- ind(one_u, d + 1)

    potential_savings <- get_potential_savings(u, nbs)
    for (i in c('s', 'd')) {
      B[[i]][change_inds] <- B[[i]][ind(u, d)] + potential_savings - beta[[i]]
      B[[i]][no_change_inds] <- B[[i]][ind(u, d)]
      B[[i]][ind(1, d + 1)] <- max(B[[i]][change_inds])
      B[[i]][ind(0, d + 1)] <- max(B[[i]][ind(0, d)], B[[i]][ind(1, d)])
    }

    if (d < p) {
      next_nbs <- all_next_nbs[[d + 1]]
      if (length(next_nbs) <= length(nbs)) {
        u_next <- extended_indicator_set(d + 1, next_nbs, band_nbs + 1)
        # print(u_next)
        u_extended <- sort_indicators(rbind(zero_u, one_u))
        # print(rbind(zero_u, one_u))
        # print(u_extended)
        corresponding_rows <- which_rows_correspond(u_extended, u_next)
        # print(corresponding_rows)
        for (i in 1:nrow(u_next)) {
          u_i <- u_next[i, , drop = FALSE]
          u_corresponds <- u_extended[corresponding_rows[i, ], ]
          for (i in c('s', 'd')) {
            B[[i]][ind(u_i, d + 1)] <- max(B[[i]][ind(u_corresponds, d + 1)])
          }
        }
        maxing_steps[d] <- TRUE
      }
    }
  }

  B_max <- list()
  for (i in c('s', 'd')) {
    B_max[[i]] <- max(B[[i]][ind(0, p + 1)], B[[i]][ind(1, p + 1)]) - alpha[[i]]
  }
  B_max <- unlist(B_max)

  s_or_d <- c('s', 'd')[which.max(B_max)]
  u_max <- backtrack_u_max(maxing_steps, B[[s_or_d]], lower_nbs)
  B_max <- max(B_max)

  list('B_max' = B_max, 'u_max' = u_max, 'B' = B, 'Q' = Q)
  # list('B_max' = B_max, 'B' = B, 'Q' = Q)
}

penalised_savings_car_aMLE <- function(x, n_full, precision_mat_obj, b = 1) {
  init_B <- function(beta) {
    B <- array(NA, dim = c(p + 1, rep(3, band_nbs + 1)))
    B[ind(0, 2)] <- 0
    B[ind(1, 2)] <- n * (2 * w[1] - Q[1, 1]) - beta
    B[ind(2, 2)] <- max(B[ind(0, 2)], B[ind(1, 2)])
    # B[ind(c(0, 0), 2)] <- 0
    # B[ind(c(0, 1), 2)] <- 0
    # B[ind(c(1, 0), 2)] <- 0
    # B[ind(c(1, 1), 2)] <- n * Q[1, 1] - beta
    B
  }

  get_potential_savings <- function(u, nbs) {
    nbs_powerset <- subsets_from_indicators(u, nbs)
    unlist(lapply(nbs_powerset, function(nb_subset) {
      n * (2 * w[d] - Q[d, d] - 2 * sum(Q[d, nb_subset]))
    }))
  }

  p <- ncol(x)
  n <- nrow(x)
  A <- precision_mat_obj$A
  all_next_nbs <- precision_mat_obj$extended_nbs

  # Restrict n to be greater than p?

  # Initialise penalty.
  sparse_penalty <- get_penalty('sparse', n_full, p, b)
  dense_penalty <- get_penalty('dense', n_full, p, b)
  beta <- list('s' = sparse_penalty$beta, 'd' = dense_penalty$beta)
  alpha <- list('s' = sparse_penalty$alpha, 'd' = dense_penalty$alpha)

  # Initialise Q, B and J and neibhbourshood structure.
  Q <- get_Q(x, A)
  w <- get_w(x, A)
  band_nbs <- band_from_list(all_next_nbs)
  ind <- init_indexing_func(band_nbs + 1)
  B <- list('s' = init_B(beta$s), 'd' = init_B(beta$d))

  maxing_steps <- c(TRUE, rep(FALSE, p - 2), TRUE)
  for (d in 2:p) {
    # print(paste0('d = ', d))
    nbs <- all_next_nbs[[d]]
    # print(nbs)
    # print(band_nbs)
    u <- extended_indicator_set(d, nbs, band_nbs)
    zero_u <- cbind(rep(0, nrow(u)), u)
    no_change_inds <- ind(zero_u, d + 1)
    one_u <- cbind(rep(1, nrow(u)), u)
    change_inds <- ind(one_u, d + 1)

    potential_savings <- get_potential_savings(u, nbs)
    for (i in c('s', 'd')) {
      B[[i]][change_inds] <- B[[i]][ind(u, d)] + potential_savings - beta[[i]]
      B[[i]][no_change_inds] <- B[[i]][ind(u, d)]
      B[[i]][ind(1, d + 1)] <- max(B[[i]][change_inds])
      B[[i]][ind(0, d + 1)] <- max(B[[i]][ind(0, d)], B[[i]][ind(1, d)])
    }

    if (d < p) {
      next_nbs <- all_next_nbs[[d + 1]]
      if (length(next_nbs) <= length(nbs)) {
        u_next <- extended_indicator_set(d + 1, next_nbs, band_nbs + 1)
        # print(u_next)
        u_extended <- sort_indicators(rbind(zero_u, one_u))
        # print(rbind(zero_u, one_u))
        # print(u_extended)
        corresponding_rows <- which_rows_correspond(u_extended, u_next)
        # print(corresponding_rows)
        for (i in 1:nrow(u_next)) {
          u_i <- u_next[i, , drop = FALSE]
          u_corresponds <- u_extended[corresponding_rows[i, ], ]
          for (i in c('s', 'd')) {
            B[[i]][ind(u_i, d + 1)] <- max(B[[i]][ind(u_corresponds, d + 1)])
          }
        }
        maxing_steps[d] <- TRUE
      }
    }
  }

  B_max <- list()
  for (i in c('s', 'd')) {
    B_max[[i]] <- max(B[[i]][ind(0, p + 1)], B[[i]][ind(1, p + 1)]) - alpha[[i]]
  }
  B_max <- unlist(B_max)

  s_or_d <- c('s', 'd')[which.max(B_max)]
  u_max <- backtrack_u_max(maxing_steps, B[[s_or_d]], all_next_nbs)
  J_max <- (1:p)[as.logical(u_max)]
  B_max <- max(B_max)

  list('B_max' = B_max, 'u_max' = u_max, 'J_max' = J_max, 'B' = B, 'Q' = Q)
  # list('B_max' = B_max, 'B' = B, 'Q' = Q)
}

penalised_savings_ar <- function(x, n_full, A, penalty_regime = 'sparse', b = 1) {
  init_B <- function() {
    B <- array(NA, dim = c(p + 1, rep(3, band_A + 1)))
    B[ind(0, 2)] <- 0
    B[ind(1, 2)] <- n * Q[1, 1] - beta
    # B[ind(c(0, 0), 2)] <- 0
    # B[ind(c(0, 1), 2)] <- 0
    # B[ind(c(1, 0), 2)] <- 0
    # B[ind(c(1, 1), 2)] <- n * Q[1, 1] - beta
    B
  }

  potential_savings <- function(u, nbs) {
    nbs_powerset <- subsets_from_indicators(u, nbs)
    unlist(lapply(nbs_powerset, function(nb_subset) {
      n * (Q[d, d] + 2 * sum(Q[d, nb_subset]))
    }))
  }

  p <- ncol(x)
  n <- nrow(x)

  assert_cov_mat(A)
  assert_equal_ncol(x, A)
  # Restrict n to be greater than p?

  # Initialise penalty.
  penalty <- get_penalty(penalty_regime, n_full, p, b)
  beta <- penalty$beta
  alpha <- penalty$alpha

  # Initialise Q, B and J and neibhbourshood structure.
  Q <- get_Q(x, A)
  band_A <- band(A)
  ind <- init_indexing_func(band_A + 1)
  B <- init_B()
  lower_nbs <- lapply(1:p, get_neighbours_less_than, sparse_mat = A)
  for (d in 2:p) {
    # print(paste0('d = ', d))
    u <- indicator_power_set(length(lower_nbs[[d]]))
    zero_u <- cbind(rep(0, nrow(u)), u)
    no_change_inds <- ind(zero_u, d + 1)
    one_u <- cbind(rep(1, nrow(u)), u)
    change_inds <- ind(one_u, d + 1)

    B[change_inds] <- B[ind(u, d)] + potential_savings(u, lower_nbs[[d]]) - beta
    B[no_change_inds] <- B[ind(u, d)]
    B[ind(1, d + 1)] <- max(B[change_inds])
    B[ind(0, d + 1)] <- max(B[ind(0, d)], B[ind(1, d)])

    u0 <- cbind(u, rep(0, nrow(u)))
    u1 <- cbind(u, rep(1, nrow(u)))
    B[ind(u, d + 1)] <- pmax(B[ind(u0, d + 1)], B[ind(u1, d + 1)])
  }

  B_max <- max(B[ind(0, p + 1)], B[ind(1, p + 1)]) - alpha
  maxing_inds <- c(1, band_A + 1, (band_A + 2):p)
  maxing_steps <- rep(FALSE, p)
  maxing_steps[maxing_inds] <- TRUE
  u_max <- backtrack_u_max(maxing_steps, B, lower_nbs)

  list('B_max' = B_max, 'u_max' = u_max, 'B' = B, 'Q' = Q)
}

penalised_savings_BF_old <- function(x, n_full, A, penalty_regime = 'sparse', b = 1) {
  p <- ncol(x)
  n <- nrow(x)

  assert_cov_mat(A)
  assert_equal_ncol(x, A)
  # Restrict n to be greater than p?

  # Initialise penalty.
  penalty <- get_penalty(penalty_regime, n_full, p, b)
  beta <- penalty$beta
  alpha <- penalty$alpha

  Q <- get_Q(x, A)
  u_mat <- indicator_power_set(p)
  B <- rep(0, nrow(u_mat))
  for (i in 1:nrow(u_mat)) {
    u <- t(u_mat[i, ])
    B[i] <- n * u %*% Q %*% t(u) - sum(u) * beta
  }

  B_max <- max(B) - alpha
  B_which_max <- which.max(B)
  u_max <- u_mat[B_which_max, ]

  list('B_max' = B_max, 'u_max' = u_max, 'B' = B, 'u' = u_mat)
}

penalised_savings_BF <- function(x, n_full, A, b = 1) {
  p <- ncol(x)
  n <- nrow(x)

  assert_cov_mat(A)
  assert_equal_ncol(x, A)
  # Restrict n to be greater than p?

  # Initialise penalty.
  penalty_vec <- combined_penalty_vec(n_full, p, b)

  Q <- get_Q(x, A)
  u_mat <- indicator_power_set(p)[-1, ]
  B <- rep(0, nrow(u_mat))
  for (i in 1:nrow(u_mat)) {
    u <- t(u_mat[i, ])
    B[i] <- n * u %*% Q %*% t(u) - penalty_vec[sum(u)]
  }

  B_max <- max(B)
  B_which_max <- which.max(B)
  u_max <- u_mat[B_which_max, ]

  list('B_max' = B_max, 'u_max' = u_max, 'B' = B, 'u' = u_mat)
}

optim_penalised_savings_c_old <- function(x, n_full, precision_mat_obj, b = 1, vec = FALSE) {
  optim_func <- function(penalty, vec) {
    if (vec) return(optimise_savings_vec(Q, w, penalty, precision_mat_obj$extended_nbs))
    else return(optimise_savings(Q, w, penalty, precision_mat_obj$extended_nbs))
  }

  n <- nrow(x)
  p <- ncol(x)
  Q <- n * get_Q(x, precision_mat_obj$A)
  w <- n * get_w(x, precision_mat_obj$A)

  penalty_regimes <- c('sparse', 'dense')
  res <- lapply(penalty_regimes, function(regime) {
    penalty <- get_penalty(regime, n_full, p, b)
    optim_func(penalty, vec)
  })
  both_B_max <- lapply(res, function(res_i) res_i$B_max)
  res[[which.max(both_B_max)]]
}

# dense_savings <- function(x, n_full, A, b = 1) {
#   p <- ncol(x)
#   n <- nrow(x)
#   Q <- get_Q(x, A)
#   u <- rep(1, p)
#   beta <- const_penalty(b, 2, n_full, p)
#   n * u %*% Q %*% u - beta
# }

subset_from_indicator <- function(u, set) {
  u <- u[u == 0 | u == 1]
  if (length(u) != length(set))
    stop(' the number of 0s and 1s in u and set must be of equal length')
  set[as.logical(u)]
}

subsets_from_indicators <- function(u, set) {
  lapply(split(u, seq(nrow(u))), subset_from_indicator, set)
}

band_from_list <- function(nbs_list) {
  p <- length(nbs_list)
  max(unlist(lapply(2:p, function(d) {
    if (length(nbs_list[[d]]) >= 1) return(d - min(nbs_list[[d]]))
    else return(0)
  })))
}

extended_indicator_set <- function(d, nbs, n_col) {
  l_nbs <- length(nbs)
  u_mat <- indicator_power_set(l_nbs)
  if (is.null(ncol(u_mat))) return(matrix(2, nrow = 1, ncol = n_col))
  else if (ncol(u_mat) == n_col) return(u_mat)
  else {
    extended_u_mat <- matrix(2, ncol = n_col, nrow = nrow(u_mat))
    extended_u_mat[, d - nbs] <- u_mat
    return(extended_u_mat)
  }
}

expand_indicator <- function(u_root, d, nbs) {
  if (sum(u_root == 2) == 0) return(u_root)
  else {
    expand_which = intersect(which(u_root == 2), c(1, d - nbs + 1))
    # if (u_root[1] == 2 && !any(expand_which == 1))
    #   expand_which <- c(1, expand_which)
    ind_power_set <- indicator_power_set(length(expand_which))
    m <- nrow(ind_power_set)
    u_expanded <- matrix(rep(u_root, m), nrow = m, byrow = TRUE)
    u_expanded[, expand_which] <- ind_power_set
    return(u_expanded)
  }
}

count_leading <- function(number, x) {
  i <- 1
  n <- 0
  while (x[i] == number) {
    n <- n + 1
    i <- i + 1
  }
  n
}

which_rows_correspond <- function(u_mat_big, u_mat) {

  if(!is.matrix(u_mat)) u_mat <- matrix(u_mat, ncol = length(u_mat))

  if (nrow(u_mat) == 1) {
    m <- nrow(u_mat_big)
    return(matrix(1:m, nrow = 1, ncol = m))
  } else {
    col2 <- which(u_mat[1, ] == 2)
    u_mat_big2 <- u_mat_big[, - col2, drop = FALSE]
    u_mat2 <- u_mat[, - col2, drop = FALSE]
    n <- nrow(u_mat2)
    m <- nrow(u_mat_big2)
    # print(n)
    # print(m)
    # print(m / n)
    which_mat <- matrix(0, nrow = n, ncol = m / n)
    for (j in 1:m) {
      k <- n
      candidates <- 1:k
      for (i in ncol(u_mat2):1) {
        if (u_mat_big2[j, i] == 0) candidates <- candidates[1:(k/2)]
        else                       candidates <- candidates[(k/2 + 1):k]
        k <- k/2
      }
      which_mat[candidates, which.min(which_mat[candidates, ])] <- j
    }
    return(which_mat)
  }
}

sort_indicators <- function(u_mat) {
  column_order <- function(x) {
    order(apply(x, 2, function(v) count_leading(0, v)))
  }
  res <- matrix(2, nrow = nrow(u_mat), ncol = ncol(u_mat))
  which_col01 <- which(u_mat[1, ] != 2)
  u_mat01 <- u_mat[, which_col01]
  res[, which_col01] <- u_mat01[, column_order(u_mat01)]
  res
}

backtrack_u_max <- function(maxing_steps, B, extended_nbs) {
  next_wmax <- function(d, d_prev, v) {
    nbs <- extended_nbs[[d]]
    u_length <- d - min(nbs) + 1
    u_root <- c(v, rep(2, u_length - length(v)))
    u_expanded <- expand_indicator(u_root, d, nbs)
    u_max <- u_expanded[which.max(B[ind(u_expanded, d + 1)]), ]
    if (d - d_prev < u_length) v_next <- u_max[(d - d_prev + 1):u_length]
    else v_next <- integer(0)
    return(list(u_max[1:(d - d_prev)], v_next))
  }

  combine <- function(u_list) {
    u_list <- lapply(u_list, rev)
    u_wmax_vec <- do.call('c', u_list)
    if (!any(u_wmax_vec == 2)) return(u_wmax_vec)
    else {
      u_vec <- u_list[[1]]
      for (i in 2:length(u_list)) {
        u_part <- u_list[[i]]
        if (any(u_vec == 2)) {
          which_2 <- which(u_vec == 2)
          l_2 <- length(which_2)
          u_vec[which_2] <- u_part[1:l_2]
          u_vec <- c(u_vec, u_part[(l_2 + 1):length(u_part)])
        } else {
          u_vec <- c(u_vec, u_part)
        }
      }
      return(u_vec)
    }
  }

  p <- length(maxing_steps)
  d_max <- which(maxing_steps)
  # print(d_max)
  ind <- init_indexing_func(length(dim(B)) - 1)
  v_next <- which.max(c(B[ind(0, p + 1)], B[ind(1, p + 1)])) - 1
  u_max_list <- list()
  for (i in length(d_max):2) {
    # print(paste0('d_max = ', d_max[i]))
    res <- next_wmax(d_max[i], d_max[i - 1], v_next)
    # print(res)
    v_next <- res[[2]]
    u_max_list[[i]] <- res[[1]]
  }
  u_max_list[[1]] <- v_next[length(v_next)]
  combine(u_max_list)
}

init_indexing_func <- function(n_col) {
  function(u, d) {
    if (!is.matrix(u)) u <- matrix(u, ncol = length(u))
    if (any(u > 2 || u < 0)) stop('u must be 0-1-2-matrix or vector.')
    if (ncol(u) < n_col)
      u <- cbind(u, matrix(2, ncol = n_col - ncol(u), nrow = nrow(u)))
    cbind(rep(d, nrow(u)), u + 1)
  }
}

