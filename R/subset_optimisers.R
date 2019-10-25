####
#### Helper function. ----------------------------------------------------------
####

dense_savings <- function(x, n_full, A, b = 1) {
  p <- ncol(x)
  n <- nrow(x)
  Q <- get_Q(x, A)
  u <- rep(1, p)
  beta <- const_penalty(b, 2, n_full, p)
  n * u %*% Q %*% u - beta
}

linear_penalty <- function(n, p, b = 1,  a = 2, C = 2) {
  a <- as.numeric(a)
  psi <- a * log(n)
  k_star <- (p + C * sqrt(p * psi)) / (2 * log(p))
  list('beta' = b * (C * psi / k_star + C * log(p)), 'k_star' = k_star)
}

const_penalty <- function(n, p, b = 1, a = 2, C = 2) {
  a <- as.numeric(a)
  psi <- a * log(n)
  b * (p + C * psi + C * sqrt(p * psi))
}

get_Q <- function(x, A) {
  mean_x <- colMeans(x, na.rm = TRUE)  # Row vector
  t(t(mean_x)) %*% t(mean_x) * A
}

get_neighbours_less_than <- function(i, sparse_mat) {
  nbs <- which(!is_zero(sparse_mat[i, ]))
  rev(nbs[nbs < i])
}

band <- function(x) {
  nnz_lower_diag <- rep(0, nrow(x) - 1)
  for (i in 2:nrow(x)) {
    nnz_lower_diag[i] <- sum(!is_zero(x[i, 1:(i - 1)]))
  }
  max(nnz_lower_diag)
}

indicator_power_set <- function(size, as_list = FALSE) {
  set <- 1:size
  ind_power_set <- vapply(2^(1:size - 1), function(k) {
    rep(c(0, 1), each = k, times = 2^size / (2 * k))
  }, numeric(2^size))
  if (as_list) ind_power_set <- split(ind_power_set, seq(nrow(ind_power_set)))
  ind_power_set
}

subset_from_indicator <- function(u, set) {
  u <- u[u == 0 | u == 1]
  if (length(u) != length(set))
    stop(' the number of 0s and 1s in u and set must be of equal length')
  set[as.logical(u)]
}

subsets_from_indicators <- function(u, set) {
  lapply(split(u, seq(nrow(u))), subset_from_indicator, set)
}

extended_indicator_set <- function(d, nbs, n_col) {
  l_nbs <- length(nbs)
  u_mat <- indicator_power_set(l_nbs)
  if (ncol(u_mat) == n_col) return(u_mat)
  else {
    extended_u_mat <- matrix(2, ncol = n_col, nrow = nrow(u_mat))
    extended_u_mat[, d - nbs] <- u_mat
    return(extended_u_mat)
  }
}

expand_indicator <- function(u_root) {
  if (sum(u_root == 2) == 0) return(u_root)
  else {
    which_is_2 <- which(u_root == 2)
    ind_power_set <- indicator_power_set(length(which_is_2))
    m <- nrow(ind_power_set)
    u_expanded <- matrix(rep(u_root, m), nrow = m, byrow = TRUE)
    u_expanded[, which_is_2] <- ind_power_set
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

  col2 <- which(u_mat[1, ] == 2)
  u_mat_big2 <- u_mat_big[, - col2, drop = FALSE]
  u_mat2 <- u_mat[, - col2, drop = FALSE]
  n <- nrow(u_mat2)
  m <- nrow(u_mat_big2)
  which_mat <- matrix(0, nrow = n, ncol = 2 * log(m/n, base = 2))
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
  which_mat
}

sort_indicators <- function(u_mat) {
  column_order <- function(x) {
    order(apply(x, 2, function(v) count_leading(0, v)))
  }
  res <- matrix(2, nrow = nrow(u_mat), ncol = ncol(u_mat))
  which_col01 <- which(u_mat[1, ] != 2)
  res[, which_col01] <- u_mat[, column_order(u_mat[, which_col01])]
  res
}

get_penalty <- function(penalty_regime, n, p, b = 1, a = 2) {
  if (penalty_regime == 'sparse') {
    beta <- linear_penalty(n, p, b = b, a = a)$beta
    alpha <- 0
  } else if (penalty_regime == 'dense') {
    beta <- 0
    alpha <- const_penalty(n, p, b = b, a = a)
  }
  list('alpha' = alpha, 'beta' = beta)
}

remaining_neighbours_below <- function(d, nb_list, p) {
  rev(intersect(1:(d - 1), unlist(nb_list[d:p])))
}

backtrack_J_max <- function(maxing_steps, B, lower_nbs) {
  next_wmax <- function(d, d_prev, v) {
    nbs <- remaining_neighbours_below(d, lower_nbs, p)
    u_length <- length(nbs) + 1
    u_root <- c(v, rep(2, u_length - length(v)))

    u_expanded <- expand_indicator(u_root)
    u_max <- u_expanded[which.max(B[ind(u_expanded, d + 1)]), ]
    v_next <- u_max[(d - d_prev + 1):u_length]
    return(list(u_max[1:(d - d_prev)], v_next))
  }

  combine <- function(u_list) {
    u_list <- lapply(u_list, rev)
    u_wmax_vec <- do.call('c', u_list)
    if (!any(u_wmax_vec == 2)) return(u_wmax_vec)
    else {

    }
  }

  p <- length(maxing_steps)
  d_max <- which(maxing_steps)
  ind <- init_ind(length(dim(B)) - 1)
  v_next <- which.max(c(B[ind(0, p + 1)], B[ind(1, p + 1)])) - 1
  u_max_list <- list()
  for (i in length(d_max):2) {
    res <- next_wmax(d_max[i], d_max[i - 1], v_next)
    v_next <- res[[2]]
    u_max_list[[i]] <- res[[1]]
  }
  u_max_list[[1]] <- v_next[length(v_next)]
  combine(u_max_list)
}

init_ind <- function(n_col) {
  function(u, d) {
    if (!is.matrix(u)) u <- matrix(u, ncol = length(u))
    if (any(u > 1 || u < 0)) stop('u must be 0-1-matrix or vector.')
    if (ncol(u) < n_col)
      u <- cbind(u, matrix(2, ncol = n_col - ncol(u), nrow = nrow(u)))
    cbind(rep(d, nrow(u)), u + 1)
  }
}

####
#### Penalized savings optimisers. ---------------------------------------------
####

penalised_savings_car <- function(x, n_full, A, penalty_regime = 'sparse', b = 1) {
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
  ind <- init_ind(band_A + 1)
  B <- init_B()

  lower_nbs <- lapply(1:p, get_neighbours_less_than, sparse_mat = A)
  maxing_steps <- c(TRUE, rep(FALSE, p - 2), TRUE)

  for (d in 2:p) {
    # print(paste0('d = ', d))
    nbs <- remaining_neighbours_below(d, lower_nbs, p)
    u <- extended_indicator_set(d, nbs, band_A)
    zero_u <- cbind(rep(0, nrow(u)), u)
    no_change_inds <- ind(zero_u, d + 1)
    one_u <- cbind(rep(1, nrow(u)), u)
    change_inds <- ind(one_u, d + 1)

    B[change_inds] <- B[ind(u, d)] + potential_savings(u, nbs) - beta
    B[no_change_inds] <- B[ind(u, d)]
    B[ind(1, d + 1)] <- max(B[change_inds])
    B[ind(0, d + 1)] <- max(B[ind(0, d)], B[ind(1, d)])

    next_nbs <- remaining_neighbours_below(d + 1, lower_nbs, p)
    if (length(next_nbs) <= length(nbs) && d < p) {
      u_next <- extended_indicator_set(d + 1, next_nbs, band_A + 1)
      u_extended <- sort_indicators(rbind(zero_u, one_u))
      corresponding_rows <- which_rows_correspond(u_extended, u_next)
      for (i in 1:nrow(u_next)) {
        u_i <- u_next[i, , drop = FALSE]
        u_corresponds <- u_extended[corresponding_rows[i, ], ]
        B[ind(u_i, d + 1)] <- max(B[ind(u_corresponds, d + 1)])
      }
      maxing_steps[d] <- TRUE
    }
  }

  B_max <- max(B[ind(0, p + 1)], B[ind(1, p + 1)])
  # J_max <-  J[p + 1, ]
  if (penalty_regime == 'dense') B_max <- B_max - alpha

  J_max <- backtrack_J_max(maxing_steps, B, lower_nbs)

  list('B_max' = B_max, 'J_max' = J_max, 'B' = B, 'Q' = Q)
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
  ind <- init_ind(band_A + 1)
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

  B_max <- max(B[ind(0, p + 1)], B[ind(1, p + 1)])
  if (penalty_regime == 'dense') B_max <- B_max - alpha
  maxing_inds <- c(1, band_A + 1, (band_A + 2):p)
  maxing_steps <- rep(FALSE, p)
  maxing_steps[maxing_inds] <- TRUE
  J_max <- backtrack_J_max(maxing_steps, B, lower_nbs)

  list('B_max' = B_max, 'J_max' = J_max, 'B' = B, 'Q' = Q)
}

penalised_savings_BF <- function(x, n_full, A, penalty_regime = 'sparse', b = 1) {
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
  B <- rep(0, nrow(u_mat) + 1)
  for (i in 1:nrow(u_mat)) {
    u <- t(u_mat[i, ])
    B[i] <- n * u %*% Q %*% t(u) - sum(u) * beta
  }

  B_max <- max(B)
  if (penalty_regime == 'dense') B_max <- B_max - alpha
  B_which_max <- which.max(B)
  J_max <- u_mat[B_which_max, ]

  list('B_max' = B_max, 'J_max' = J_max, 'B' = B, 'u' = u_mat)
}

####
#### Testing functions ------------------------------------------
####

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
  else return(list('BF' = optim_BF$J_max, 'OP' = optim_CAR_OP$J_max))
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

B_compare <- function(method = 'ar', p = 10, r = 2, rho = 0.9, prop = 0.8) {
  ind <- function(u, d) {
    if (any(u > 1 || u < 0)) stop('u must be 0-1-vector.')
    if (length(u) < band_A + 1) u <- c(u, rep(2, band_A + 1 - length(u)))

    matrix(c(d, u + 1), ncol = band_A + 1 + length(d))
  }

  set.seed(1)
  res <- compare_OP_BF(method, p = p, r = r, rho = rho, prop = prop, return_all = TRUE)
  B_BF <- res$BF$B
  u <- res$BF$u
  B_OP <- res$OP$B
  B_mat <- matrix(0, nrow = nrow(u), ncol = 2)
  # print(B_OP)
  dim_B <- dim(B_OP)

  print(c(res$BF$B_max, res$OP$B_max))
  print(res$BF$J_max)
  print(res$OP$J_max)
  for (i in 1:dim_B[1])

  u_OP <- function(B_OP) {

  }

  B_OP_corresponding_to <- function(u) {

  }
}

# test_that('CAR OP returns equally as brute force', {
#   n <- 100
#   p <- 10
#   r <- 1:4
#   set.seed(5)
#   rho <- c(-0.4, -0.2, 0.2, 0.5, 0.9)
#   prop <- c(0.2, 0.5, 0.8, 1)
#   res_B <- array(0, dim = c(length(rho), length(prop), length(r)))
#   res_J <- array(FALSE, dim = c(length(rho), length(prop), length(r)))
#   res <- list()
#   for (k in seq_along(r)) {
#     for (j in seq_along(prop)) {
#       for (i in seq_along(rho)) {
#         res <- compare_penalised_savings_AR(n, p, r[k], rho[i], prop[j])
#         expect_equal(res$OP$B_max, res$BF$B_max)
#         expect_equal(res$OP$J_max, res$BF$J_max)
#       }
#     }
#   }
# })
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

  list('k_max' = optimal_k, 'J_max' = optimal_J, 'B_max' = optimal_savings,
       'S_obj' = optimal_subset_res$S, 'B0' = optimal_subset_res$B0,
       'B1' = optimal_subset_res$B1)
}

penalised_savings_iid <- function(x_subset, n, A = diag(1, ncol(x_subset)), b = 1) {
  n_subset <- nrow(x_subset)
  p <- ncol(x_subset)
  beta <- b * adjusted_penalty(2, n, p)

  mean_x2 <- colMeans(x_subset, na.rm = TRUE)^2
  order_mean_x2 <- order(mean_x2, decreasing = TRUE)
  sorted_mean_x2 <- mean_x2[order_mean_x2]
  savings_k <- cumsum(n_subset * sorted_mean_x2 - beta)
  optimal_savings <- max(savings_k)
  optimal_k <- which.max(savings_k)
  optimal_J <- order_mean_x2[1:optimal_k]

  return(list('k_max' = optimal_k, 'J_max' = optimal_J, 'S_max' = optimal_savings))
}

#' @param x A data matrix with p columns and n = (e - s) rows.
#' @param A A precision matrix.
#' @param l The maximum size of the subset.
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
    J_max <- rep(1, p)
  } else {
    B_max <- max(B0[p + 1], B1[p + 1], na.rm = TRUE)
    J_max <- J[p + 1, ]
  }

  list('B_max' = B_max, 'J_max' = J_max,
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
    J_max <- rep(1, p)
  } else {
    B_max <- max(B0[p + 1], B1[p + 1], na.rm = TRUE)
    J_max <- J[p + 1, ]
  }

  list('B_max' = B_max, 'J_max' = J_max,
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
    J_max <- rep(1, p)
  } else {
    B_max <- max(B[ind(0, p + 1)], B[ind(1, p + 1)], na.rm = TRUE)
    J_max <- J[ind(0, c(p + 1, 1:p))]
  }

  list('B_max' = B_max, 'J_max' = J_max, 'B' = B, 'J' = J, 'Q' = Q)
}
