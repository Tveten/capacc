robust_cov_mat <- function(x) {
  p <- ncol(x)
  n <- nrow(x)
  scores_x <- apply(x, 2, function(x_i) qnorm(rank(x_i)/(n + 1)))
  mad_x <- diag(apply(x, 2, mad))
  mad_x %*% cor(scores_x) %*% mad_x
}

zero_inds <- function(adj_mat) {
  p <- nrow(adj_mat)
  zeros <- list()
  for (i in 2:p) {
    for (j in 1:(i - 1)) {
      if (adj_mat[i, j] >= -1e-8 && adj_mat[i, j] <= 1e-8)
        zeros[[length(zeros) + 1]] <- c(i, j)
    }
  }
  do.call("rbind", zeros)
}

estimate_precision_mat <- function(x, adj_mat, sparse = TRUE) {
  if (ncol(x) >= nrow(x))
    warning("There may be convergence issues when p >= n.")

  S <- robust_cov_mat(x)
  zeros <- zero_inds(adj_mat)
  S_inv <- suppressWarnings(glasso::glasso(S, rho = 0, zero = zeros)$wi)
  if (sparse) return(Matrix::Matrix(S_inv, sparse = TRUE))
  else return(S_inv)
}

CAR_func <- function(W) {
  p <- nrow(W)
  function(sigma2, rho) {
    1 / sigma2 * (diag(1, p) - rho * W)
  }
}

CAR_with_nugget_func <- function(W) {
  p <- nrow(W)
  function(sigma2, rho, delta) {
    solve(sigma2 * solve((diag(1, p) - rho * W)) + diag(delta, p))
  }
}

estimate_CAR_precision_mat <- function(x, adj_mat, Q_theta = NULL, param_lim = NULL) {
  obj_func <- function(...) {
    log(det(Q_theta(...))) - sum(diag((S %*% Q_theta(...))))
  }

  n <- nrow(x)
  p <- ncol(x)
  S <- robust_cov_mat(x)
  eigen_values_adj_mat <- eigen(adj_mat, symmetric = TRUE, only.values = TRUE)$values
  Q_theta <- CAR_with_nugget_func(adj_mat)

  param_lim <- list("sigma2" = c(0, 1.3),
                    "rho"    = c(1 / eigen_values_adj_mat[p], 1 / eigen_values_adj_mat[1]),
                    "delta"  = c(0, 1.3))
  param_grid <- expand.grid(lapply(param_lim, function(lim) {
    seq(lim[1] + 0.01, lim[2] - 0.01, length.out = 20)
  }))
  param_grid <- param_grid[param_grid$sigma2 + param_grid$delta >= 0.7, ]
  obj_func_evals <- unlist(Map(obj_func, sigma2 = param_grid$sigma2, rho = param_grid$rho, delta = param_grid$delta))
  optim_params <- param_grid[which.max(obj_func_evals), ]
  optim_Q_theta <- Q_theta(optim_params$sigma2, optim_params$rho, optim_params$delta)
  list("Q" = optim_Q_theta, "params" = optim_params, "value" = max(obj_func_evals))
}


#==============#
##### OLD ######
#==============#
BIC_precision_mat <- function(A, S, n, rho) {
  nnz <- sum(A[lower.tri(A, diag = TRUE)] != 0)
  - log(det(A)) + sum(diag(A %*% S)) + log(n) / n * nnz
}

select_precision_mat <- function(glasso_res, S, n) {
  BIC_per_lambda <- unlist(lapply(1:length(glasso_res$lambda), function(i) {
    BIC_precision_mat(glasso_res$icov[[i]], S, n, glasso_res$lambda[i])
  }))
  glasso_res$icov[[which.min(BIC_per_lambda)]]
}

glasso_precision_mat <- function(lambda_min_ratio = 0.1, scr = FALSE, scr_num = 10, eps = 0.2) {
  if (precision_mat == 'automatic')
    precision_mat <- estimate_sparse_precision_mat(x, lambda_min_ratio, eps = eps,
                                                   scr = scr, scr_num = scr_num,
                                                   standardise = FALSE)
  precision_mat_order <- 1:p
  band_precision_mat <- band(precision_mat)
  reordered_precision_mat <- cuthill_mckee(precision_mat, return_all = TRUE)
  if (reordered_precision_mat$band < band_precision_mat) {
    # print('Reordering')
    precision_mat <- reordered_precision_mat$x
    precision_mat_order <- reordered_precision_mat$order
    band_precision_mat <- reordered_precision_mat$band
    x <- x[, precision_mat_order]
  }
  print('Estimated precision matrix:')
  print(precision_mat)
  print(paste0('Band = ', band(precision_mat)))
  if (band_precision_mat > 15)
    warning(paste0('Computation will be slow because the band of the precision matrix is ', band_precision_mat))
  if (band_precision_mat >= 20)
    stop(paste0('The band of the precision matrix is too high for computation: ', band_precision_mat))
  precision_mat
}

estimate_sparse_precision_mat <- function(x, lambda_min_ratio = 0.1, nlambda = 10,
                                          scr = FALSE, scr_num = 10, eps = 0.2,
                                          robust = TRUE, standardise = TRUE) {
  if (standardise) x <- anomaly::robustscale(x)
  if (robust) S <- robust_cov_mat(x)
  else        S <- var(x)
  glasso_res <- huge::huge(S, lambda.min.ratio = lambda_min_ratio,
                           nlambda = nlambda, scr = scr, scr.num = scr_num,
                           method = 'glasso')
  prec_mat <- select_precision_mat(glasso_res, S, nrow(x))
  prec_mat[abs(prec_mat) <= eps] <- 0
  prec_mat
}
