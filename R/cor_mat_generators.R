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

constant_cor_mat2 <- function(p, rho, standardise = TRUE, sparse = TRUE) {
  one_vec <- matrix(1, nrow = p)
  Sigma <- (1 - rho) * diag(1, p) + rho / p * one_vec %*% t(one_vec)
  print(Sigma)
  Q <- solve(Sigma)

  if (standardise) Q <- standardise_precision_mat(Q)
  if (sparse) return(Matrix::Matrix(Q, sparse = TRUE))
  else return(Q)
}

ar_precision_mat <- function(p, phi, sparse = TRUE) {
  Sigma_inv <- matrix(0, nrow = p, ncol = p)
  diag(Sigma_inv) <- 1 + phi^2
  diag(Sigma_inv)[c(1, p)] <- 1

  delta <- row(Sigma_inv) - col(Sigma_inv)
  Sigma_inv[off_diag_ind(c(-1, 1), p)] <- -phi

  Sigma_inv <- Sigma_inv / (1 - phi^2)
  if (sparse) return(Matrix::Matrix(Sigma_inv, sparse = TRUE))
  else return(Sigma_inv)
}

#' @export
car_precision_mat <- function(nbs = banded_neighbours(2, 10), rho = 0.5, sigma = 1,
                              min_nbs = 1, max_nbs = 3, standardised = TRUE, sparse = TRUE) {
  p <- length(nbs)
  W <- adjacency_mat(nbs, sparse)
  if (all(unlist(Map(length, nbs)) == 0)) {
    Sigma_inv <- diag(1, p)
  } else {
    eigen_values_W <- eigen(W)$values
    if (rho <= 1 / eigen_values_W[p]) {
      min_rho <- 1 / eigen_values_W[p]
      warning(paste0('rho has been altered to 1 / smallest_eigen_value(W) = ', min_rho))
      rho <- min_rho + sqrt(.Machine$double.eps)
    }
    Sigma_inv <- 1 / sigma^2 * (diag(as.vector(W %*% rep(1, p))) - rho * W)
    if (standardised) Sigma_inv <- standardise_precision_mat(Sigma_inv)
  }

  if (sparse) return(Matrix::Matrix(Sigma_inv, sparse = TRUE))
  else return(Sigma_inv)
}

#' @export
adjacency_mat <- function(nbs_list, sparse = TRUE) {
  p <- length(nbs_list)
  adj_mat <- do.call('rbind', lapply(nbs_list, indicator, p = p))
  if (!isSymmetric(adj_mat)) {
    adj_mat <- symmetric_from_lower_tri(adj_mat)
  }
  if (sparse) return(Matrix::Matrix(adj_mat, sparse = TRUE))
  else return(adj_mat)
}

symmetric_from_lower_tri <- function(x_tri) {
  upper_triangle <- upper.tri(x_tri)
  x_tri[upper_triangle] <- t(x_tri)[upper_triangle]
  x_tri
}

#' @export
banded_neighbours <- function(band, p) {
  if (band == 0) {
    nbs <- lapply(1:p, function(i) integer(0))
  } else {
    nbs <- list()
    for (i in 1:p) {
      if (i == 1) lower <- integer(0)
      else        lower <- max(1, i - band):(i - 1)
      if (i == p) upper <- integer(0)
      else        upper <- (i + 1):min(i + band, p)
      nbs[[i]] <- c(lower, upper)
    }
  }
  nbs
}

number_neighbours <- function(n_vec_lower) {
  p <- length(n_vec_lower) + 1
  nbs <- list(integer(0))
  for (i in 2:p) {
    if (n_vec_lower[i - 1] == 0) nbs[[i]] <- integer(0)
    else nbs[[i]] <- max(1, i - n_vec_lower[i - 1]):(i - 1)
  }
  nbs
}

# random_neighbours_banded <- function(p, min_nbs = 1, max_nbs = 3) {
#   max_nbs <- min(max_nbs, p - 1)
#   min_nbs <- min(max_nbs, max(1, min_nbs))
#
#   n_nbs <- vapply(1:(p - 1), function(i) sample(min_nbs:max_nbs, 1), numeric(1))
#   number_neighbours(n_nbs)
# }

random_neighbours <- function(p, min_nbs = 1, max_nbs = 3) {
  adj_mat <- radjacency_mat(p, min_nbs, max_nbs)
  nbs <- list()
  for (i in 1:p) {
    nbs[[i]] <- which(adj_mat[i, ] != 0)
  }
  nbs
}

#' @export
lattice_neighbours <- function(p) {
  grid_nbs <- function(i, j, n_row, n_col) {
    grid_ind <- c(i, j - 1,
                      i - 1, j)
    if (j < n_col) grid_ind <- c(grid_ind, i, j + 1)
    if (i < n_row) grid_ind <- c(grid_ind, i + 1, j)
    n_nbs <- length(grid_ind) / 2
    matrix(grid_ind, nrow = n_nbs, ncol = 2, byrow = TRUE)
  }

  width <- ceiling(sqrt(p))
  if (width * (width - 1) <= p) height <- width
  else height <- width - 1

  nbs <- list()
  grid <- matrix(c(1:p, rep(0, height * width - p)),
                 nrow = height, ncol = width, byrow = TRUE)
  for (i in 1:height) {
    for (j in 1:width) {
      if (grid[i, j] != 0) {
        curr_nbs <- grid[grid_nbs(i, j, height, width)]
        curr_nbs <- curr_nbs[curr_nbs != 0]
        nbs[[(i - 1) * width + j]] <- curr_nbs
      }
    }
  }
  nbs
}

radjacency_mat <- function(p, min_nbs = 1, max_nbs = 3) {
  get_possible_nbs <- function(i) {
    possible_nbs <- setdiff(1:p, i)
    which_already_nbs <- which(adj_mat[i, ] == 1)
    possible_nbs <- setdiff(possible_nbs, which_already_nbs)
    n_nbs_already <- as.vector(apply(adj_mat[possible_nbs, , drop = FALSE], 1, sum))
    possible_nbs[n_nbs_already < max_nbs]
  }

  adj_mat <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    # print(paste0('i = ', i))
    possible_nbs <- get_possible_nbs(i)
    n_nbs <- sum(adj_mat[i, ])
    if (n_nbs < min_nbs) {
      n_new_nbs <- sample(1:min(max_nbs - n_nbs, length(possible_nbs)))
      new_nbs <- sample(possible_nbs, n_new_nbs)
      adj_mat[i, new_nbs] <- 1
      adj_mat[new_nbs, i] <- 1
    }
  }
  adj_mat
}

standardise_precision_mat <- function(Sigma_inv) {
  Sigma <- solve(Sigma_inv)
  D <- diag(sqrt(diag(Sigma)))
  D %*% Sigma_inv %*% D
}

rcor_mat <- function(d, k0 = d, alphad = 1) {
  # K0: Sparsity level, number of correlated dimensions.

  if (k0 == 0) return(diag(rep(1, d)))

  Sigma <- effrcor::rcorrmatrix(k0, alphad = alphad)
  if (k0 != d) {
    identity_mat <- diag(rep(1, d - k0))
    zero_mat <- matrix(0, ncol = d - k0, nrow = k0)
    Sigma <- cbind(Sigma, zero_mat)
    Sigma <- rbind(Sigma, cbind(t(zero_mat), identity_mat))
  }
  structure(Sigma, 'which_dims_cor' = 1:k0)
}

car_mat_from_type <- function(p, precision_type, rho = 0.9, band = 2,
                          min_nbs = 1, max_nbs = 3) {
  if (precision_type == "lattice")
    nbs <- lattice_neighbours(p)
  else if (precision_type == "banded")
    nbs <- banded_neighbours(band, p)
  else if (precision_type == "random")
    nbs <- random_neighbours(p, min_nbs, max_nbs)
  car_precision_mat(nbs, rho = rho, standardised = FALSE, sparse = FALSE)
}

block_precision_mat <- function(p, m, within_block_type = "banded",
                                rho = 0.7, band = 2, sigma = 1,
                                min_nbs = 1, max_nbs = 3,
                                standardise = TRUE, sparse = TRUE) {
  block_Q <- function() {
    function(k) {
      car_mat_from_type(k, within_block_type, rho, band, min_nbs, max_nbs)
    }
  }

  if (length(m) == 1 && m < p) {
    remainder <- p %% m
    m <- rep(m, p %/% m)
    if (remainder > 0) m <- c(m, remainder)
  }

  Q <- diag(1, p)
  s <- cumsum(m)
  Q[1:s[1], 1:s[1]] <- block_Q()(m[1])
  if (length(m) > 1) {
    for (i in 2:length(m)) {
      if (m[i] > 1)
        Q[(s[i - 1] + 1):s[i], (s[i - 1] + 1):s[i]] <- block_Q()(m[i])
    }
  }
  if (standardise) Q <- standardise_precision_mat(Q)
  if (sparse) return(Matrix::Matrix(Q, sparse = TRUE))
  else return(Q)
}

Wishart_precision <- function(p, n = 5 * p, rho = 0.9, precision_type = "banded",
                              band = 2, min_nbs = 1, max_nbs = 3,
                              standardise = TRUE,  sparse = TRUE) {
  Sigma <- solve(car_mat_from_type(p, precision_type, rho, band, min_nbs, max_nbs))
  Q <- solve(1 / (n - 1) * rWishart(1, n, Sigma)[, , 1])
  if (standardise) Q <- standardise_precision_mat(Q)
  if (sparse) return(Matrix::Matrix(Q, sparse = TRUE))
  else return(Q)
}

#' Cuthill McKee (CM) algorithm
#'
#' Transform sparse matrix into a banded matrix. Code inspired by the netprioR package.
#'
#' @author Fabian Schmich
#' @param x Input matrix. Must be square and symmetric.
#' @return A list with the reordered matrix 'x' and the ordering 'order'.
cuthill_mckee <- function(x, reverse = FALSE, return_all = FALSE) {
  comp_nr <- 1:ncol(x)
  degree <- apply(x, 1, function(u) sum(u != 0))
  R <- rep(0, ncol(x))
  i <- 1
  R[i] <- comp_nr[which.min(degree)]
  l_R <- 1
  while (any(R == 0)) {
    current_neighbours <- which(x[R[i], ] != 0)
    Q <- setdiff(current_neighbours, R[1:l_R])
    if (length(Q) > 0) {
      R[(l_R + 1):(l_R + length(Q))] <- Q[order(degree[Q])]
      l_R <- l_R + length(Q)
    } else {
      R[l_R + 1] <- comp_nr[- R[1:l_R]][which.min(degree[- R[1:l_R]])]
      l_R <- l_R + 1
    }
    i <- i + 1
  }
  if (reverse) R <- rev(R)
  reordered_x <- x[R, R]
  if (return_all) return(list('x'       = reordered_x,
                              'order'   = R,
                              'reverse' = reverse,
                              'band'    = band(reordered_x)))
  else return(reordered_x)
}

compare_CM_alg <- function(p = 50, n_sim = 100) {
  n_computations <- function(B) {
    lower_nbs <- lapply(1:p, get_neighbours_less_than, sparse_mat = B)
    n_comp <- vapply(2:p, function(d) {
      2^length(remaining_neighbours_below(d, lower_nbs, p))
    }, numeric(1))
    sum(n_comp)
  }

  res <- matrix(0, ncol = 2, nrow = n_sim)
  for (i in 1:n_sim) {
    A <- radjacency_mat(p, max_nbs = 10)
    A_CM <- cuthill_mckee(A)
    A_reverse_CM <- cuthill_mckee(A, reverse = TRUE)
    res[i, 1] <- n_computations(A_CM)
    res[i, 2] <- n_computations(A_reverse_CM)
  }
  list('mean_computations' = res,
       'prop_CM_best'      = sum(apply(res, 1, function(x) x[1] < x[2])) / n_sim,
       'bands'             = c(band(A), band(A_CM), band(A_reverse_CM)))
}
