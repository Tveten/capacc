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

car_precision_mat <- function(nbs = list('random', 10), rho = 0.5, sigma = 1,
                              min_nbs = 1, max_nbs = 3, standardised = TRUE) {
  if (nbs[[1]][1] == 'random') {
    p <- nbs[[2]]
    W <- radjacency_mat(p, min_nbs, max_nbs)
  } else {
    p <- length(nbs)
    W <- adjacency_mat(nbs)
  }
  eigen_values_W <- eigen(W)$values
  if (rho <= 1 / eigen_values_W[p]) {
    min_rho <- 1 / eigen_values_W[p]
    warning(paste0('rho has been altered to 1 / smallest_eigen_value(W) = ', min_rho))
    rho <- min_rho + sqrt(.Machine$double.eps)
  }
  Sigma_inv <- 1 / sigma^2 * (diag(as.vector(W %*% rep(1, p))) - rho * W)
  if (standardised) return(standardise_precision_mat(Sigma_inv))
  else return(Sigma_inv)
}

car_precision_mat2 <- function(nbs, m = 1, rho = 0.5, sigma = 1) {
  p <- length(nbs)
  W <- adjacency_mat(nbs)
  rowsum_W <- rowSums(W)
  W_pluss <- W / rowsum_W
  M_pluss_inv <- diag(rowsum_W)
  eigen_values_W <- eigen(W)$values
  if (rho <= 1 / eigen_values_W[p])
    stop(paste0('rho must be greater than 1 / smallest_eigen_value(W) = ', 1 / eigen_values_W[p]))
  1 / sigma^2 * M_pluss_inv %*% (diag(1, p) - rho * W_pluss)
}

adjacency_mat <- function(nbs_list) {
  p <- length(nbs_list)
  adj_mat <- do.call('rbind', lapply(nbs_list, indicator, p = p))
  if (!isSymmetric(adj_mat)) {
    adj_mat <- symmetric_from_lower_tri(adj_mat)
  }
  adj_mat
}

symmetric_from_lower_tri <- function(x_tri) {
  upper_triangle <- upper.tri(x_tri)
  x_tri[upper_triangle] <- t(x_tri)[upper_triangle]
  x_tri
}

banded_neighbours <- function(band, p) {
  nbs <- list()
  for (i in 1:p) {
    if (i == 1) lower <- integer(0)
    else        lower <- max(1, i - band):(i - 1)
    if (i == p) upper <- integer(0)
    else        upper <- (i + 1):min(i + band, p)
    nbs[[i]] <- c(lower, upper)
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

random_neighbours_banded <- function(p, min_nbs = 1, max_nbs = 3) {
  max_nbs <- min(max_nbs, p - 1)
  min_nbs <- min(max_nbs, max(1, min_nbs))

  n_nbs <- vapply(1:(p - 1), function(i) sample(min_nbs:max_nbs, 1), numeric(1))
  number_neighbours(n_nbs)
}

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

  if (return_all) return(list('x' = x[R, R], 'order' = R, 'reverse' = reverse))
  else return(x[R, R])
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
    A <- radjacency_mat(p, max_nbs = 5)
    A_CM <- cuthill_mckee(A)
    A_reverse_CM <- cuthill_mckee(A, reverse = TRUE)
    res[i, 1] <- n_computations(A_CM)
    res[i, 2] <- n_computations(A_reverse_CM)
  }
  list('mean_computations' = res,
       'prop_CM_best'       = sum(apply(res, 1, function(x) x[1] < x[2])) / n_sim)
}
