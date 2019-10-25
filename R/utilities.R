setup_parallel <- function() {
  n_cores <- parallel::detectCores()
  c <- parallel::makeCluster(n_cores - 1, outfile = '', type = 'PSOCK')
  doParallel::registerDoParallel(c)
  c
}

stop_parallel <- function(c) {
  parallel::stopCluster(c)
}

extract_nested_element <- function(i, list_of_lists) {
  lapply(list_of_lists, function(list_in_list) list_in_list[[i]])
}

off_diag_ind <- function(nr, p) {
  A <- matrix(FALSE, nrow = p, ncol = p)
  delta <- row(A) - col(A)
  for (n in nr) {
    A[delta == n] <- TRUE
  }
  A
}

off_diag <- function(nr, x) {
  x[off_diag_ind(nr, nrow(x))]
}

#' Cuthill McKee (CM) algorithm
#'
#' Transform sparse matrix into a banded matrix. Code inspired by the netprioR package.
#'
#' @author Fabian Schmich
#' @param x Input matrix. Must be square and symmetric.
#' @return A list with the reordered matrix 'x' and the ordering 'order'.
cuthill_mckee <- function(x, reverse = FALSE) {
  comp_nr <- 1:ncol(x)
  degree <- apply(x, 1, function(u) sum(u != 0))
  R <- rep(0, ncol(x))
  R[1] <- comp_nr[which.min(degree)]
  i <- 1
  while (any(R == 0)) {
    i <- i + 1
    current_neighbours <- which(x[R[(i - 1)], ] != 0)
    Q <- setdiff(current_neighbours, R[1:(i - 1)])
    if (length(Q) > 0) R[i:(i + length(Q) - 1)] <- Q[order(degree[Q])]
    else R[i] <- comp_nr[- R[1:(i - 1)]][which.min(degree[- R[1:(i - 1)]])]
  }
  if (reverse) R <- rev(R)

  list('x' = x[R, R], 'order' = R, 'reverse' = reverse)
}

count_nz_columns <- function(x) {
  m <- nrow(x)
  sum(apply(x, 2, function(u) isTRUE(all.equal(u, rep(0, m)))))
}

is_zero <- function(u) {
  vapply(u, function(j) isTRUE(all.equal(j, 0)), logical(1))
}

indicator <- function(J, p) {
  x <- rep(0, p)
  x[J] <- 1
  x
}
