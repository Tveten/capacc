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
