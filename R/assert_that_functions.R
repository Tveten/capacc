is_positive_definite <- function(cov_mat, tol = .Machine$double.eps) {
  eigen_values <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  all(eigen_values >= tol)
}

assert_cov_mat <- function(cov_mat) {
  cov_mat_name <- deparse(substitute(cov_mat))
  na_msg <- paste0(cov_mat_name, ' cannot contain NAs.')
  assertthat::assert_that(!any(is.na(cov_mat)), msg = na_msg)
  numeric_msg <- paste0(cov_mat_name, ' must be numeric.')
  assertthat::assert_that(is.numeric(cov_mat), msg = numeric_msg)
  matrix_msg <- paste0(cov_mat_name, ' must be of class "matrix".')
  assertthat::assert_that(class(cov_mat) == 'matrix', msg = matrix_msg)
  symmetric_msg <- paste0(cov_mat_name, ' is not a symmetric matrix.')
  assertthat::assert_that(isSymmetric(cov_mat), msg = symmetric_msg)
  posdef_msg <- paste0(cov_mat_name, ' is not a positive definite matrix (some eigenvalues are < 1e-16).')
  assertthat::assert_that(is_positive_definite(cov_mat), msg = posdef_msg)
}

assert_equal_ncol <- function(x, y) {
  x_name <- deparse(substitute(x))
  y_name <- deparse(substitute(y))
  ncol_msg <- paste0(x_name, ' and ', y_name, ' must have an equal number of columns.')
  assertthat::assert_that(ncol(x) == ncol(y), msg = ncol_msg)
}

is_whole_number <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}

is_in_interval <- function(x, interval) {
  x >= interval[1] & x <= interval[2]
}

assert_integer_in_interval <- function(x, interval) {
  x_name <- deparse(substitute(x))
  int_interval_msg <- paste0(x_name, ' must be an integer between ', interval[1],
                         ' and ', interval[2], '.')
  assertthat::assert_that(all(is_in_interval(x, interval)),
                          all(is_whole_number(x)),
                          msg = int_interval_msg)
}
