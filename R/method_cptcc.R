#' Changepoints in cross-correlated data---CAPA-CC
#'
#' @name cpt.cc
#'
#' @description
#' A method for detecting a single changepoint in cross-correlated data based on CPT-CC by Tveten, Eckley, Fearnhead (2020).
#' To detect multiple changepoints, combine this function with a version of Binary Segmentation (BS), like Wild BS or Seeded BS.
#'
#' @param x An n x p data matrix where each row is an observation vector.
#' @param Q An estimate of the precision matrix. See \code{\link{robust_sparse_precision}}. Must be a sparse matrix from the Matrix package.
#' @param b The scaling factor for the collective anomaly penalty. Defaults to 1.
#' @param min_seg_len The minimum segment length. Defaults to 2.
#'
#' @return An S3 class of type cptcc with the following components:
#' \describe{
#' \item{\code{x}}{The input data matrix.}
#' \item{\code{cpts}}{A data frame with two columns: variate (which variable is affected) and cpt (the location of the estimated changepoint).}
#' }
#'
#' @references
#'
#' @examples
#' library(capacc)
#' x <- simulate_cor(n = 100, locations = 50, durations = 50, proportions = 0.2)$x
#' diff_x <- x[2:nrow(x), ] - x[1:(nrow(x) - 1), ]
#' Q <- 2 * robust_sparse_precision(diff_x, adjacency_mat(banded_neighbours(2, ncol(x)), sparse = FALSE))
#' res <- cpt.cc(x, Q, b = 1, min_seg_len = 5)
#'
#' @export
#'

cpt.cc<-function(x,Q,b=1,min_seg_len=2)
{
    # make sure the callable transform object is of type function
    # data needs to be in the form of an array
    x<-to_array(x)
    if(!is.array(x))
    {
        stop("cannot convert x to an array")
    }
    if(any(vapply(x, is.na, logical(1))))
    {
        stop("x contains NA values")
    }
    if(any(vapply(x, is.null, logical(1))))
    {
        stop("x contains NULL values")
    }
    if(!is.numeric(x))
    {
        stop("x must be of type numeric")
    }
    # check dimensions
    if(length(dim(x)) != 2)
    {
        stop("cannot convert x to a two dimensional array")
    }
    if(dim(x)[1] == 1)
    {
        x<-t(x)
    }

    # check Q
    if (ncol(x) != ncol(Q))
    {
      stop("The number of columns in x are not equal to the number of columns in Q")
    }
    if (ncol(Q) != nrow(Q))
    {
      stop("Q must be a square matrix")
    }
    if (!is(Q, 'sparseMatrix'))
    {
      stop("Q is not a sparse matrix from the Matrix package.")
    }
    if (!Matrix::isSymmetric(Q))
    {
      stop("Q must be a symmetrix matrix.")
    }

    # check min_seg_len
    if(!is.numeric(min_seg_len))
    {
        stop("min_seg_len must be numeric")
    }
    if(min_seg_len < 1)
    {
        stop("min_seg_len must be greater than zero")
    }
    if(min_seg_len > dim(x)[1])
    {
        stop("min_seg_len must be less tha the number of observations in x")
    }

    # check b
    if(!is.numeric(b))
    {
        stop("b must be numeric")
    }
    if(length(b) != 1)
    {
        stop("b must be a single scalar value")
    }
    if(b < 0)
    {
        stop("b values must be >= 0")
    }

    cpts <- cptcc(x, Q, b, min_seg_len)
    cpts_res <- data.frame("variate" = cpts$J, "cpt" = cpts$cpt)
    invisible(structure(list("x" = x, "cpts" = cpts_res), class = 'cptcc'))
}
