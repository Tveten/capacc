#' Noise standardisation for multivariate time series.
#' @description Each row of the input matrix is normalised by the estimated standard deviation computed through the median absolute deviation of increments.
#' @param x An input matrix of real values.
#' @details This is an auxiliary function used by the \code{InspectChangepoint} package.
#' @return A rescaled matrix of the same size is returned.
#' @examples
#' x <- matrix(rnorm(40),5,8) * (1:5)
#' x.rescaled <- rescale.variance(x)
#' x.rescaled
#' @export
rescale.variance <- function(x, by_row = FALSE){
    p <- dim(x)[1]
    n <- dim(x)[2]
    for (j in 1:p){
        scale <- mad(diff(x[j,]))/sqrt(2)
        x[j,] <- x[j,] / scale
    }
    return(x)
}


#' CUSUM transformation
#' @description Performing CUSUM transformation to the input matrix of multivariate time series. If the input is a vector, it is treated as a matrix of one row.
#' @details For any integers p and n, the CUSUM transformation \eqn{T_{p,n}: R^{p\times n}\to R^{p\times (n-1)}} is defined by
#' \deqn{
#'    [T_{p,n}(M)]_{j,t} := \sqrt{t(n-t)/n}\biggl(\frac{1}{n-t}\sum_{r=t+1}^n M_{j,r} - \frac{1}{t}\sum_{r=1}^t M_{j,r}\biggr).
#' }
#'
#' @param x input matrix
#' @return The transformed matrix is returned. Note that the returned matrix has the same number of rows but one fewer columns compared with the input matrix.
#' @examples
#' x <- matrix(rnorm(20),4,5)
#' cusum.transform(x)
#' @export
cusum.transform <- function(x){
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation

    leftsums <- t(apply(x, 1, cumsum))
    rightsums <- leftsums[, n] - leftsums
    leftsums <- t(leftsums)
    rightsums <- t(rightsums)
    t <- 1:(n - 1)

    # constructing CUSUM matrix
    return(t((rightsums[t,] / (n-t) - leftsums[t,] / t) * sqrt(t * (n-t) / n)))
}


#' Projection onto the standard simplex
#' @description The input vector is projected onto the standard simplex, i.e. the set of vectors of the same length as the input vector with non-negative entries that sum to 1.
#' @param v Input vector
#' @details This is an auxiliary function used by the \code{InspectChangepoint} package.
#' @return A vector in the standard simplex that is closest to the input vector is returned.
#' @references Chen, Y. and Ye, X. (2011) Projection onto a simplex. arXiv preprint, arxiv:1101.6081.
#' @examples
#' v <- rnorm(10)
#' PiW(v)
#' @export
PiW <- function(v){
    # projection of a vector to the standard simplex (non-negative entries
    # summing to 1)
    ord = order(v, decreasing = TRUE);
    v = v[ord]
    s = 0
    for (i in 1:length(v)){
        s = s + v[i]
        if (i == length(v)) {delta = (s-1)/i; break}
        d = s - v[i+1]*i
        if (d >= 1) {delta = (s - 1)/i; break}
    }
    v[1:i] = v[1:i] - delta; v[-(1:i)] = 0;
    w = rep(0, length(v))
    w[ord] = v
    w
}


#' Matrix projection onto the nuclear norm unit sphere
#' @description Projection (with respect to the inner product defined by the Frobenius norm) of a matrix onto the unit sphere defined by the nuclear norm.
#' @details This is an auxiliary function used by the \code{InspectChangepoint} package. The projection is achieved by first performing a singular value decomposition, then projecting the vector of singular values onto the standard simplex, and finally using singular value decomposition in reverse to build the projected matrix.
#' @param M Input matrix
#' @return A matrix of the same dimension as the input is returned.
#' @examples
#' M <- matrix(rnorm(20),4,5)
#' PiS(M)
#' @export
PiS <- function(M){
    # projection of a matrix (in Frobenius norm) to the ball of nuclear norm 1
    tmp = svd(M);
    tmp$u%*%diag(PiW(tmp$d))%*%t(tmp$v)
}


#' Soft thresholding a vector
#' @param x a vector of real numbers
#' @param lambda soft thresholding value
#' @return a vector of the same length
#' @description entries of v are moved towards 0 by the amount lambda until they hit 0.
#' @export
vector.soft.thresh <- function(x, lambda){
  sign(x)*pmax(0,(abs(x) - lambda))
}


#' Norm of a vector
#' @description Calculate the entrywise L_q norm of a vector or a matrix
#' @param v a vector of real numbers
#' @param q a nonnegative real number or Inf
#' @param na.rm boolean, whether to remove NA before calculation
#' @return the entrywise L_q norm of a vector or a matrix
#' @export
vector.norm <- function(v, q=2, na.rm=FALSE){
    if (na.rm) v <- na.omit(v)
    if (q > 0) (sum(abs(v)^q))^(1/q)
    else if (q == 0) sum(v!=0)
    else if (q == Inf) max(abs(v))
    else NaN
}


#' Power method for finding the leading eigenvector of a symmetric matrix
#' @param A a square symmetric matrix
#' @param eps tolerance for convergence (in Frobenius norm)
#' @param maxiter maximum iteration
#' @return a unit-length leading eigenvector of A
#' @export
power.method <- function(A, eps = 1e-10, maxiter = 10000){
  if (nrow(A) != ncol(A)) stop('powerMethod requires a square matrix')
  if (!all.equal(A, t(A))) stop('powerMethod requires a symmetric matrix')

  d <- nrow(A)
  v <- rnorm(d); v <- v/vector.norm(v)

  for (i in seq_len(maxiter)){
    v_new <- A%*%v; v_new <- v_new/vector.norm(v_new)
    if ((vector.norm(v_new - v)) < eps) break
    v <- v_new
  }
  if (i == maxiter) warning('max iter reached without convergence.')
  return(as.numeric(v_new))
}


#' Computing the sparse leading left singular vector of a matrix
#'
#' @description Estimating the sparse left leading singular vector by first computing a maximiser Mhat of the convex problem
#' \deqn{<Z, M> - \lambda |M|_1}
#' subject to the Schatten norm constraint |M|_schatten <= 1 using alternating direction method of multipliers (ADMM). Then the leading left singular vector of Mhat is returned.
#'
#' @details In case of schatten = 2, a closed-form solution for Mhat using matrix soft thresholding is possible. We use the closed-form solution instead of the ADMM algorithm to speed up the computation.
#'
#' @param Z Input matrix whose left leading singular vector is to be estimated.
#' @param lambda Regularisation parameter
#' @param schatten Schatten norm constraint to be used. Default uses Schatten-2-norm, i.e. the Frobenius norm. Also possible to use Schatten-1-norm, the nuclear norm.
#' @param tolerance Tolerance criterion for convergence of the ADMM algorithm. Not used when shatten=2.
#' @param max.iter Maximum number of iteration in the ADMM algorithm. Not used when shatten=2.
#'
#' @return A vector that has the same length as nrow(Z) is returned.
#' @examples
#' Z <- matrix(rnorm(20),4,5)
#' lambda <- 0.5
#' sparse.svd(Z, lambda)
#' @export
sparse.svd <- function(Z, lambda, schatten=c(1, 2), tolerance=1e-5, max.iter=10000){
    if (missing(schatten)) schatten <- 2
    if (schatten == 2){
        # with Frobenius norm constraint, the sparse vector is obtained by soft
        # thresholding
        Mhat <- vector.soft.thresh(Z, lambda)
    } else {
        # with nuclear norm constraint, the sparse vector is obtained by ADMM
        p <- dim(Z)[1]; n <- dim(Z)[2]; gamma <- 1;
        X <- matrix(0,p,n); Y <- matrix(0,p,n); U <- matrix(0,p,n)
        iter <- 0
        while ((iter < max.iter) | (max.iter == 0)) {
            iter <- iter + 1
            X <- PiS(Y - U + gamma * Z)
            Y <- vector.soft.thresh(X + U, lambda * gamma)
            U <- U + (X - Y)
            if (vector.norm(X - Y) < tolerance) break
        }
        Mhat <- X
    }

    if (nrow(Mhat) < ncol(Mhat)){
        vector.proj <- power.method(Mhat%*%t(Mhat), 1e-5)
    } else {
        tmp <- Mhat %*% power.method(t(Mhat)%*%Mhat, 1e-5)
        vector.proj <- tmp/vector.norm(tmp)
    }

    return(vector.proj)
}


#' Single changepoint estimation
#' @description Estimate the location of one changepoint in a multivariate time
#' series. It uses the function \code{\link{sparse.svd}} to estimate the best
#' projection direction, then using univariate CUSUM statistics of the projected
#' time series to estimate the changepoint location.
#' @param x A (p x n) data matrix of multivariate time series, each column
#' represents a data point
#' @param lambda Regularisation parameter. If no value is supplied, the dafault
#' value is chosen to be sqrt(log(log(n)*p/2)) for p and n number of rows and
#' columns of the data matrix x respectively.
#' @param schatten The Schatten norm constraint to use in the \code{\link{sparse.svd}}
#'  function. Default is schatten = 2, i.e. a Frobenius norm constraint.
#' @param sample.splitting Whether the changepoint should be estimated via
#' sample splitting. The theoretical result is proven only for the sample
#' splitted version of the algorithm. However, the default setting in practice
#' is without sample splitting.
#' @param standardize.series Whether the given time series should be
#' standardised before estimating the projection direction. Default is FALSE,
#' i.e. the input series is assume to have variance 1 in each coordinate.
#' @param view.cusum Whether to show a plot of the projected CUSUM series
#' @return A list of two items:
#' \itemize{
#'   \item changepoint - A single integer value estimate of the changepoint
#'   location is returned. If the estimated changepoint is z, it means that the
#'   multivariate time series is piecewise constant up to z and from z+1
#'   onwards.
#'   \item cusum - The maximum absolute CUSUM statistic of the projected
#'   univariate time series associated with the estimated changepoint.
#'   \item vector.proj - the vector of projection, which is proportional to an estimate of the vector of change.
#' }
#' @references Wang, T., Samworth, R. J. (2016) High-dimensional changepoint estimation via sparse projection. Arxiv preprint: arxiv1606.06246.
#' @examples
#' n <- 2000; p <- 1000; k <- 32; z <- 400; vartheta <- 0.12; sigma <- 1; shape <- 3
#' noise <- 0; corr <- 0
#' obj <- single.change(n,p,k,z,vartheta,sigma,shape,noise,corr)
#' x <- obj$x
#' locate.change(x)
#' @export

single_cor_inspect <- function(x, Q, b = 1, schatten=2, sample.splitting=FALSE,
                               standardize.series=FALSE, view.cusum=FALSE)
{
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation
    lambda <- b * sqrt(log(log(n)*p)/2)
    threshold <- lambda
    if (standardize.series) x <- rescale.variance(x)
    if (sample.splitting){
        x1 <- x[,seq(1,n,by=2)]
        x2 <- x[,seq(2,n,by=2)]
    } else {
        x1 <- x
        x2 <- x
    }

    # construct cusum matrix of x
    cusum.matrix1 <- cusum.transform(x1)
    if (sample.splitting) {
        cusum.matrix2 <- cusum.transform(x2)
    } else {
        cusum.matrix2 <- cusum.matrix1
    }

    # estimate changepoint
    if (lambda >= max(abs(cusum.matrix1))) lambda <- max(abs(cusum.matrix1)) - 1e-10

    vector.proj <- as.numeric(Q %*% sparse.svd(cusum.matrix1, lambda, schatten))
    cusum.proj <- t(cusum.matrix2)%*%vector.proj

    if (view.cusum) plot(as.numeric(cusum.proj), ylab='projected cusum', pch=20)

    ret <- NULL
    ret$cpt <- which.max(abs(cusum.proj))
    if (sample.splitting) ret$changepoint <- ret$changepoint * 2
    ret$value <- max(abs(cusum.proj))
    ret$vector.proj <- vector.proj

    return(ret)
}
