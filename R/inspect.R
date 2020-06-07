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
  return(t((rightsums[t, , drop = FALSE] / (n-t) - leftsums[t, , drop = FALSE] / t) *
             sqrt(t * (n-t) / n)))
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
    vector.proj <- power.method(Mhat%*%t(Mhat), 1e-8)
  } else {
    tmp <- Mhat %*% power.method(t(Mhat)%*%Mhat, 1e-8)
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
single_cor_inspect <- function(x, Q, b = 1, cpt = NA, schatten=2, standardize.series=TRUE)
{
  x <- as.matrix(x)
  if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
  p <- dim(x)[1] # dimensionality of the time series
  n <- dim(x)[2] # time length of the observation
  if (standardize.series) {
    x <- rescale.variance(x)
    Q <- standardise_precision_mat(Q)
  }
  lambda <- sqrt(log(log(n) * p) / 2)
  threshold <- b * 4 * sqrt(log(n * p))  # Assumes x and Q are standardised.

  cusum.matrix <- cusum.transform(x)
  if (lambda >= max(abs(cusum.matrix))) lambda <- max(abs(cusum.matrix)) - 1e-10
  vector.proj <- as.numeric(Q %*% sparse.svd(cusum.matrix, lambda, schatten))
  vector.proj <- vector.proj / vector.norm(vector.proj)
  cusum.proj <- t(cusum.matrix)%*%vector.proj

  if (is.na(cpt)) {
    list("cpt"   = which.max(abs(cusum.proj)),
         "value" = max(abs(cusum.proj)) - threshold,
         "proj"  = vector.proj)
  } else {
    list("cpt"   = cpt,
         "S_max" = abs(cusum.proj[cpt]) - threshold,
         "proj"  = vector.proj)
  }
}

#' Computing threshold used in \code{inspect}
#' @description The threshold level to be used in \code{inspect} is computed via Monte Carlo simulation of multivariate time series that do not contain any changepoints.
#' @param n Time length of the observation.
#' @param p Dimension of the multivariate time series.
#' @param nrep Number of Monte Carlo repetition to be used.
#' @return A numeric value indicating the threshold level that should be used based on the Monte Carlo simulation.
#' @examples
#' compute.threshold(n=200, p=50)
compute.threshold <- function(n, Q, alpha = 0.05, nrep=100){
    p <- ncol(Q)
    cusum.stats = rep(0,nrep)
    for (i in 1:nrep) {
        x <- t(simulate_cor(n, p, Sigma = solve(Q), vartheta = 0)$x)
        cusum.stats[i] = locate.change(x, Q)$cusum
    }
    as.numeric(quantile(cusum.stats, probs = 1 - alpha))
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
locate.change <- function(x, Q, lambda, schatten=2, sample.splitting=FALSE,
                          standardize.series=FALSE, view.cusum=FALSE)
{
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation
    if (missing(lambda)) lambda <- sqrt(log(log(n)*p)/2)
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
    vector.proj <- vector.proj / vector.norm(vector.proj)
    cusum.proj <- t(cusum.matrix2)%*%vector.proj

    if (view.cusum) plot(as.numeric(cusum.proj), ylab='projected cusum', pch=20)

    ret <- NULL
    ret$changepoint <- which.max(abs(cusum.proj))
    if (sample.splitting) ret$changepoint <- ret$changepoint * 2
    ret$cusum <- max(abs(cusum.proj))
    ret$vector.proj <- vector.proj

    return(ret)
}

#' Informative sparse projection for estimation of changepoints (inspect)
#' @description This is the main function of the package InspectChangepoint. The function \code{inspect} estimates the locations of multiple changepoints in the mean structure of a multivariate time series. Multiple changepoints are estimated using a (wild) binary segmentation scheme, whereas each segmentation step uses the \code{\link{locate.change}} function.
#'
#' @usage inspect(x, lambda, threshold, schatten=c(1,2), M)
#'
#' @param x The input data matrix of a high-dimensional time series, with each component time series stored as a row.
#' @param lambda Regularisation parameter used in \code{\link{locate.change}}.  If no value is supplied, the dafault value is chosen to be log(log(n)*p/2), where p and n are the number of rows and columns of the data matrix x respectively.
#' @param threshold Threshold level for testing whether an identified changepoint is a true changepoint. If no value is supplied, the threshold level is computed via Monte Carlo simulation of 100 repetitions from the null model.
#' @param schatten The Schatten norm constraint to use in the \code{\link{locate.change}} function. Default is schatten = 2, i.e. a Frobenius norm constraint.
#' @param M The Monte Carlo parameter used for wild binary segmentation. Default is M = 0, which means a classical binary segmentation scheme is used.
#'
#' @details The input time series is first standardised using the \code{\link{rescale.variance}} function. Recursive calls of the \code{\link{locate.change}} function then segments the multivariate time series using (wild) binary segmentation. A changepoint at time z is defined here to mean that the time series has constant mean structure for time up to and including z and constant mean structure for time from z+1 onwards.
#'
#' More details about model assumption and theoretical guarantees can be found in Wang and Samworth (2016). Note that Monte Carlo computation of the threshold value can be slow, especially for large p. If \code{inspect} is to be used multiple times with the same (or similar) data matrix size, it is better to precompute the threshold level via Monte Carlo simulation by calling the \code{\link{compute.threshold}} function.
#'
#' @return The return value is an S3 object of class 'inspect'. It contains a list of two objeccts:
#' \itemize{
#' \item{x }{The input data matrix}
#' \item{changepoints }{A matrix with three columns. The first column contains the locations of estimated changepoints sorted in increasing order; the second column contains the maximum CUSUM statistics of the projected univariate time series associated with each estimated changepoint; the third column contains the depth of binary segmentation for each detected changepoint.}
#' }
#'
#' @references Wang, T. and Samworth, R. J. (2018) High dimensional changepoint estimation via sparse projection. \emph{J. Roy. Statist. Soc., Ser. B}, \strong{80}, 57--83.
#'
#' @examples
#' n <- 500; p <- 100; ks <- 30; zs <- c(125,250,375)
#' varthetas <- c(0.1,0.15,0.2); overlap <- 0.5
#' obj <- multi.change(n, p, ks, zs, varthetas, overlap)
#' x <- obj$x
#' threshold <- compute.threshold(n,p)
#' ret <- inspect(x, threshold = threshold)
#' ret
#' summary(ret)
#' plot(ret)
#' @import stats
#' @import graphics
#' @export
cor_inspect <- function(x, Q, b, lambda, schatten=c(1, 2), M){
    # basic parameters and initialise
    x <- as.matrix(x)
    if (dim(x)[2] == 1) x <- t(x) # treat univariate time series as a row vector
    p <- dim(x)[1] # dimensionality of the time series
    n <- dim(x)[2] # time length of the observation
    if (missing(lambda)) lambda <- sqrt(log(log(n)*p)/2)
    if (missing(b) || is.na(b) || is.null(b)) threshold <- compute.threshold(n, Q)
    else threshold <- b * 4 * sqrt(log(n * p))  # Assumes x and Q are standardised.
    if (missing(schatten)) schatten <- 2
    if (missing(M)) M <- 100
    x <- rescale.variance(x)
    Q <- standardise_precision_mat(Q)

    # generate random time windows of length at least 2
    rnd1 <- sample(0:(n-2), M, replace = TRUE)
    rnd2 <- sample(0:(n-2), M, replace = TRUE)
    window_s <- pmin(rnd1, rnd2)
    window_e <- pmax(rnd1, rnd2) + 2

    # recursive function for binary segmentation
    BinSeg <- function(x, s, e, depth, parent.val){
        if (e - s <= 2) return(NULL) # stop when the segment has only one point
        ind <- (window_s >= s) & (window_e <= e) # \mathcal{M}_{s,e}
        max.val <- -1
        cp <- 0

        for (m in c(0,((1:M)[ind]))) {
            if (m == 0) {
                s_m <- s
                e_m <- e
            } else {
                s_m <- window_s[m]
                e_m <- window_e[m]
            }
            obj <- locate.change(x[,(s_m+1):e_m, drop = FALSE], Q)
            if (obj$cusum > max.val) {
                max.val <- obj$cusum
                cp <- s_m + obj$changepoint
            }
        }

        # recurse
        ret <- NULL
        ret$location <- cp
        ret$max.proj.cusum <- max.val #min(parent.val, max.val)
        ret$depth <- depth
        if (ret$max.proj.cusum < threshold) {
            return(NULL)
        } else {
            return(cbind(BinSeg(x, s, cp, depth + 1, ret$max.proj.cusum),
                         ret,
                         BinSeg(x, cp, e, depth + 1, ret$max.proj.cusum)))
        }
    }

    # return all changepoints of x
    ret <- NULL
    ret$x <- x
    ret$changepoints <- BinSeg(x, 0, n, depth=1, parent.val=.Machine$double.xmax)
    ret$changepoints <- t(matrix(as.numeric(ret$changepoints), nrow = 3))
    colnames(ret$changepoints) = c('location', 'max.proj.cusum', 'depth')
    class(ret) <- 'inspect'
    return(ret)
}

anomalies_inspect <- function(res, x, tol = 1) {
  res <- as.data.table(res)
  n <- nrow(x)
  res <- rbind(data.table(location = 0, max.proj.cusum = 0, depth = 0),
               res,
               data.table(location = n, max.proj.cusum = 0, depth = 0))
  res$mean_size <- 0
  for (i in 2:length(res$location)) {
    seg_mean <- colMeans(x[(res$location[i - 1] + 1):res$location[i], , drop = FALSE])
    res$mean_size[i] <- sign(sum(seg_mean)) * sqrt(sum(seg_mean^2))
  }
  starts <- integer(0)
  ends <- integer(0)
  in_anom <- FALSE
  curr_start_ind <- 0
  i <- 2
  while (i <= nrow(res)) {
    if (!in_anom && abs(res$mean_size[i]) >= tol) {
      curr_start_ind <- i - 1
      starts <- c(starts, res$location[curr_start_ind] + 1)
      in_anom <- TRUE
      i <- i + 1
    } else {
      if (in_anom) {
        end_anom <- is_in_interval(res$mean_size[i], c(-tol, tol))
        switch_anom <- (res$mean_size[curr_start_ind + 1] < 0 && res$mean_size[i] > tol) ||
          (res$mean_size[curr_start_ind + 1] > 0 && res$mean_size[i] < - tol)
        if (end_anom) {
          ends <- c(ends, res$location[i - 1])
          in_anom <- FALSE
        } else if (switch_anom) {
          ends <- c(ends, res$location[i - 1])
          curr_start_ind <- i - 1
          starts <- c(starts, res$location[curr_start_ind] + 1)
        }
      }
      i <- i + 1
    }
  }
  if (in_anom) ends <- c(ends, res[.N, location])
  if (length(starts) != length(ends)) {
    print(starts)
    print(ends)
    stop("Bug when extracting inspect anomalies. Unequal number of start and end points.")
  }
  anoms <- data.table(start = starts, end = ends)
  return(list("collective" = anoms[start != end],
              "point"      = data.table(location = anoms[start == end, start])))
}

