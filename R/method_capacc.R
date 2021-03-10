#' Collective and point anomalies in cross-correlated data---CAPA-CC
#'
#' @name capa.cc
#'
#' @description
#' A method for detecting anomalous segments and points based on CAPA-CC by Tveten, Eckley, Fearnhead (2020).
#'
#' @param x An n x p data matrix where each row is an observation vector.
#' @param Q An estimate of the precision matrix. See \code{\link{robust_sparse_precision}}. Must be a sparse matrix from the Matrix package.
#' @param b The scaling factor for the collective anomaly penalty. Defaults to 1.
#' @param b_point The scaling factor for the point anomaly penalty. Defaults to 1.
#'
#' @param min_seg_len The minimum segment length. Defaults to 2.
#' @param max_seg_len The maximum segment length. Defaults to 10^8.
#' @param transform A function used to centre the data prior to analysis by \code{capa.cc}. This can, for example, be used to compensate for the effects of autocorrelation in the data.
#' The default value is \code{transform=centralise}, which centrailises the data by the median.
#' Other choices are available in the anomaly package, for example \code{anomaly::robustscale} and \code{anomaly::ac_corrected}), and a user defined function can also be specified.
#'
#' @return An S3 class of type capacc with the following components:
#' \describe{
#' \item{\code{x}}{The input data matrix.}
#' \item{\code{anoms}}{A data frame with four columns: start (start-point of the anomaly), end (end-point of an anomaly), variate (which variable is affected) and size (the estimated size of the mean component for the given variate).}
#' }
#'
#' @references
#'
#' @examples
#' library(capacc)
#' x <- simulate_cor()$x
#' Q <- robust_sparse_precision(x, adjacency_mat(banded_neighbours(2, ncol(x)), sparse = FALSE))
#' res <- capa.cc(x, Q, b = 1, min_seg_len = 5)
#' plot(res)
#' collective_anomalies(res)
#'
#' @export
#'
capa.cc<-function(x,Q,b=1,b_point=1,min_seg_len=2,max_seg_len=10^8,transform=centralise)
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
    # try transforming the data
    if(!purrr::is_function(transform))
    {
        stop("transform must be a function")
    }
    x<-transform(x)
    # and check the transformed data
    if(!is.array(x))
    {
        stop("cannot convert the transformed x to an array - check the transform function")
    }
    if(any(vapply(x, is.na, logical(1))))
    {
        stop("the transformed x contains NA values - check the transform function")
    }
    if(any(vapply(x, is.null, logical(1))))
    {
        stop("the transformed x contains NULL values - check the transform function")
    }
    if(!is.numeric(x))
    {
        stop("the transformed x must be of type numeric - check the transform function")
    }
    if(length(dim(x)) != 2)
    {
        stop("cannot convert transformed x to a two dimensional array - check the transform function")
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

    # check max_seg_len
    if(is.null(max_seg_len))
    {
        max_seg_len=dim(x)[1]
    }
    if(!is.numeric(max_seg_len))
    {
        stop("max_seg_len must be numeric")
    }
    if(max_seg_len < 1)
    {
        stop("max_seg_len must be greater than zero")
    }
    # check relative values of min and max segment length
    if(max_seg_len < min_seg_len)
    {
        stop("max_seg_len must be greater than min_seg_len")
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
    # check b_point
    if(!is.numeric(b_point))
    {
        stop("b_point must be numeric")
    }
    if(length(b_point) != 1)
    {
        stop("b_point must be a single scalar value")
    }
    if(b_point < 0)
    {
        stop("b_point must be >= 0")
    }

    anoms <- capacc(x, Q, b, b_point, min_seg_len, max_seg_len)
    invisible(structure(list("x" = x, "anoms" = anoms), class = 'capacc'))
}

#### C++ helpers #####

#' @export
collective_anomalies <- function(capacc_res) {
  res_dt <- as.data.table(capacc_res$anoms)
  res_dt <- res_dt[start != end]
  names(res_dt)[4] <- "mean_change"
  res_dt
}

#' @export
point_anomalies <- function(capacc_res) {
  res_dt <- as.data.table(capacc_res$anoms)
  res_dt <- res_dt[start == end][, 2:4]
  names(res_dt)[c(1, 3)] <- c("location", "strength")
  res_dt
}



######################
#### R functions #####
######################

init_precision_mat <- function(precision_mat, tol = sqrt(.Machine$double.eps)) {
  p <- ncol(precision_mat)
  precision_mat[abs(precision_mat) <= tol] <- 0
  lower_nbs <- lower_nbs(precision_mat)
  extended_nbs <- extended_lower_nbs(lower_nbs)
  list('Q' = precision_mat, 'nbs' = lower_nbs, 'extended_nbs' = extended_nbs)
}

capaccR <- function(x, Q, b = 1,
                        min_seg_len = 3, max_seg_len = round(nrow(x) / 3),
                        prune = TRUE, print_progress = FALSE) {
  # TODO: Restrict min_seg_len and max_seg_len
  # Robustly normalise x
  n <- nrow(x)
  p <- ncol(x)
  # x <- anomaly::robustscale(x)
  # TODO: Implement standardisation based on Q.
  # Now: Assumes x and Q are standardised before input.
  # precision_mat <- robust_sparse_precision(x, adj_model)
  if (!any(class(Q) == "dsCMatrix")) Q <- Matrix::Matrix(Q, sparse = TRUE)
  Q_obj <- init_precision_mat(Q)
  # Q_order <- 1:p

  x <- rbind(rep(0, p), x)  # So the index m always means m - 1.

  # Setting up penalties
  penalty_names <- c("sparse", "dense", "point")
  penalty <- lapply(penalty_names, get_penalty, n = n, p = p, b = b)
  names(penalty) <- penalty_names

  # Initialising the DP for the optimal cost.
  C <- rep(0, n + 1)
  S <- matrix(0, nrow = n + 1, ncol = n + 1)
  J_max <- list()
  S_point <- rep(0, n + 1)
  J_point <- list()
  anom <- matrix(0, nrow = n + 1, ncol = 2)
  prune_t <- rep(0, n + 1)

  # Running OP
  for (m in (min_seg_len + 1):(n + 1)) {
    if (print_progress)
      if (m %% 10 == 0) print(paste0(round(100 * m/(n + 1), 2), '% complete.'))
    ## Collective anomalies:
    ts <- max(m - min_seg_len - max_seg_len + 1, 1):(m - min_seg_len)
    ts <- setdiff(ts, unique(prune_t))
    J_m <- list()
    for (t in ts) {
      saving <- optimise_mvnormal_saving(x[(t + 1):m, , drop = FALSE], Q,
                                         Q_obj$nbs,
                                         Q_obj$extended_nbs,
                                         penalty$dense$alpha,
                                         penalty$sparse$beta,
                                         penalty$sparse$alpha)

      # J_m[[t]] <- Q_order[saving$J_max]  # To refer to the original variates.
      J_m[[t]] <- saving$J_max
      S[m, t] <- saving$S_max
    }
    C1s <- C[ts] + S[m, ts]
    C1_max <- max(C1s)
    t_max <- ts[which.max(C1s)]
    J_max[[m]] <- J_m[[t_max]]

    ## Point anomalies:
    point_saving <- optimise_mvnormal_saving(x[m, , drop = FALSE], Q,
                                             Q_obj$nbs,
                                             Q_obj$extended_nbs,
                                             Inf,
                                             penalty$point$beta,
                                             penalty$point$alpha)
    # J_point[[m]] <- Q_order[point_saving$J_max]  # To refer to the original variates.
    J_point[[m]] <- point_saving$J_max
    S_point[m] <- point_saving$S_max
    C2 <- C[m - 1] + S_point[m]

    C0 <- C[m - 1]
    C[m] <- max(C0, C1_max, C2)

    if (prune) {
      t_to_prune <- ts[C0 >= C1s + penalty$dense$alpha]
      # print("Prune:")
      # print(t_to_prune)
      prune_t[t_to_prune] <- t_to_prune
    }

    # anomaly_type == 0: Previous cost is maximal. Points one time-step back.
    # anomaly_type == 1: Current cost for a collective anomaly is maximal. Points to the maximising time-step.
    # anomaly_type == 2: Current cost for a point anomaly is maximal. Points to the current time-step.
    anomaly_type <- which.max(c(C0, C1_max, C2)) - 1
    if (anomaly_type == 0) {
      anom[m, ] <- c(m - 1, anomaly_type)
    } else if (anomaly_type == 1) {
      anom[m, ] <- c(t_max, anomaly_type)
    } else if (anomaly_type == 2) {
      anom[m, ] <- c(m - 1, anomaly_type)
    }
  }
  # return(list('x'       = x[-1, order(Q_order)],
  return(list('x'       = x[-1, ],
              'S'       = S,
              'J'       = J_max,
              'S_point' = S_point,
              'J_point' = J_point,
              'C'       = C,
              'Q'       = Q,
              'anom'    = anom))
}

capacc_exact <- function(x, Q, b = 1, b_point = 1,
                             min_seg_len = 2, max_seg_len = 100,
                             prune = TRUE, print_progress = FALSE) {
  # TODO: Restrict min_seg_len and max_seg_len
  # Robustly normalise x
  n <- nrow(x)
  p <- ncol(x)
  # x <- anomaly::robustscale(x)
  # TODO: Implement standardisation based on Q.
  # Now: Assumes x and Q are standardised before input.
  if (!any(class(Q) == "dsCMatrix")) Q <- Matrix::Matrix(Q, sparse = TRUE)
  x <- rbind(rep(0, p), x)  # So the index m always means m - 1.
  mu_est <- mu_MLE(Q)

  # # Setting up penalties
  penalty <- get_penalty_vec("combined", n, p, b)
  point_penalty <- get_penalty_vec("point", n, p, b_point)

  # Initialising the DP for the optimal cost.
  C <- rep(0, n + 1)
  S <- matrix(0, nrow = n + 1, ncol = n + 1)
  J_max <- list()
  S_point <- rep(0, n + 1)
  J_point <- list()
  anom <- matrix(0, nrow = n + 1, ncol = 2)
  prune_t <- rep(0, n + 1)

  # Running OP
  for (m in (min_seg_len + 1):(n + 1)) {
    if (print_progress)
      if (m %% 10 == 0) print(paste0(round(100 * m/(n + 1), 2), '% complete.'))
    ## Collective anomalies:
    ts <- max(m - min_seg_len - max_seg_len + 1, 1):(m - min_seg_len)
    ts <- setdiff(ts, unique(prune_t))
    J_m <- list()
    for (t in ts) {
      saving <- optimise_mvnormal_saving_BF(x[(t + 1):m, ], Q, penalty, mu_est)
      J_m[[t]] <- saving$J_max
      S[m, t] <- saving$S_max
    }
    C1s <- C[ts] + S[m, ts]
    C1_max <- max(C1s)
    t_max <- ts[which.max(C1s)]
    J_max[[m]] <- J_m[[t_max]]

    ## Point anomalies:
    point_saving <- optimise_mvnormal_saving_BF(x[m, , drop = FALSE], Q,
                                                point_penalty, mu_est)
    J_point[[m]] <- point_saving$J_max
    S_point[m] <- point_saving$S_max
    C2 <- C[m - 1] + S_point[m]

    C0 <- C[m - 1]
    C[m] <- max(C0, C1_max, C2)

    if (prune) {
      t_to_prune <- ts[C0 >= C1s + penalty$dense$alpha]
      prune_t[t_to_prune] <- t_to_prune
    }

    # anomaly_type == 0: Previous cost is maximal. Points one time-step back.
    # anomaly_type == 1: Current cost for a collective anomaly is maximal. Points to the maximising time-step.
    # anomaly_type == 2: Current cost for a point anomaly is maximal. Points to the current time-step.
    anomaly_type <- which.max(c(C0, C1_max, C2)) - 1
    if (anomaly_type == 0) {
      anom[m, ] <- c(m - 1, anomaly_type)
    } else if (anomaly_type == 1) {
      anom[m, ] <- c(t_max, anomaly_type)
    } else if (anomaly_type == 2) {
      anom[m, ] <- c(m - 1, anomaly_type)
    }
  }
  return(list('x'       = x[-1, ],
              'S'       = S,
              'J'       = J_max,
              'S_point' = S_point,
              'J_point' = J_point,
              'C'       = C,
              'Q'       = Q,
              'anom'    = anom))
}


collective_anomaliesR <- function(capacc_res) {
  anom_list <- list()
  m <- nrow(capacc_res$anom)
  # Indices m correspond to m - 1 in the data.
  while (m >= 1) {
    if (capacc_res$anom[m, 2] == 0) m <- capacc_res$anom[m, 1]
    else if (capacc_res$anom[m, 2] == 1) {
      end <- m - 1
      start <- capacc_res$anom[m, 1]
      J <- capacc_res$J[[m]]
      means <- colMeans(capacc_res$x[start:end, J, drop = FALSE])
      anom_df <- data.frame('start'       = rep(start, length(J)),
                            'end'         = rep(end, length(J)),
                            'variate'     = J,
                            'mean_change' = means)
      anom_list[[length(anom_list) + 1]] <- anom_df
      m <- start
    } else if (capacc_res$anom[m, 2] == 2) {
      m <- capacc_res$anom[m, 1]
    }
  }
  if (length(anom_list) == 0)
    return(data.frame('start' = integer(0), 'end' = integer(0),
                      'variate' = integer(0), 'mean_change' = integer(0)))
  else {
    anom_dt <- data.table::as.data.table(do.call('rbind', anom_list))
    anom_dt <- anom_dt[order(start)][, .SD[order(variate)], by = start]
    return(as.data.frame(anom_dt))
  }
}

point_anomaliesR <- function(capacc_res) {
  anom_list <- list()
  m <- nrow(capacc_res$anom)
  while (m >= 1) {
    # if (capacc_res$anom[m, 2] == 0) m <- capacc_res$anom[m, 1]
    # else if (capacc_res$anom[m, 2] == 1) {
    #   m <- capacc_res$anom[m, 1]
    # } else if (capacc_res$anom[m, 2] == 2) {
    if (capacc_res$anom[m, 2] == 2) {
      J_point <- capacc_res$J_point[[m]]
      anom_df <- data.frame('location' = rep(m - 1, length(J_point)),
                            'variate'  = J_point,
                            'strength' = capacc_res$x[m - 1, J_point])
      anom_list[[length(anom_list) + 1]] <- anom_df
      # m <- capacc_res$anom[m, 1]
    }
    m <- capacc_res$anom[m, 1]
  }
  if (length(anom_list) == 0)
    return(data.frame('location' = integer(0), 'variate' = integer(0),
                      'mean_change' = numeric(0)))
  else {
    anom_dt <- data.table::as.data.table(do.call('rbind', anom_list))
    anom_dt <- anom_dt[order(location)][, .SD[order(variate)], by = location]
    return(as.data.frame(anom_dt))
  }
}
