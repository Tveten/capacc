
#' A function for generating simulated multivariate data
#'
#' @name simulate_cor
#'
#' @description Generates cross-correlated multivariate simulated data having n observations and p variates. The data have a Gaussian distribution with the specified covariance matrix except at
#' a specified number of locations where there is a change in mean in a proportion of the variates. The function is useful for generating data to demonstrate and assess
#' multivariate anomaly detection methods such as \code{capa.cc}, \code{capa.mv} \code{pass} and \code{inspect}.
#'
#' @param n The number of observations. The default is \code{n=100}.
#' @param p The number of variates. The default is \code{p=10}.
#' @param vartheta The size of the mean change vector in L2 distance. Defaults to 5.
#' @param shape An integer between 0 and 10 specifying the shape of a change. Defaults to 0, which means equally changing components. 5 gives mean components drawn from an i.i.d. Gaussian distriubtion, while 6 draws changes from the data distribution. See the function \code{generate_change} for more details.
#' @param change_seed The seed of the drawn mean change.
#' @param Sigma The data covariance matrix. The default is the identity matrix.
#' @param locations A vector of locations (or scalar for a single location) where the change in mean occurs. The default is \code{locations=40}.
#' @param durations A scalar or vector (the same length as \code{locations}) of values indicating the duration for the change in mean. If the durations are all
#' of the same length then a scalar value can be used. The default is \code{durations=20}.
#' @param proportions A scalar or vector (the same length as \code{locations}) of values in the range (0,1] indicating the proportion of variates at each location that are affected by
#' the change in mean. If the proportions are all same than a scalar value can be used. The default is \code{proportions=0.1}.
#' @param change_type A string specifying which variables are affected. Options include "adjacent", "adjacent_lattice", "scattered", "block_scattered", "custom" and "random". See the function \code{get_affected_dims} for more details.
#' @param changing_vars If \code{change_type="custom"}, which variables are anomalous?
#' @param point_locations A vector with locations of point anomalies. Defaults to NA.
#' @param point_proportions A vector of the same length as \code{point_locations} specifying the the proportion of variables affected by each point anomaly.
#' @param point_mu A vector of the same length as \code{point_locations} specifying the mean of all variables affected by each point anomaly.
#'
#' @return A matrix with n rows and p columns
#'
#'
#' @examples
#' library(anomaly)
#' sim.data<-simulate(500,200,2,c(100,200,300),6,c(0.04,0.06,0.08))
#'
#' @export
simulate_cor <-function(n=100,p=10,vartheta=5,shape=0,change_seed=NA,
                        Sigma=diag(1, p),locations=40,durations=20,proportions=0.1,
                        change_type = 'adjacent', changing_vars = NA,
                        point_locations = NA, point_proportions = NA, point_mu = NA,
                        n_sd_changes = 0)
{
    if(length(n) > 1)
    {
        stop("n must be a scalar")
    }
    if(length(p) > 1)
    {
        stop("p must be a scalar")
    }
    if(n < 1)
    {
        stop("n must be > 0")
    }
    if(p < 1)
    {
        stop("p must be > 0")
    }
    if(!Reduce("&&",locations > 0))
    {
        stop("location values must all be > 0")
    }
    if(!Reduce("&&",durations > 0))
    {
        stop("durations must all be > 0")
    }
    if(!Reduce("&&",proportions <= 1))
    {
        stop("proportion values  must all be <= 1")
    }
    if(length(proportions) != 1 && length(proportions) != length(locations))
    {
        stop("proportions must be a scalar or a vector the same size as locations")
    }
    if(length(durations) != 1 && length(durations) != length(locations))
    {
        stop("durations must be a scalar or a vector the same size as locations")
    }
    if(length(vartheta) != 1 && length(vartheta) != length(locations))
    {
        stop("vartheta must be a scalar or a vector the same size as locations")
    }
    if(length(change_type) != 1 && length(change_type) != length(locations))
    {
        stop("change_type must be a scalar or a vector the same size as locations")
    }
    if(length(durations) == 1)
    {
        durations <- rep(durations,length(locations))
    }
    if(length(proportions) == 1)
    {
        proportions <- rep(proportions,length(locations))
    }
    if(length(vartheta) == 1)
    {
        vartheta <- rep(vartheta, length(locations))
    }
    if(length(change_type) == 1)
    {
        change_type <- rep(change_type, length(locations))
    }
    if(!Reduce("&",locations+durations <= n))
    {
        stop("locations+durations must be <= n")
    }
    if (length(point_proportions) == 1)
    {
      point_proportions <- rep(point_proportions, length(point_locations))
    }
    if (length(point_mu) == 1)
    {
      point_mu <- rep(point_mu, length(point_locations))
    }

    s <- locations
    e <- locations + durations
    X = MASS::mvrnorm(n, rep(0, p), Sigma)
    mu <- matrix(0, ncol = p, nrow = length(s))
    for (k in 1:length(s))
    {
      if (vartheta[k] > 0 && proportions[k] > 0) {
        J <- get_affected_dims(change_type[k], proportions[k], p, changing_vars)
        mu[k, ] <- generate_change(vartheta[k], J, shape, Sigma, change_seed)
        mu_mat <- t(replicate(e[k] - s[k], mu[k, ]))
        X[(s[k] + 1):e[k], ] <- X[(s[k] + 1):e[k], ] + mu_mat
      }
    }
    if (!any(is.na(point_locations))) {
      for (k in 1:length(point_locations)) {
        J <- get_affected_dims("random", point_proportions[k], p, changing_vars)
        X[point_locations[k], J] <- X[point_locations[k], J] + point_mu[k]
      }
    }
    if (n_sd_changes > 0) X <- add_changing_variance(X, n_sd_changes)
    list("x" = X, "mu" = mu)
}

simulate_cor_ <- function(data = init_data()) {
  simulate_cor(n                 = data$n,
               p                 = data$p,
               vartheta          = data$vartheta,
               shape             = data$shape,
               change_seed       = data$change_seed,
               Sigma             = data$Sigma,
               locations         = data$locations,
               durations         = data$durations,
               proportions       = data$proportions,
               change_type       = data$change_type,
               changing_vars     = data$changing_vars,
               point_locations   = data$point_locations,
               point_proportions = data$point_proportions,
               point_mu          = data$point_mu,
               n_sd_changes      = data$n_sd_changes)
}

# - Start from location[j] + 1 to locations[j] + s[j] to be consistent with
# notation in paper.
# - Locations, durations and proportions that scale automatically with p and n.
# - Description of elements in capa.mv.class (object returned from capa.mv) or
#   how to access them in documentation of capa.mv.
# - Source of many false positives:
#     * Contaminated baseline estimate.
#     * Error in penalty? psi must have proportionality constant > 2 or > 3?
# - Possibility to simulate data with no anomalies.
# - High-dim plot which also shows J?

get_affected_dims <- function(change_type, prop, p, changing_vars) {
  k <- ceiling(prop * p)
  if (k < 1) return(integer(0))
  if (change_type == "custom") {
    if (any(is.na(changing_vars)))
      stop("If change_type is 'custom', you must provide a changing_vars vector")
    return(changing_vars)
  } else if (change_type == 'adjacent')
    return(1:k)
  else if (change_type == 'adjacent_lattice')
    return(unique(c(1, unlist(lattice_neighbours(p))))[1:k])
  else if (change_type == 'scattered')
    return(round(seq(1, p, length.out = k)))
  else if (change_type == "block_scattered") {
    n_blocks <- min(k, 3)
    block_sizes <- rep(k %/% n_blocks, n_blocks)
    n_residual <- k %% n_blocks
    if (n_residual > 0)
      block_sizes[1:n_residual] <- block_sizes[1:n_residual] + 1
    block_starts <- round(seq(1, p - block_sizes[n_blocks] + 1,
                          length.out = n_blocks))
    J <- do.call("c", lapply(1:n_blocks, function(i) {
      block_starts[i]:(block_starts[i] + block_sizes[i] - 1)
    }))
    return(J)
  } else if (change_type == 'random')
    return(sample(1:p, k))
}

generate_change <- function(vartheta, J, shape, Sigma, seed = NA) {
  if (!is.na(seed)) set.seed(seed)
  p <- ncol(Sigma)
  k <- length(J)
  if      (shape == 0) mu_J <- rep(1,k)
  else if (shape == 1) mu_J <- k:1
  else if (shape == 2) mu_J <- (k:1)^2
  else if (shape == 3) mu_J <- (1:k)^(-1/2)
  else if (shape == 4) mu_J <- rchisq(k, 2)
  else if (shape == 5) mu_J <- rnorm(k)
  else {
    if      (shape == 6) Sigma_c <- Sigma
    else if (shape == 7) Sigma_c <- constant_cor_mat(ncol(Sigma), 0.7)
    else if (shape == 8) Sigma_c <- constant_cor_mat(ncol(Sigma), 0.8)
    else if (shape == 9) Sigma_c <- constant_cor_mat(ncol(Sigma), 0.9)
    else if (shape == 10) Sigma_c <- constant_cor_mat(ncol(Sigma), 0.99)
    mu_J <- MASS::mvrnorm(1, rep(0, length(J)), Sigma_c[J, J, drop = FALSE])
  }
  mu <- rep(0, p)
  mu[J] <- mu_J
  mu / vector.norm(mu) * vartheta
}

add_changing_variance <- function(x, n_changes = 10, sd_range = c(0.5, 2)) {
  p <- ncol(x)
  n <- nrow(x)
  cpt <- round(seq(0, n, length.out = n_changes + 1))
  sd <- runif(n_changes, sd_range[1], sd_range[2])
  for (i in 2:length(cpt)) {
    x[(cpt[i - 1] + 1):cpt[i], ] <- sd[i - 1] * x[(cpt[i - 1] + 1):cpt[i], ]
  }
  x
}






