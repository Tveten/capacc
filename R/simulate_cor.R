
#' A function for generating simulated multivariate data
#'
#' @name simulate_cor
#'
#' @description Generates multivariate simulated data having n observations and p variates. The data have a standard Gaussian distribution except at
#' a specified number of locations where there is a change in mean in a proportion of the variates. The function is useful for generating data to demonstrate and assess
#' multivariate anomaly detection methods such as \code{capa.mv} and \code{pass}.
#'
#' @param n The number of observations. The default is \code{n=100}.
#' @param p The number of variates. The default is \code{p=10}.
#' @param mu The change in mean. Default is \code{mu=1}.
#' @param locations A vector of locations (or scalar for a single location) where the change in mean occurs. The default is \code{locations=20}.
#' @param durations A scalar or vector (the same length as \code{locations}) of values indicating the duration for the change in mean. If the durations are all
#' of the same length then a scalar value can be used. The default is \code{durations=20}.
#' @param proportions A scalar or vector (the same length as \code{locations}) of values in the range (0,1] indicating the proportion of variates at each location that are affected by
#' the change in mean. If the proportions are all same than a scalar value can be used. The default is \code{proportions=0.1}.
#'
#' @return A matrix with n rows and p columns
#'
#'
#' @examples
#' library(anomaly)
#' sim.data<-simulate(500,200,2,c(100,200,300),6,c(0.04,0.06,0.08))
#'
#' @export
simulate_cor <-function(n=100,p=10,vartheta=1,shape=0,change_seed=NA,
                        Sigma=diag(1, p),locations=40,durations=20,proportions=0.1,
                        change_type = 'adjacent', changing_vars = NA,
                        point_locations = NA, point_proportions = NA, point_mu = NA)
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
    if(!Reduce("&",locations+durations <= n))
    {
        stop("locations+durations must be <= n")
    }

    s <- locations
    e <- locations + durations
    X = MASS::mvrnorm(n, rep(0, p), Sigma)
    mu <- matrix(0, ncol = p, nrow = length(s))
    for (k in 1:length(s))
    {
      if (vartheta[k] > 0 && proportions[k] > 0) {
        J <- get_affected_dims(change_type, proportions[k], p, changing_vars)
        mu[k, ] <- generate_change(vartheta, J, shape, Sigma, change_seed)
        # mu[k, J] <- mu_J
        mu_mat <- t(replicate(e[k] - s[k], mu[k, ]))
        X[(s[k] + 1):e[k], ] <- X[(s[k] + 1):e[k], ] + mu_mat
      }
    }
    if (!any(is.na(point_locations))) {
      for (k in 1:length(point_locations)) {
        J <- get_affected_dims(change_type, point_proportions[k])
        X[point_locations[k], J] <- X[point_locations[k], J] + point_mu[k]
      }
    }
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
               point_mu          = data$point_mu)
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
    if (is.na(changing_vars))
      stop("If change_type is 'custom', you must provide a changing_vars vector")
    return(changing_vars)
  } else if (change_type == 'adjacent')
    return(1:k)
  else if (change_type == 'adjacent_lattice')
    return(unique(c(1, unlist(lattice_neighbours(p))))[1:k])
  else if (change_type == 'scattered')
    return(round(seq(2, p - 1, length.out = k)))
  else if (change_type == 'randomised')
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
  else if (shape == 6)
    mu_J <- MASS::mvrnorm(1, rep(0, length(J)), Sigma[J, J, drop = FALSE])
  mu <- rep(0, p)
  mu[J] <- mu_J
  mu / vector.norm(mu) * vartheta
}

