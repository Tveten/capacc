linear_penalty <- function(n, p, b = 1, a = 2, C = 2) {
  a <- as.numeric(a)
  psi <- a * log(n)
  k_star <- (p + C * sqrt(p * psi)) / (2 * log(p))
  # list('beta' = b * (C * psi / k_star + C * log(p)), 'k_star' = k_star)
  list('alpha' = b * C * psi, 'beta' = b * C * log(p), 'k_star' = k_star)
}

const_penalty <- function(n, p, b = 1, a = 2, C = 2) {
  a <- as.numeric(a)
  psi <- a * log(n)
  b * (p + C * psi + C * sqrt(p * psi))
}

get_penalty <- function(penalty_regime, n, p, b = 1, a = 2) {
  if (penalty_regime == 'sparse') {
    penalty_obj <- linear_penalty(n, p, b = b, a = a)
    beta <- penalty_obj$beta
    alpha <- penalty_obj$alpha
    k_star <- penalty_obj$k_star
  } else if (penalty_regime == 'dense') {
    beta <- 0
    alpha <- const_penalty(n, p, b = b, a = a)
    k_star <- p + 1
  } else if (penalty_regime == 'point') {
    penalty_obj <- linear_penalty(n, p, b = b, a = a)
    beta <- penalty_obj$beta + penalty_obj$alpha
    alpha <- 0
    k_star <- p + 1
  }
  list('alpha' = alpha, 'beta' = beta, 'k_star' = k_star)
}

combined_penalty_vec <- function(n, p, b = 1, a = 2) {
  sparse_penalty <- get_penalty('sparse', n, p, b, a)
  dense_penalty <- get_penalty('dense', n, p, b, a)
  penalty_vec <- 0:p * sparse_penalty$beta + sparse_penalty$alpha
  penalty_vec[penalty_vec > dense_penalty$alpha] <- dense_penalty$alpha
  penalty_vec
}

iid_penalty <- function(n, p, b = 1) {
    s = 2 * log(n)
    a_vector = qchisq(seq(p-1,0)/p, 1)

    penalty_1 = 2*s + 2*1:p*log(p)
    penalty_2 = rep(p + 2*s + 2*sqrt(p*s),p)

    penalty_3    = 2*(s+log(p)) + 1:p + 2*p*a_vector*dchisq(a_vector, 1) + 2*sqrt((1:p+2*p*a_vector*dchisq(a_vector, 1))*(s+log(p)))
    penalty_3[p] = 2*(s+log(p)) + p + 2*sqrt(p*(s+log(p)))

    beta = b * diff( c(0,pmin(penalty_1,penalty_2,penalty_3)))
    beta
}

iid_point_penalty <- function(n, p, b = 1) {
  point_pen <- get_penalty("point", n, p, b)
  point_pen$beta
}

