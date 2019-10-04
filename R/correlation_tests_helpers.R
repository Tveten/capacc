adjusted_penalty <- function(a, n, p, C = 2) {
  a <- as.numeric(a)
  psi <- a * log(n)
  # k_star <- sqrt(p) * psi / log(p)
  k_star <- (p + C * sqrt(p * psi)) / (2 * log(p))
  sum_beta <- C * psi + C * 1:p * log(p)
  sum_beta[1:p > k_star] <- p + C * psi + C * sqrt(p * psi)
  diff(c(0, sum_beta))
}

signal_strength <- function(mu, proportion, p, n, a) {
  psi <- a * log(n)
  k_star <- sqrt(p) * psi / log(p)
  size_J <- round(proportion * p)
  squared_norm_mu <- mu^2 * size_J

  if(size_J <= k_star) return(squared_norm_mu / (log(p) + psi / size_J))
  if(size_J > k_star) return(squared_norm_mu / (sqrt(p * psi) / size_J + psi / size_J))
}

min_duration <- function(a = 4, n = 10^4, p = 4, proportion = 0.5, mu = 1, C = 2) {
  delta_squared <- signal_strength(mu, proportion, p, n, a)
  40 * C / delta_squared
}
