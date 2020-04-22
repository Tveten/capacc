
saving_iid <- function(J, x) {
  means <- colMeans(x[, J])
  nrow(x) * sum(means^2)
}

find_threshold <- function(data = init_data(), alpha = 0.99, n_sim = 10^4) {
  current_thresholds <- fread("./results/alpha_dense.csv")
  a <- alpha
  current_thresholds <- current_thresholds[n == data$n & p == data$p & alpha == a & precision_type == data$precision_type & band == data$band & rho == data$rho]
  if (nrow(current_thresholds) > 0) stop("A threshold for this data setup already exists.")

  res_dt <- data.table("value" = rep(0, 2 * n_sim),
                       "cost" = c(rep("cor", n_sim),
                                  rep("iid", n_sim)))
  for (i in 1:n_sim) {
    x <- simulate_cor(n = data$n, p = data$p, mu = 0, Sigma = data$Sigma)
    res_dt$value[i] <- dense_mvnormal_savings(matrix(colMeans(x[1:2, ]), nrow = data$p),
                                              data$Sigma_inv, 2)
    res_dt$value[n_sim + i] <- saving_iid(1:data$p, x[1:2, ])
  }
  threshold_dt <- res_dt[, .(threshold      = quantile(value, alpha),
                             alpha          = alpha,
                             n              = data$n,
                             p              = data$p,
                             precision_type = data$precision_type,
                             band           = data$band,
                             rho            = data$rho),
                         by = cost]
  fwrite(threshold_dt, "./results/alpha_dense.csv", append = TRUE)
}


#' @export
saving_distr <- function(s, e, data = init_data(n = 100, p = 10, mu = 0, rho = 0.9),
                         b = 1, n_sim = 10^4, a = 0.99, seed = NULL) {
  get_title <- function(s, e) {
    if (is_in_interval(data$locations + 1, c(s, e)) | is_in_interval(data$locations + data$durations, c(s, e))) {
      mu <- data$mu
      prop <- data$proportions
    } else {
      mu <- 0
      prop <- 0
    }
    paste0("2-banded, rho=", data$rho, ", p=", data$p, ", n=", data$n,
           ", s=", s, ", e=", e, ", mu=", data$mu, ", prop=", data$proportions)
  }

  get_threshold <- function(data, a) {
    alpha <- fread("./results/alpha_dense.csv")
    alpha <- alpha[n == data$n & p == data$p & alpha == a & precision_type == data$precision_type & band == data$band & rho == data$rho]
    if (nrow(alpha) == 0) {
      message("Finding alpha_dense. May take time.")
      find_threshold(data, alpha = a)
      return(get_threshold(data, a))
    } else return(alpha)
  }

  if (!is.null(seed)) set.seed(seed)
  res_dt <- data.table("value" = rep(0, 3 * n_sim),
                       "cost" = c(rep("cor", n_sim),
                                  rep("iid", n_sim),
                                  rep("chisq_p", n_sim)))
  for (i in 1:n_sim) {
    x <- simulate_cor_(data)
    res_dt$value[i] <- dense_mvnormal_savings(matrix(colMeans(x[s:e, ]), nrow = data$p),
                                              data$Sigma_inv, e - s + 1)
    res_dt$value[n_sim + i] <- saving_iid(1:data$p, x[s:e, ])
  }
  res_dt$value[(2 * n_sim + 1):(3 * n_sim)] <- rchisq(n_sim, data$p)
  alpha <- get_threshold(data, a)
  print(paste0("cor detection rate: ", res_dt[cost == "cor", sum(value > alpha[cost == "cor", threshold]) / n_sim]))
  print(paste0("iid detection rate: ", res_dt[cost == "iid", sum(value > alpha[cost == "iid", threshold]) / n_sim]))

  ggplot2::ggplot(data = res_dt, ggplot2::aes(value, colour = cost)) +
    ggplot2::geom_density() +
    ggplot2::scale_x_continuous("Saving") +
    ggplot2::ggtitle(get_title(s, e)) +
    ggplot2::geom_vline(xintercept = alpha[, threshold],
                        color = c("green", "blue"), size=1.5)
}
