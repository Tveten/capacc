
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
    if (is_in_interval(data$locations, c(s, e)) | is_in_interval(data$locations + data$durations, c(s, e))) {
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
    x <- simulate_cor(n                 = data$n,
                      p                 = data$p,
                      mu                = data$mu,
                      Sigma             = data$Sigma,
                      locations         = data$locations,
                      durations         = data$durations,
                      proportions       = data$proportions,
                      change_type       = data$change_type,
                      point_locations   = data$point_locations,
                      point_proportions = data$point_proportions,
                      point_mu          = data$point_mu)
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

# find_penalty_scale <- function(setup = init_data(), alpha = 0.01, n_sim = 10^4) {
#   current_scales <- fread("./results/penalty_scales.csv")
#   current_scales <- current_scales[n == setup$n & p == setup$p & alpha == a & cor_mat_type == setup$cor_mat_type & band == setup$band & rho == setup$rho]
#   if (nrow(current_scales) > 0) stop("A penalty scale for this setup already exists.")
#
#   res_dt <- data.table("value" = rep(0, 2 * n_sim),
#                        "type" = c(rep("cor", n_sim), rep("iid", n_sim)))
#   for (i in 1:n_sim) {
#     x <- simulate_cor(n = 2, p = setup$p, mu = 0, Sigma = setup$Sigma, locations = 1, durations = 1)
#     res_dt$value[i] <- dense_mvnormal_savings(matrix(colMeans(x), nrow = setup$p),
#                                               setup$Sigma_inv, 2)
#     res_dt$value[n_sim + i] <- saving_iid(1:setup$p, x)
#   }
#   alpha_dense <- get_penalty("dense", setup$n, setup$p)$alpha
#   q <- 1 - alpha / setup$n
#   scale_dt <- res_dt[, .(b            = quantile(value, q) / alpha_dense,
#                          alpha        = alpha,
#                          n            = setup$n,
#                          p            = setup$p,
#                          cor_mat_type = setup$cor_mat_type,
#                          band         = setup$band,
#                          rho          = setup$rho),
#                      by = type]
#   fwrite(scale_dt, "./results/penalty_scales.csv", append = TRUE)
# }
#
# find_penalty_scale_sparse <- function(setup = init_data(), alpha = 0.01, n_sim = 10^4) {
#   current_scales <- fread("./results/penalty_scales_sparse.csv")
#   current_scales <- current_scales[n == setup$n & p == setup$p & alpha == a & cor_mat_type == setup$cor_mat_type & band == setup$band & rho == setup$rho]
#   if (nrow(current_scales) > 0) stop("A penalty scale for this setup already exists.")
#
#   res_dt <- data.table("value" = rep(0, 2 * n_sim),
#                        "type" = c(rep("cor", n_sim), rep("iid", n_sim)))
#   Q <- setup$Sigma_inv
#   for (i in 1:n_sim) {
#     x <- simulate_cor(n = 2, p = setup$p, mu = 0, Sigma = setup$Sigma, locations = 1, durations = 1)
#     mean_x <- colMeans(x)
#     res_dt$value[i] <- max(2 * ((2 * mean_x * Q %*% mean_x) - diag(as.matrix(Q)) * mean_x^2))
#     res_dt$value[n_sim + i] <- max(2 * mean_x^2)
#   }
#   penalty <- get_penalty("sparse", setup$n, setup$p)
#   q <- 1 - alpha / setup$n
#   scale_dt <- res_dt[, .(b            = quantile(value, q) / (penalty$alpha + penalty$beta),
#                          alpha        = alpha,
#                          n            = setup$n,
#                          p            = setup$p,
#                          cor_mat_type = setup$cor_mat_type,
#                          band         = setup$band,
#                          rho          = setup$rho),
#                      by = type]
#   fwrite(scale_dt, "./results/penalty_scales_sparse.csv", append = TRUE)
# }
#
#
# check_fpr <- function(setup = init_data(), scale_by = "sparse", b = NULL,
#                       n_sim = 10^2, a = 0.01, seed = NULL) {
#   get_b <- function(setup, a) {
#     if (scale_by == "dense") file_name <- "./results/penalty_scales.csv"
#     if (scale_by == "sparse") file_name <- "./results/penalty_scales_sparse.csv"
#     bs <- fread(file_name)
#     bs <- bs[n == setup$n & p == setup$p & alpha == a & cor_mat_type == setup$cor_mat_type & band == setup$band & rho == setup$rho]
#     if (nrow(bs) == 0) {
#       message("Finding penalty scale. May take time.")
#       if (scale_by == "dense") find_penalty_scale(setup, alpha = a)
#       if (scale_by == "sparse") find_penalty_scale_sparse(setup, alpha = a)
#       return(get_b(setup, a))
#     } else return(bs)
#   }
#
#   costs <- c("cor", "iid")
#   if (is.null(b)) bs <- get_b(setup, a)
#   else bs <- data.table("type" = costs, "b" = b)
#   print(bs)
#   seeds <- sample(1:10^7, n_sim)
#   res <- do.call("rbind", lapply(seeds, function(seed) {
#     do.call("rbind", lapply(costs, function(cost) {
#       sim_res <- simulate_mvcapa(setup,
#                                  cost_type        = cost,
#                                  b                = bs[type == cost, b],
#                                  seed             = seed,
#                                  return_anom_only = TRUE)
#       data.table("type" = cost, "fp" = !is.na(sim_res$collective$start[1]))
#     }))
#   }))
#   res[, .(fp = mean(fp)), by = type]
# }
