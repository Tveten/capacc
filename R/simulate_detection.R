#' @export
simulate_detection <- function(data = init_data(), method = method_params(),
                               seed = NULL, standardise_output = TRUE) {
  format_cpt_out <- function(cost, res) {
    if (res$value <= 0)
      return(data.table("cpt" = integer(0), "variate" = integer(0),
                        "mean1" = numeric(0), "mean2" = numeric(0)))
    else {
      if (cost == "sinspect") res$J <- which(res$proj != 0)
      mean1 <- colMeans(x$x[1:res$cpt, res$J, drop = FALSE])
      mean2 <- colMeans(x$x[(res$cpt + 1):nrow(x$x), res$J, drop = FALSE])
      return(data.table("cpt" = res$cpt, "variate" = res$J,
                        "mean1" = mean1, "mean2" = mean2))
    }
  }

  format_output <- function(cost, res) {
    if (cost == "cor_exact")
      return(list("collective" = collective_anomaliesR(res),
                  "point"      = point_anomaliesR(res)))
    else if (cost == "cor" || cost == "cor_sparse")
      return(list("collective" = collective_anomalies(list("anoms" = res)),
                  "point"      = point_anomalies(list("anoms" = res))))
    else if (cost == "iid" || cost == "decor")
      return(list("collective" = anomaly::collective_anomalies(res),
                  "point"      = anomaly::point_anomalies(res)))
    else if (cost %in% c("mvlrt", "sinspect"))
      return(format_cpt_out(cost, res))
    else if (cost == "inspect")
      return(anomalies_from_cpt(as.data.table(res)$location, x$x))
    else if (cost == "gflars")
      return(anomalies_from_cpt(res, x$x))
    else if (cost == "var_pgl")
      return(anomalies_from_cpt(res, x$x))
  }

  if (!is.null(seed)) set.seed(seed)
  x <- simulate_cor_(data)
  if (method$cost == "iid" || method$cost == "decor") {
    if (method$cost == "decor") x$x <- robust_whitening(x$x)
    else x$x <- robust_scale(x$x)
    # print(x$x[, 1])
    # print(round(robust_cov_mat(x$x), 2))
    beta <- iid_penalty(data$n, data$p, method$b)
    beta_tilde <- iid_point_penalty(data$n, data$p, max(0.05, method$b))
    res <- anomaly::capa.mv(x$x,
                            beta        = beta,
                            beta_tilde  = beta_tilde,
                            min_seg_len = method$minsl,
                            max_seg_len = method$maxsl,
                            type        = "mean")
  } else if (method$cost == "gflars") {
    res <- jointseg::jointSeg(x$x, method = "GFLars", K = round(data$n / 10))$bestBkp
  } else if (method$cost == "var_pgl") {
    res <- var_pgl(x$x, minsl = method$minsl, maxsl = method$maxsl)
  } else {
    Q_hat <- get_Q_hat(x$x, data, method)
    if (grepl("cor", method$cost)) {
      x$x <- centralise(x$x)
      if (method$cost == "cor_exact") capacc_func <- capacc_exact
      else if(method$cost == "cor")   capacc_func <- capacc
      else if(method$cost == "cor_sparse")   capacc_func <- capacc_sparse
      res <- capacc_func(x$x, Q_hat,
                         b                = method$b,
                         b_point          = max(0.05, method$b),
                         min_seg_len      = method$minsl,
                         max_seg_len      = method$maxsl)
    } else if (method$cost == "mvlrt")
      res <- cptcc(x$x, Q_hat,
                   b = method$b,
                   min_seg_len = method$minsl)
    else if (method$cost == "sinspect")
      res <- single_cor_inspect(t(x$x), Q_hat, method$b)
    else if (method$cost == "inspect")
      res <- cor_inspect(t(x$x), Q_hat, method$b)$changepoints
  }
  if (standardise_output) return(format_output(method$cost, res))
  else return(res)
}

#' @export
simulate_detection_known <- function(data = init_data(), method = method_params(),
                                     seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  x <- simulate_cor_(data)
  if (grepl("iid", method$cost) || method$cost == "decor") {
    if (method$cost == "decor") x$x <- robust_whitening(x$x)
    else x$x <- robust_scale(x$x)
    x_anom <- x$x[(data$locations + 1):(data$locations + data$durations), ]
    if (method$cost == "iid" || method$cost == "decor") {
      penalty_vec <- cumsum(iid_penalty(data$n, data$p, method$b))
      return(optimise_mvnormal_iid_saving(x_anom, penalty_vec))
    } else if (method$cost == "iid_dense") {
      alpha <- get_penalty("dense", data$n, data$p, method$b)$alpha_const
      return(list(S_max = saving_iid(1:data$p, x_anom) - alpha))
    }
  } else if (grepl("cor", method$cost) && method$cost != "decor") {
    Q_hat <- get_Q_hat(x$x, data, method)
    x$x <- centralise(x$x)
    x_anom <- x$x[(data$locations + 1):(data$locations + data$durations), ]
    if (method$cost == "cor" || method$cost == "cor_sparse") {
      if (method$cost == "cor") penalty <- get_penalty("combined", data$n, data$p, method$b)
      if (method$cost == "cor_sparse") penalty <- get_penalty("sparse", data$n, data$p, method$b)
      lower_nbs <- lower_nbs(Q_hat)
      extended_nbs <- extended_lower_nbs(lower_nbs)
      return(optimise_mvnormal_saving(x_anom, Q_hat, lower_nbs, extended_nbs,
                                      penalty$alpha_const, penalty$beta,
                                      penalty$alpha_lin))
    } else if (method$cost == "cor_exact") {
      penalty <- get_penalty_vec("combined", data$n, data$p, method$b)
      return(optimise_mvnormal_saving_BF(x_anom, Q_hat, penalty, mu_MLE(Q_hat)))
    } else if (method$cost == "cor_BF") {
      penalty <- get_penalty_vec("combined", data$n, data$p, method$b)
      return(optimise_mvnormal_saving_BF(x_anom, Q_hat, penalty, mu_aMLE()))
    } else if (method$cost == "cor_dense") {
      alpha <- get_penalty("dense", data$n, data$p, method$b)$alpha_const
      mean_x <- matrix(colMeans(x_anom), nrow = data$p)
      return(list(S_max = dense_mvnormal_savings(mean_x, Q_hat, nrow(x_anom)) - alpha))
    } else if (method$cost == "cor_dense_optimal") {
      alpha <- get_penalty("dense", data$n, data$p, method$b)$alpha_const
      mu <- generate_change(method$size_mu, data$p, 0)
      return(list(S_max = optimal_mvnormal_dense_saving(x_anom, Q_hat, mu, alpha)))
    }
  } else if (method$cost == "sinspect") {
    Q_hat <- get_Q_hat(x$x, data, method)
    single_cor_inspect(t(x$x), Q_hat, method$b, data$locations, standardize.series = TRUE)
  } else if (method$cost == "mvlrt") {
    Q_hat <- get_Q_hat(x$x, data, method)
    return(optimise_mvnormal_lr(data$locations, x$x, Q_hat, method$b))
  }
}

test_optim <- function(data, method) {
  x <- simulate_cor_(data)
  Q_hat <- get_Q_hat(x$x, data, method)
  x$x <- centralise(x$x)
  x_anom <- x$x[(data$locations + 1):(data$locations + data$durations), ]
  penalty <- get_penalty('combined', data$n, data$p, method$b)
  lower_nbs <- lower_nbs(Q_hat)
  extended_nbs <- extended_lower_nbs(lower_nbs)
  optimise_mvnormal_saving(x_anom, Q_hat, lower_nbs, extended_nbs,
                           penalty$alpha_const, penalty$beta,
                           penalty$alpha_lin)
}

#' @export
test_runtime <- function(cost, b = 1, est_band = 2, times = 5) {
  data <- init_data(p = 100, n = 200, vartheta = 0, rho = 0.9)
  method <- method_params(cost, b, precision_est_struct = "banded", est_band = est_band)
  print(microbenchmark::microbenchmark(
    simulate_detection(data, method, standardise_output = TRUE),
    times = times
  ))
}

