
#' @export
tuning_params <- function(alpha = 0.05, tol = 0.02, max_iter = 50,
                          init_b = c(0.05, 1, 10), n_sim = 200) {
  list("alpha"           = alpha,
       "alpha_tol"       = tol,
       "tuning_max_iter" = max_iter,
       "init_b"          = init_b,
       "tuning_n_sim"    = n_sim)
}

tune_penalty <- function(data = init_data(mu = 0), params = mvcapa_params(),
                         tuning = tuning_params(), seed = NA) {

  add_setup_info <- function(res) {
    res <- cbind(res,
                 as.data.table(data[!grepl("Sigma", names(data))]),
                 as.data.table(params[names(params) != "b"]),
                 as.data.table(tuning[names(tuning) != "init_b"]))
    res$seed <-  seed
    res
  }

  fp_rate <- function(b) {
    fps <- unlist(lapply(1:tuning$tuning_n_sim, function(i) {
      params$b <- b
      !is.na(simulate_mvcapa(data, params, return_anom_only = TRUE)$collective$start[1])
    }))
    data.table("b" = b, "fp" = mean(fps), "diff" = mean(fps) - tuning$alpha)
  }

  new_b <- function(res) {
    res <- res[order(diff)]
    b_lower <- as.numeric(res[diff <= 0, min(b)])
    b_upper <- as.numeric(res[diff > 0][diff == min(diff), max(b)])
    if (is.na(b_lower) || is.na(b_upper))
      warning("Will not converge. init_b does not cover the selected alpha level.")
    exp((log(b_lower) + log(b_upper)) / 2)
  }

  if (!is.na(seed)) set.seed(seed)
  data$mu <- 0
  message(paste0("Tuning ", params$cost,  " penalty for n = ", data$n, ", p = ", data$p,
                 ", cor_mat_type = ", data$cor_mat_type,
                 ", band = ", data$band,
                 " and rho = ", data$rho))
  res <- do.call('rbind', Map(fp_rate, tuning$init_b))
  print(res)
  while (!any(abs(res$diff) <= tuning$alpha_tol + 1e-6) && nrow(res) <= tuning$tuning_max_iter) {
    res <- rbind(res, fp_rate(new_b(res)))
    print(res[.N])
  }
  res <- add_setup_info(res[which.min(abs(diff))])
  fwrite(res, "./results/penalties.csv", append = TRUE)
}

#' @export
get_tuned_penalty <- function(data = init_data(mu = 0), params = mvcapa_params(),
                              tuning = tuning_params(), seed = NA) {

  pen <- fread("./results/penalties.csv")
  pen <- pen[n == data$n &
               p == data$p &
               rho == data$rho &
               precision_type == data$precision_type &
               band == data$band &
               block_size == data$block_size &
               cost == params$cost &
               precision_est_struct == params$precision_est_struct &
               minsl == params$minsl &
               maxsl == params$maxsl &
               alpha == tuning$alpha &
               alpha_tol == tuning$alpha_tol &
               tuning_n_sim == tuning$tuning_n_sim]
  if (is.na(params$est_band)) pen <- pen[is.na(est_band)]
  else pen <- pen[est_band == params$est_band]
  if (nrow(pen) == 0) {
    message("The penalty for this setup has not been tuned yet.")
    tune_penalty(data, params, tuning, seed)
    return(get_tuned_penalty(data, params, tuning, seed))
  } else return(pen)
}

#' @export
tune_many <- function(n = 100, p = 10, band = 2) {
  init_b <- c(0.5, 1, 4)
  cor_mat_types <- c("banded", "block_banded", "lattice")
  costs <- c("iid", "cor")
  rhos = c(-0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99)
  lapply(cor_mat_types, function(cor_mat_type) {
    lapply(rhos, function(rho) {
      lapply(costs, function(cost) {
        if (cor_mat_type == "block_banded") m <- p - 1
        else m <- p
        data <- init_data(n = n, p = p, cor_mat_type = cor_mat_type,
                          band = band, rho = rho, block_size = m)
        get_tuned_penalty(data, mvcapa_params(cost), tuning_params(init_b = init_b))
      })
    })
  })
}
