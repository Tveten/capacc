
#' @export
tuning_params <- function(alpha = 0.05, tol = 0.02, max_iter = 50,
                          init_b = c(0.05, 1, 10), n_sim = 200) {
  list("alpha"           = alpha,
       "alpha_tol"       = tol,
       "tuning_max_iter" = max_iter,
       "init_b"          = init_b,
       "tuning_n_sim"    = n_sim)
}

#' @export
tune_penalty <- function(data = init_data(mu = 0), method = method_params(),
                         tuning = tuning_params(), seed = NA) {

  add_setup_info <- function(res) {
    which_data <- !(grepl("Sigma", names(data)) | names(data) == c("changing_vars"))
    res <- cbind(res,
                 as.data.table(data[which_data]),
                 as.data.table(method[names(method) != "b"]),
                 as.data.table(tuning[names(tuning) != "init_b"]))
    res$seed <-  seed
    res
  }

  fp_rate <- function(b) {
    fps <- unlist(lapply(1:tuning$tuning_n_sim, function(i) {
      method$b <- b
      !is.na(simulate_mvcapa(data, method, return_anom_only = TRUE)$collective$start[1])
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
  data$vartheta <- 0
  message(paste0("Tuning ", method$cost,
                 " penalty for n=", data$n,
                 ", p=", data$p,
                 ", precision_type=", data$precision_type,
                 ", block_size=", data$block_size,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", precision_est_struct=", method$precision_est_struct,
                 ", est_band=", method$est_band))
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
get_tuned_penalty <- function(data = init_data(mu = 0), method = method_params(),
                              tuning = tuning_params(), seed = NA) {
  if (method$cost == "cor_exact") method$cost <- "cor"
  if (!file.exists("./results/penalties.csv")) {
    message(paste0("File does not exist. Making file penalties.csv in ./results/"))
    tune_penalty(data, method, tuning, seed)
    return(get_tuned_penalty(data, method, tuning, seed))
  } else {
    res <- read_penalties("penalties", c(data, method, tuning))
    if (nrow(res) == 0) {
      tune_penalty(data, method, tuning, seed)
      return(get_tuned_penalty(data, method, tuning, seed))
    } else return(res)
  }
}

#' @export
tune_many <- function(n = 100, p = 10, band = 2) {
  init_b <- c(0.5, 1, 4)
  precision_types <- c("banded", "block_banded", "lattice")
  costs <- c("iid", "cor")
  rhos = c(0.5, 0.7, 0.8, 0.9, 0.95, 0.99)
  lapply(precision_types, function(precision_type) {
    lapply(rhos, function(rho) {
      lapply(costs, function(cost) {
        if (precision_type == "block_banded") m <- p - 1
        else m <- p
        data <- init_data(n = n, p = p, precision_type = precision_type,
                          band = band, rho = rho, block_size = m)
        get_tuned_penalty(data, method_params(cost), tuning_params(init_b = init_b))
      })
    })
  })
}

plot_penalty <- function(data, param) {

}
