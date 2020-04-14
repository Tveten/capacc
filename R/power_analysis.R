
tpfp_rate <- function(anom_list) {
  true_anoms <- data.frame("start" = change_setup$locations,
                           "end"   = change_setup$locations + change_setup$durations)
  est_anoms <- data.frame("start" = unique(anom_list$collective$start),
                          "end"   = unique(anom_list$collective$end))
  n_anom <- nrow(true_anoms)
  if (is.na(est_anoms$start[1])) {
    n_est_anom <- 0
    n_tp <- 0
    tp <- 0
  } else {
    n_est_anom <- nrow(est_anoms)
    n_tp <- 0
    for (i in 1:nrow(true_anoms)) {
      correct_start <- is_in_interval(est_anoms$start, true_anoms$start[i] + c(- tol, tol))
      correct_end <- is_in_interval(est_anoms$end, true_anoms$end[i] + c(- tol, tol))
      which_correct <- which(correct_start & correct_end)
      if (length(which_correct) > 0) {
        n_tp <- n_tp + 1
        if (which_correct[1] < nrow(est_anoms))
          est_anoms <- est_anoms[(which_correct[1] + 1):nrow(est_anoms), ]
        else break
      }
    }
    tp <- n_tp / n_anom
  }
  n_fp <- n_est_anom - n_tp
  min_seg_len <- 2
  n_negatives <- (change_setup$n  - n_anom * min_seg_len) / min_seg_len
  fp <- n_fp / n_negatives
  return(data.frame("fp" = fp, "tp" = tp))
}

curve_params <- function(max_dist = 0.2, max_iter = 50, n_sim = 100,
                         init_values = c(0.)) {
  list("curve_max_dist" = max_dist,
       "curve_max_iter" = max_iter,
       "curve_n_sim"    = n_sim,
       "init_value"     = init_values)
}

power_curve <- function(data = init_data(), params = mvcapa_params(),
                        tuning = tuning_params(), curve = curve_params(),
                        loc_tol = 10, seed = NA) {
  add_setup_info <- function(res_dt) {
    res <- cbind(res,
                 as.data.table(data[!grepl("Sigma", names(data))]),
                 as.data.table(params[names(params) != "b"]),
                 as.data.table(tuning[names(tuning) != "init_b"]),
                 as.data.table(curve))
    res$loc_tol <- loc_tol
    res$seed <-  seed
    res
  }

  est_power <- function(mu) {
    power <- mean(unlist(lapply(1:n_sim, function(i) {
      params$mu <- mu
      tpfp_rate(simulate_mvcapa(data, params, return_anom_only = TRUE))$tp
    })))
    data.table("mu" = mu, "power" = power)
  }

  split_inds <- function(res) {
    exceeds_max_dist <- adjacent_dist(res[, .(FP_rate, TP_rate)]) > max_dist
    ind <- which(exceeds_max_dist) + 1
    ind <- ind[!(res[ind - 1]$TP_rate >= 0.99 & res[ind]$TP_rate >= 0.99)]
    ind <- ind[!(res[ind - 1]$FP_rate <= 0.01 & res[ind]$FP_rate <= 0.01)]
    ind
  }

  params$b <- get_tuned_penalty(data, params, tuning, seed = sample(1:1000), 1)$b
  if (!is.na(seed)) set.seed(seed)
  bs <- sort(c(10^(-4), 1, 10, exp(c(-3, 0.5))))
  res <- do.call('rbind', Map(est_power, mus))
  ind <- split_inds(res)
  while (length(ind) > 0 && nrow(res) <= max_iter) {
    # print(res)
    new_bs <- exp((log(res$b[ind - 1]) + log(res$b[ind])) / 2)
    res <- rbind(res, do.call('rbind', Map(est_power, new_mus)))
    res <- res[order(b)]
    ind <- split_inds(res)
  }
  add_setup_info(res)
}
