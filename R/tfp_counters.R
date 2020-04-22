count_tfp_anom <- function(anom_list, tol, data) {
  true_anoms <- data.frame("start" = data$locations + 1,
                           "end"   = data$locations + data$durations)
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
  n_negatives <- (data$n  - n_anom * min_seg_len) / min_seg_len
  fp <- n_fp / n_negatives
  return(data.frame("fp" = fp, "tp" = tp, "n_fp" = n_fp, "n_tp" = n_tp))
}

