count_tfp_anom <- function(anom_list, tol, data) {
  true_anoms <- data.frame("start" = data$locations + 1,
                           "end"   = data$locations + data$durations)
  est_anoms <- data.frame("start" = unique(anom_list$collective$start),
                          "end"   = unique(anom_list$collective$end))
  n_anom <- nrow(true_anoms)
  n_est_anom <- nrow(est_anoms)
  if (n_est_anom > 0) {
    n_est_anom <- 0
    n_tp <- 0
    tp <- 0
  } else {
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

count_collective_anomalies <- function(anom_list) {
  length(unique(anom_list$collective$start))
}

rand_counts_mvcapa <- function(anom_list, data) {
  true_anoms <- data.frame("start" = data$locations + 1,
                           "end"   = data$locations + data$durations)
  est_anoms <- data.frame("start" = unique(anom_list$collective$start),
                          "end"   = unique(anom_list$collective$end))
}

label_anom_est <- function(anom_list, data) {
  labs <- rep(0, data$n)
  starts <- unique(anom_list$collective$start)
  ends <- unique(anom_list$collective$end)
  anom_inds <- integer(0)
  if (length(starts) > 0)
    anom_inds <- c(anom_inds, unlist(lapply(1:length(starts), function(i) {
      starts[i]:ends[i]
    })))
  anom_inds <- c(anom_inds, anom_list$point$location)
  if (length(anom_inds) > 0)
    labs[anom_inds] <- 1
  labs
}

true_labels <- function(data) {
  labs <- rep(0, data$n)
  anom_inds <- c(unlist(lapply(1:length(data$locations), function(i) {
    (data$locations[i] + 1):(data$locations[i] + data$durations[i])
  })), data$point_locations)
  labs[anom_inds] <- 1
  labs
}

count_binary_tfp <- function(x, y) {
  tp <- sum(x & y)
  tn <- sum(!x & !y)
  fp <- sum(x & !y)
  fn <- sum(!x & y)
  list(tp = tp, tn = tn, fp = fp, fn = fn)
}

count_rand_tfp <- function(x, y) {
  # x: Estimated classification labels.
  # y: True classification labels.

  equal_label <- function(u) {
    force(u)
    function(i, j) {
      u[i] == u[j]
    }
  }

  n <- length(unlist(x))
  inds <- combn(1:n, 2)
  equal_lab_x <- unlist(Map(equal_label(x), i = inds[1, ], j = inds[2, ]))
  equal_lab_y <- unlist(Map(equal_label(y), i = inds[1, ], j = inds[2, ]))
  tp <- sum(equal_lab_x & equal_lab_y)
  tn <- sum(!equal_lab_x & !equal_lab_y)
  fp <- sum(equal_lab_x & !equal_lab_y)
  fn <- sum(!equal_lab_x & equal_lab_y)
  list(tp_rand = tp, tn_rand = tn, fp_rand = fp, fn_rand = fn)
}

adjusted_rand_index <- function(x, y) {
  stopifnot(length(x) == length(y))
  x_labs <- unique(x)
  y_labs <- unique(y)
  n <- matrix(0, nrow = length(x_labs), ncol = length(y_labs))
  for (i in seq_along(x_labs)) {
    for (j in seq_along(y_labs)) {
      n[i, j] <- sum(x == x_labs[i] & y == y_labs[j])
    }
  }
  a <- rowSums(n)
  b <- colSums(n)
  m <- length(x)
  ab_sum <- sum(choose(a, 2)) + sum(choose(b, 2))
  ab_prod <- sum(choose(a, 2)) * sum(choose(b, 2)) / choose(m, 2)
  (sum(choose(n, 2)) - ab_prod) / (ab_sum / 2 - ab_prod)
}


perf_metric_texts <- function(postfix = NULL) {
  list(tp = paste0("tp", postfix),
       fp = paste0("fp", postfix),
       tn = paste0("tn", postfix),
       fn = paste0("fn", postfix))
}

precision <- function(a, postfix = NULL) {
  e <- perf_metric_texts(postfix)
  if ((a[[e$tp]] + a[[e$fp]]) > 0) return(a[[e$tp]] / (a[[e$tp]] + a[[e$fp]]))
  else return(0)
}

recall <- function(a, postfix = NULL) {
  e <- perf_metric_texts(postfix)
  if (a[[e$tp]] + a[[e$fn]] > 0) return(a[[e$tp]] / (a[[e$tp]] + a[[e$fn]]))
  else return(0)
}

acc <- function(a, postfix = NULL) {
  e <- perf_metric_texts(postfix)
  (a[[e$tp]] + a[[e$tn]]) / (a[[e$tp]] + a[[e$tn]] + a[[e$fn]] + a[[e$fp]])
}

bacc <- function(a, postfix = NULL) {
  e <- perf_metric_texts(postfix)
  (a[[e$tp]] / (a[[e$tp]] + a[[e$fn]]) + a[[e$tn]] / (a[[e$tn]] + a[[e$fp]])) / 2
}



