####
#### Single anomaly
####
mvcapa_1anom2 <- function(x, A, b = 1, l = 5, M = nrow(x), cost_type = 'ar1') {
  get_savings_func <- function(cost_type) {
    if (cost_type == 'iid') return(penalised_savings_iid)
    if (cost_type == 'ar1') return(penalised_savings_ar1)
  }

  if (M < l) stop('M (maximum segment length) must be greater than or equal to l (minimum segment length).')

  p <- ncol(x)
  n <- nrow(x)

  penalised_savings <- get_savings_func(cost_type)
  candidate_cps <- (n - M):(n - l)
  C <- rep(0, n + 1)
  J <- matrix(0, ncol = p, nrow = n)
  S <- list()
  for (s in candidate_cps) {
    savings_res <- penalised_savings(x[(s + 1):n, , drop = FALSE], n, A, b)
    C[s + 1] <- savings_res$B_max
    J[s + 1, savings_res$J_max] <- 1
    S[[s + 1]] <- savings_res$S_obj
  }

  if (all(C <= 0)) max_s <- NA
  else max_s <- which.max(C) - 1

  list('max' = max(C), 'max_s' = max_s, 'J_seq' = J, 'S_obj_seq' = S)
}

simulate_mvcapa_1anom2 <- function(setup = init_setup(200, 10, proportions = 0),
                                  cost_type = 'ar1', phi = 0.5, cor_mat_type = 'ar1',
                                  b = 1, l = 5, M = nrow(sim_data)) {
  get_Sigma <- function(cor_mat_type) {
    if (cor_mat_type == 'iid')
      return(list('mat'     = diag(1, setup$p),
                  'inverse' = diag(1, setup$p)))
    if (cor_mat_type == 'ar1')
      return(list('mat'     = ar_cor_mat(setup$p, phi),
                  'inverse' = ar_precision_mat(setup$p, phi)))
  }

  Sigma <- get_Sigma(cor_mat_type)
  sim_data <- simulate_cor(n = setup$n, p = setup$p, mu = setup$mu, Sigma = Sigma$mat,
                           locations = setup$locations, durations = setup$durations,
                           proportions = setup$proportions, change_type = setup$change_type)
  mvcapa_1anom2(sim_data, Sigma$inverse, b, l, M, cost_type)
}

test_conditions <- function() {
  no_diff <- function(x) {
    diff_x <- apply(x, 2, diff)
    diff_change_vec <- rep(FALSE, nrow(diff_x))
    for (i in 1:nrow(diff_x)) {
      if (all(diff_x[i, ] == 0)) diff_change_vec[i] <- TRUE
    }
    diff_change_vec
  }

  positive_diff <- function(x, eps = 0.1) {
    diff_x <- apply(x, 2, diff)
    condition_met <- rep(FALSE, nrow(diff_x))
    for (i in 1:nrow(diff_x)) {
      if (all(diff_x[i, 1:3] > 0) && all(abs(diff_x[i, 4:10]) < eps))
        condition_met[i] <- TRUE
    }
    condition_met
  }

  res <- simulate_mvcapa_1anom(init_setup(n = 100, p = 10, proportions = 0.3), b = 0.5)
  J <- res$J[1:96, ]
  S_single <- matrix(0, nrow = length(res$S_obj_seq), ncol = ncol(res$J))
  S_inter <- matrix(0, nrow = length(res$S_obj_seq), ncol = ncol(res$J))
  for (i in 1:length(res$S_obj_seq)) {
    S_single[i, ] <- res$S_obj_seq[[i]]$single
    S_inter[i, ] <- res$S_obj_seq[[i]]$interaction
  }

  print(J)

  # S_single_order <- t(apply(S_single, 1, order, decreasing = TRUE))
  # S_inter_order <- t(apply(S_inter[, -1] - S_single[, -1], 1, order, decreasing = TRUE))
  # no_diff_single <- no_diff(S_single_order)
  # no_diff_inter <- no_diff(S_inter_order)

  pdiff_single <- positive_diff(S_single)
  no_diff_J <- no_diff(J)
  print(apply(S_single, 2, diff))

  res <- cbind(no_diff_J, pdiff_single, pdiff_single & !no_diff_J)
  print(J[res[, 3], ])
  print(res)
  sum(res[, 3])
}

p_false2 <- function(cost_type, setup = init_setup(200, 10), cor_mat_type = 'ar1',
                    phi = 0.5, n_sim = 100) {
  get_bs <- function(cost_type, phi) {
    if (cost_type == 'iid' || abs(phi) <= 0.5) return(seq(0.1, 2.5, 0.1))
    else {
      scaling_factor <- (1 + phi^2) / (1 - phi^2)
      return(seq(scaling_factor / 2, 1.2 * scaling_factor, length.out = 50))
    }
  }

  bs <- get_bs(cost_type, phi)
  n_b <- length(bs)

  I_false <- matrix(0, ncol = n_b, nrow = n_sim)
  mean_runtimes <- rep(0, n_b)

  comp_cluster <- setup_parallel()
  `%dopar%` <- foreach::`%dopar%`
  res_mat <- foreach::foreach(i = 1:length(bs),
                                   .export = c('simulate_mvcapa_1anom',
                                               'ar_cor_mat', 'ar_precision_mat',
                                               'simulate_cor', 'mvcapa_1anom',
                                               'penalised_savings_iid',
                                               'penalised_savings_ar1'),
                                   .combine = 'rbind') %dopar% {
    time_start <- proc.time()[3]
    I_false <- rep(0, n_sim)
    for (j in 1:n_sim) {
      I_false[j] <- as.numeric(simulate_mvcapa_1anom(setup, cost_type, phi, cor_mat_type, b = bs[i], l = 5, M = 50)$max > 0)
    }
    mean_runtime <- (proc.time()[3] - time_start) / n_sim
    c(bs[i], mean(I_false), mean_runtime)
  }
  stop_parallel(comp_cluster)

  res_df <- as.data.frame(res_mat)
  names(res_df) <- c('b', 'p_false', 'mean_runtime')

  return(cbind(res_df, data.frame('n' = rep(setup$n, n_b), 'p' = rep(setup$p, n_b),
                                  'cost_type' = rep(cost_type, n_b),
                                  'cor_mat_type' = rep(cor_mat_type, n_b),
                                  'phi' = rep(phi, n_b), 'n_sim' = rep(n_sim, n_b))))
}

get_p_false_file_name2 <- function(n, p, sign_phi) {
  if (sign_phi == 1) return(paste0('p_false_results_positive_n', n, 'p', p, '.RData'))
  if (sign_phi == -1) return(paste0('p_false_results_negative_n', n, 'p', p, '.RData'))
}

run_and_save_p_false_dt2 <- function(n = 200, p = 10, n_sim = 500, sign_phi = 1) {
  cost_types <- c('iid', 'ar1')
  phis <- sign_phi * c(0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99)

  p_false_list <- list(p_false2('iid', cor_mat_type = 'iid', phi = 0, n_sim = n_sim))
  for (phi in phis) {
    for (cost_type in cost_types) {
      print(c(phi, cost_type))
      p_false_list[[length(p_false_list) + 1]] <- p_false2(cost_type,
                                                           setup = init_setup(n, p),
                                                           phi = phi, n_sim = n_sim)
      print(p_false_list[[length(p_false_list)]])
    }
  }
  p_false_dt <- data.table::as.data.table(do.call('rbind', p_false_list))
  file_name <- get_p_false_file_name(n, p, sign_phi)
  save(p_false_dt, file = file_name)
}

get_penalty_scales2 <- function(alpha, n, p, sign_phi = 1) {
  file_name <- get_p_false_file_name2(n, p, sign_phi)
  load(file_name)
  bs <- p_false_dt[, b[which.min(abs(p_false - alpha))], by = .(cost_type, phi)]
  names(bs)[3] <- 'b'
  bs$diff <- p_false_dt[, min(abs(p_false - alpha)), by = .(cost_type, phi)]$V1
  if (any(bs$diff > 0.02))
    sapply(bs[diff > 0.02, diff], function(d) warning(paste0('p_false - alpha = ', d)))
  bs
}

run_1anom_test2 <- function(setup = init_setup(200, 10, 0.3, 1, durations = 20,
                                              change_type = 'adjacent'),
                           sign_phi = 1, alpha = 0.05, n_sim = 100) {

  penalty_dt <- get_penalty_scales(alpha, setup$n, setup$p, sign_phi)
  comp_cluster <- setup_parallel()
  `%dopar%` <- foreach::`%dopar%`
  res_mat <- foreach::foreach(i = 1:nrow(penalty_dt),
                                   .export = c('simulate_mvcapa_1anom',
                                               'ar_cor_mat', 'ar_precision_mat',
                                               'simulate_cor', 'mvcapa_1anom',
                                               'penalised_savings_iid',
                                               'penalised_savings_ar1'),
                                   .combine = 'rbind') %dopar% {
    I_detect <- rep(0, n_sim)
    cp_locations <- rep(0, n_sim)
    for (k in 1:n_sim) {
      mvcapa_sim_res <- simulate_mvcapa_1anom(setup, penalty_dt[i, cost_type],
                                              penalty_dt[i, phi],
                                              b = penalty_dt[i, b], l = 5, M = 50)
      I_detect[k] <- as.numeric(mvcapa_sim_res$max > 0)
      cp_locations[k] <- mvcapa_sim_res$max_s
    }
    c(mean(I_detect),
      mean(cp_locations, na.rm = TRUE),
      quantile(cp_locations, 0.025, na.rm = TRUE),
      quantile(cp_locations, 0.975, na.rm = TRUE))
  }
  stop_parallel(comp_cluster)
  res_df <- as.data.frame(res_mat)
  names(res_df) <- c('p_detect', 'cp_est', 'cp_lower_ci', 'cp_upper_ci')

  penalty_dt <- cbind(penalty_dt, res_df)
  penalty_dt$n <- rep(setup$n, nrow(penalty_dt))
  penalty_dt$p <- rep(setup$p, nrow(penalty_dt))
  penalty_dt$prop <- rep(setup$proportions, nrow(penalty_dt))
  penalty_dt$mu <- rep(setup$mu, nrow(penalty_dt))
  penalty_dt$dur <- rep(setup$durations, nrow(penalty_dt))
  penalty_dt$change_type <- rep(setup$change_type, nrow(penalty_dt))
  penalty_dt$n_sim <- rep(n_sim, nrow(penalty_dt))

  penalty_dt
}

get_res_1anom_file_name2 <- function(n, p, change_type, sign_phi) {
  if (sign_phi == 1)
    return(paste0('res_1anom_', change_type, '_positive_n', 200, 'p', 10, '.RData'))
  if (sign_phi == -1)
    return(paste0('res_1anom_', change_type, '_negative_n', 200, 'p', 10, '.RData'))
}

run_all_1anom_tests2 <- function(change_type = 'adjacent', sign_phi = 1, n_sim = 500) {
  props <- c(0.1, 0.3, 0.6)
  mean_changes <- c(0.25, 0.5, 0.75, 1, 1.25)
  print(change_type)

  anom_res <- do.call('rbind', lapply(props, function(prop) {
    print(paste0('Proportion = ', prop))
    do.call('rbind', lapply(mean_changes, function(mu) {
      print(paste0('Mean change = ', mu))
      run_1anom_test2(init_setup(200, 10, prop, mu, change_type = change_type),
                     sign_phi = sign_phi, n_sim = n_sim)
    }))
  }))
  file_name <- get_res_1anom_file_name(anom_res$n[1], anom_res$p[1], change_type, sign_phi)
  save(anom_res, file = file_name)
  anom_res
}

plot_1anom_res2 <- function(change_type = 'adjacent', sign_phi = 1, proportion = 0.1) {
  file_name <- get_res_1anom_file_name(200, 10, change_type, sign_phi)
  load(file_name)
  curr_res_dt <- anom_res[prop == proportion & mu > 0]
  curr_res_dt[, 'phi' := as.factor(phi)]
  ggplot2::ggplot(data = curr_res_dt, ggplot2::aes(x = mu, y = p_detect,
                                                   colour = phi, linetype = cost_type,
                                                   group = interaction(cost_type, phi))) +
    ggplot2::geom_line() +
    ggplot2::ggtitle(paste0(change_type, ', prop = ', proportion))
}

plot_1anom_diff2 <- function(change_type = 'adjacent', sign_phi = 1, proportion = 0.1) {
  file_name <- get_res_1anom_file_name(200, 10, change_type, sign_phi)
  load(file_name)
  curr_res_dt <- anom_res[prop == proportion & mu > 0]
  curr_res_dt[, 'phi' := as.factor(phi)]
  # print(curr_res_dt[cost_type == 'ar1', .p_detect] - curr_res_dt[cost_type == 'iid'])
  diff_dt <- curr_res_dt[, diff(p_detect), .(phi, mu)]
  names(diff_dt)[3] <- 'p_detect_diff'
  ggplot2::ggplot(data = diff_dt, ggplot2::aes(x = mu, y = p_detect_diff, color = phi)) +
    ggplot2::geom_line() +
    ggplot2::ggtitle(paste0(change_type, ', prop = ', proportion)) +
    ggplot2::coord_cartesian(ylim = c(-0.2, 1))
}

group_plot_1anom2 <- function(plot_func = plot_1anom_diff, sign_phi = 1) {
  proportions <- c(0.1, 0.3, 0.6)
  change_type <- c('adjacent', 'scattered')

  plots <- lapply(change_type, function(type) {
    lapply(proportions, plot_func, change_type = type, sign_phi = sign_phi)
  })
  plots <- c(plots[[1]], plots[[2]])
  gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = plots, nrow = length(change_type)))

}
