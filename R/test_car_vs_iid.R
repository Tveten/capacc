init_setup <- function(n = 10^3, p = 4, proportions = 0, mu = 1,
                       locations = n - durations - 1, durations = 20,
                       change_type = 'adjacent') {
  return(list('n'           = n,
              'p'           = p,
              'proportions' = proportions,
              'mu'          = mu,
              'locations'   = locations,
              'durations'   = durations,
              'change_type' = change_type))
}

mvcapa_1anom <- function(x, A, b = 1, l = 5, M = nrow(x), cost_type = 'car') {
  get_savings_func <- function(cost_type) {
    if (cost_type == 'iid') return(penalised_savings_iid)
    if (cost_type == 'ar1') return(penalised_savings_ar1)
    if (cost_type == 'car') return(penalised_savings_car)
  }

  if (M < l) stop('M (maximum segment length) must be greater than or equal to l (minimum segment length).')

  p <- ncol(x)
  n <- nrow(x)

  penalised_savings <- get_savings_func(cost_type)
  candidate_cps <- (n - M):(n - l)
  C <- rep(0, n + 1)
  J <- matrix(0, ncol = p, nrow = n)
  for (s in candidate_cps) {
    curr_x <- x[(s + 1):n, , drop = FALSE]
    savings_res <- penalised_savings(curr_x, n, A, b)
    C[s + 1] <- savings_res$B_max
    J[s + 1, ] <- savings_res$J_max
  }

  if (all(C <= 0)) max_s <- NA
  else max_s <- which.max(C) - 1

  # list('max' = max(C), 'max_s' = max_s, 'J_seq' = J)
  list('max' = max(C), 'max_s' = max_s)
}

simulate_mvcapa_1anom <- function(setup = init_setup(200, 10, proportions = 0),
                                  cost_type = 'car', rho = 0.5,
                                  cor_mat_type = 'lattice_car',
                                  b = 1, l = 5, M = nrow(sim_data)) {
  get_Sigma <- function(cor_mat_type) {
    if (cor_mat_type == 'iid')
      return(list('mat'     = diag(1, setup$p),
                  'inverse' = diag(1, setup$p)))
    if (cor_mat_type == 'ar1')
      return(list('mat'     = ar_cor_mat(setup$p, rho),
                  'inverse' = ar_precision_mat(setup$p, rho)))
    if (cor_mat_type == 'lattice_car') {
      precision_mat <- car_precision_mat(lattice_neighbours(setup$p), rho)
      return(list('mat'     = solve(precision_mat),
                  'inverse' = precision_mat))
    }
  }

  Sigma <- get_Sigma(cor_mat_type)
  sim_data <- simulate_cor(n = setup$n, p = setup$p, mu = setup$mu, Sigma = Sigma$mat,
                           locations = setup$locations, durations = setup$durations,
                           proportions = setup$proportions, change_type = setup$change_type)
  mvcapa_1anom(sim_data, Sigma$inverse, b, l, M, cost_type)
}

p_false <- function(cost_type, setup = init_setup(200, 10),
                    cor_mat_type = 'lattice_car', rho = 0.5, n_sim = 100,
                    run_in_parallel = FALSE) {

  get_bs <- function(cost_type, rho) {
    if (cost_type == 'iid' || abs(rho) <= 0.5) {
      return(seq(0.5, 2, length.out = 40))
    } else if (cost_type == 'ar1') {
      scaling_factor <- (1 + rho^2) / (1 - rho^2)
      return(seq(scaling_factor / 2, 1.2 * scaling_factor, length.out = 50))
    } else if (cost_type == 'car') {
      scaling_factor <- 1 / (1 - rho^8)
      return(seq(scaling_factor / 2.5, 1.1 * scaling_factor, length.out = 40))
    }
  }

  bs <- get_bs(cost_type, rho)
  n_b <- length(bs)

  I_false <- matrix(0, ncol = n_b, nrow = n_sim)

  if (run_in_parallel) {
    comp_cluster <- setup_parallel()
    `%dopar%` <- foreach::`%dopar%`
    res_mat <- foreach::foreach(i = 1:length(bs),
                                     .export = c('simulate_mvcapa_1anom',
                                                 'ar_cor_mat', 'ar_precision_mat',
                                                 'car_precision_mat', 'adjacency_mat',
                                                 'lattice_neighbours',
                                                 'simulate_cor', 'mvcapa_1anom',
                                                 'penalised_savings_iid',
                                                 'penalised_savings_car'),
                                     .combine = 'rbind') %dopar% {
      time_start <- proc.time()[3]
      I_false <- rep(0, n_sim)
      for (j in 1:n_sim) {
        mvcapa_sim_res <- simulate_mvcapa_1anom(setup, cost_type, rho, cor_mat_type,
                                                b = bs[i], l = 5, M = 50)
        I_false[j] <- as.numeric(mvcapa_sim_res$max > 0)
      }
      mean_runtime <- (proc.time()[3] - time_start) / n_sim
      c(bs[i], mean(I_false), mean_runtime)
    }
    stop_parallel(comp_cluster)
  } else {
    mean_runtimes <- rep(0, n_b)
    res_mat <- matrix(NA, nrow = n_b, ncol = 3)
    for (i in 1:n_b) {
      time_start <- proc.time()[3]
      I_false <- rep(0, n_sim)
      for (j in 1:n_sim) {
        mvcapa_sim_res <- simulate_mvcapa_1anom(setup, cost_type, rho, cor_mat_type,
                                                b = bs[i], l = 5, M = 50)
        I_false[j] <- as.numeric(mvcapa_sim_res$max > 0)
      }
      mean_runtime <- (proc.time()[3] - time_start) / n_sim
      res_mat[i, ] <- c(bs[i], mean(I_false), mean_runtime)
    }
  }

  res_df <- as.data.frame(res_mat)
  names(res_df) <- c('b', 'p_false', 'mean_runtime')

  return(cbind(res_df, data.frame('n' = rep(setup$n, n_b), 'p' = rep(setup$p, n_b),
                                  'cost_type' = rep(cost_type, n_b),
                                  'cor_mat_type' = rep(cor_mat_type, n_b),
                                  'rho' = rep(rho, n_b), 'n_sim' = rep(n_sim, n_b))))
}

get_p_false_file_name <- function(n, p, cor_mat_type) {
  paste0('p_false_results_', cor_mat_type, '_n', n, 'p', p, '.RData')
}

run_and_save_p_false_dt <- function(n = 200, p = 10, n_sim = 300,
                                    cor_mat_type = 'lattice_car',
                                    run_in_parallel = TRUE) {
  cost_types <- c('iid', 'car')
  rhos <- c(-0.35, - 0.2, - 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99)

  p_false_list <- list(p_false('iid', cor_mat_type = 'iid', rho = 0, n_sim = n_sim))
  for (rho in rhos) {
    for (cost_type in cost_types) {
      print(c(rho, cost_type))
      next_i <- length(p_false_list) + 1
      curr_time <- round(proc.time()[3] / 60, 2)
      p_false_list[[next_i]] <- p_false(cost_type, setup = init_setup(n, p),
                                        cor_mat_type = cor_mat_type,
                                        rho = rho, n_sim = n_sim,
                                        run_in_parallel = run_in_parallel)
      print(p_false_list[[length(p_false_list)]])
      print(paste0('Time used on current iteration: ',
                   round(proc.time()[3] / 60 - curr_time, 2), ' min.'))
    }
  }
  p_false_dt <- data.table::as.data.table(do.call('rbind', p_false_list))
  file_name <- get_p_false_file_name(n, p, cor_mat_type)
  save(p_false_dt, file = file_name)
}

get_penalty_scales <- function(alpha, n, p, cor_mat_type = 'lattice_car') {
  file_name <- get_p_false_file_name(n, p, cor_mat_type)
  load(file_name)
  # print(p_false_dt[, .(cost_type, mean_runtime)])
  bs <- p_false_dt[, b[which.min(abs(p_false - alpha))], by = .(cost_type, rho)]
  names(bs)[3] <- 'b'
  bs$diff <- p_false_dt[, min(abs(p_false - alpha)), by = .(cost_type, rho)]$V1
  tol <- 0.03
  if (any(bs$diff > tol))
    sapply(bs[diff > tol, diff], function(d) warning(paste0('p_false - alpha = ', d)))
  bs
}

run_1anom_test <- function(setup = init_setup(200, 10, 0.3, 1, change_type = 'adjacent'),
                           cor_mat_type = 'lattice_car', alpha = 0.05, n_sim = 100,
                           run_in_parallel = FALSE) {

  penalty_dt <- get_penalty_scales(alpha, setup$n, setup$p, cor_mat_type)
  if (run_in_parallel) {
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
                                                penalty_dt[i, rho], cor_mat_type,
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
  } else {
    res_mat <- matrix(0, ncol = 4, nrow = nrow(penalty_dt))
    for (i in 1:nrow(penalty_dt)) {
      I_detect <- rep(0, n_sim)
      cp_locations <- rep(0, n_sim)
      for (k in 1:n_sim) {
        mvcapa_sim_res <- simulate_mvcapa_1anom(setup, penalty_dt[i, cost_type],
                                                penalty_dt[i, rho], cor_mat_type,
                                                b = penalty_dt[i, b], l = 5, M = 50)
        I_detect[k] <- as.numeric(mvcapa_sim_res$max > 0)
        cp_locations[k] <- mvcapa_sim_res$max_s
      }
      res_mat[i, ] <- c(mean(I_detect),
                        mean(cp_locations, na.rm = TRUE),
                        quantile(cp_locations, 0.025, na.rm = TRUE),
                        quantile(cp_locations, 0.975, na.rm = TRUE))
    }
  }
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

get_res_1anom_file_name <- function(n, p, change_type, cor_mat_type) {
  paste0('res_1anom_', change_type,'_', cor_mat_type, '_n', 200, 'p', 10, '.RData')
}

run_all_1anom_tests <- function(change_type = 'adjacent',
                                cor_mat_type = 'lattice_car', n_sim = 300,
                                props = c(0.1, 0.3, 0.6, 1),
                                mean_changes = c(0.25, 0.5, 0.75, 1, 1.25),
                                append_to_file = FALSE) {
  print(change_type)
  anom_res <- do.call('rbind', lapply(props, function(prop) {
    print(paste0('Proportion = ', prop))
    do.call('rbind', lapply(mean_changes, function(mu) {
      print(paste0('Mean change = ', mu))
      run_1anom_test(init_setup(200, 10, prop, mu, change_type = change_type),
                     cor_mat_type = cor_mat_type, n_sim = n_sim,
                     run_in_parallel = TRUE)
    }))
  }))

  file_name <- get_res_1anom_file_name(anom_res$n[1], anom_res$p[1], change_type, cor_mat_type)
  if (append_to_file) {
    anom_res_new <- anom_res
    load(file_name)
    anom_res <- rbind(anom_res, anom_res_new)
  }
  save(anom_res, file = file_name)
  anom_res
}

plot_1anom_res <- function(change_type = 'adjacent', cor_mat_type = 'lattice_car',
                           proportion = 0.1) {
  file_name <- get_res_1anom_file_name(200, 10, change_type, cor_mat_type)
  load(file_name)
  curr_res_dt <- anom_res[prop == proportion & mu > 0]
  curr_res_dt[, 'rho' := as.factor(rho)]
  ggplot2::ggplot(data = curr_res_dt, ggplot2::aes(x = mu, y = p_detect,
                                                   colour = rho, linetype = cost_type,
                                                   group = interaction(cost_type, rho))) +
    ggplot2::geom_line() +
    ggplot2::ggtitle(paste0(change_type, ', prop = ', proportion))
}

plot_1anom_diff <- function(change_type = 'adjacent', cor_mat_type = 'lattice_car',
                            proportion = 0.1) {
  file_name <- get_res_1anom_file_name(200, 10, change_type, cor_mat_type)
  load(file_name)
  curr_res_dt <- anom_res[prop == proportion & mu > 0]
  curr_res_dt[, 'rho' := as.factor(rho)]
  # print(curr_res_dt[cost_type == 'ar1', .p_detect] - curr_res_dt[cost_type == 'iid'])
  diff_dt <- curr_res_dt[, diff(p_detect), .(rho, mu)]
  names(diff_dt)[3] <- 'p_detect_diff'
  ggplot2::ggplot(data = diff_dt, ggplot2::aes(x = mu, y = p_detect_diff, color = rho)) +
    ggplot2::geom_line() +
    ggplot2::ggtitle(paste0(change_type, ', prop = ', proportion)) +
    ggplot2::coord_cartesian(ylim = c(-0.2, 1))
}

group_plot_1anom <- function(plot_func = plot_1anom_diff,
                             proportions = c(0.1, 0.3, 0.6, 1),
                             cor_mat_type = 'lattice_car') {
  change_type <- c('adjacent', 'scattered')

  plots <- lapply(change_type, function(type) {
    lapply(proportions, plot_func, change_type = type, cor_mat_type = cor_mat_type)
  })
  plots <- c(plots[[1]], plots[[2]])
  gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = plots, nrow = length(change_type)))

}

