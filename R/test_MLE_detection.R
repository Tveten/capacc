init_setup_MLE <- function(proportions = 0.3, mu = 0.7, n = 100, p = 8,
                           locations = n - durations - 1, durations = 10,
                           change_type = 'adjacent') {
  return(list('n'           = n,
              'p'           = p,
              'proportions' = proportions,
              'mu'          = mu,
              'locations'   = locations,
              'durations'   = durations,
              'change_type' = change_type))
}

mvcapa_1anom_MLE <- function(x, A, b = 1, l = 5, M = nrow(x),
                             cost_type = 'cor_MLE') {
  get_savings_func <- function(cost_type) {
    if (cost_type == 'iid') return(penalised_savings_iid(n, A))
    if (cost_type == 'cor_MLE') return(optim_penalised_savings_BF(n, A, mu_MLE(A)))
    if (cost_type == 'cor_aMLE') return(optim_penalised_savings_BF(n, A, mu_aMLE()))
  }

  if (M < l) stop('M (maximum segment length) must be greater than or equal to l (minimum segment length).')

  p <- ncol(x)
  n <- nrow(x)

  penalised_savings <- get_savings_func(cost_type)
  candidate_cps <- (n - M):(n - l)
  C <- rep(0, n + 1)
  u <- matrix(0, ncol = p, nrow = n)
  for (s in candidate_cps) {
    curr_x <- x[(s + 1):n, , drop = FALSE]
    savings_res <- penalised_savings(curr_x, b)
    C[s + 1] <- savings_res$B_max
    u[s + 1, ] <- savings_res$u_max
  }

  if (all(C <= 0)) max_s <- NA
  else max_s <- which.max(C) - 1

  # list('max' = max(C), 'max_s' = max_s, 'u_seq' = u)
  list('max' = max(C), 'max_s' = max_s)
}

simulate_mvcapa_1anom_MLE <- function(setup = init_setup_MLE(),
                                      cost_type = 'cor_MLE', rho = 0.5,
                                      cor_mat_type = 'lattice_car',
                                      b = 1, l = 5, M = 50) {
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
  sim_data <- simulate_cor(n = setup$n, p = setup$p, vartheta = vartheta_from_mu(setup$mu),
                           Sigma = Sigma$mat,
                           locations = setup$locations, durations = setup$durations,
                           proportions = setup$proportions, change_type = setup$change_type)
  mvcapa_1anom_MLE(sim_data, Sigma$inverse, b, l, M, cost_type)
}

true_false_positives <- function(cost_type, change_setup = init_setup_MLE(),
                                 rho = 0.5, root_seeds = 1:200,
                                 run_in_parallel = FALSE) {
  add_setup_info <- function(res_dt) {
    res_dt$cost_type <- rep(cost_type, nrow(res_dt))
    res_dt$rho <- rep(rho, nrow(res_dt))
    res_dt$n <- rep(change_setup$n, nrow(res_dt))
    res_dt$p <- rep(change_setup$p, nrow(res_dt))
    res_dt$mu <- rep(change_setup$mu, nrow(res_dt))
    res_dt$durations <- rep(change_setup$durations, nrow(res_dt))
    res_dt$proportions <- rep(change_setup$proportions, nrow(res_dt))
    res_dt$change_type <- rep(change_setup$change_type, nrow(res_dt))
    res_dt
  }

  get_bs <- function(cost_type) {
    # if(cost_type == 'iid') return(seq(0.1, 2, length.out = 16))
    # else {
    #   if (rho <= 0.8) return(seq(0.1, 2, length.out = 16))
    #   else return(c(seq(0.1, 2, length.out = 8), seq(2.2, 4, length.out = 8)))
    # }
    seq(0.1, 2, length.out = 16)
  }

  sim_TP_FP <- function(b, seeds) {
    split_point <- round(length(seeds) / 2)
    FP_seeds <- seeds[1:split_point]
    FP <- vapply(FP_seeds, function(i) {
      set.seed(i)
      no_change_setup <- init_setup_MLE(0)
      cp_est <- simulate_mvcapa_1anom_MLE(no_change_setup, cost_type, rho, b = b)$max_s
      !is.na(cp_est)
    }, logical(1))
    FP_rate <- sum(FP, na.rm = TRUE) / length(FP_seeds)

    TP_seeds <- seeds[(split_point + 1):length(seeds)]
    TP <- vapply(TP_seeds, function(i) {
      set.seed(i)
      cp_est <- simulate_mvcapa_1anom_MLE(change_setup, cost_type, rho, b = b)$max_s
      is_in_interval(cp_est, c(70, 100))
    }, logical(1))
    TP_rate <- sum(TP, na.rm = TRUE) / length(TP_seeds)
    c(b, FP_rate, TP_rate)
  }

  bs <- get_bs(cost_type)
  seeds <- lapply(seq.int(1, 10^6, length.out = length(bs)),
                  function(i) i + root_seeds)
  if (run_in_parallel) {
    comp_cluster <- setup_parallel()
    `%dopar%` <- foreach::`%dopar%`
    res <- foreach::foreach(i = 1:length(bs),
                            .export = c('simulate_mvcapa_1anom_MLE',
                                        'ar_cor_mat', 'ar_precision_mat',
                                        'car_precision_mat', 'adjacency_mat',
                                        'lattice_neighbours',
                                        'simulate_cor', 'mvcapa_1anom_MLE',
                                        'penalised_savings_iid',
                                        'optim_penalised_savings_BF')) %dopar% {
       sim_TP_FP(bs[[i]], seeds[[i]])
    }
    stop_parallel(comp_cluster)
    res <- do.call('rbind', res)
  } else {
    res <- do.call('rbind', Map(sim_TP_FP, bs, seeds))
  }
  res_dt <- as.data.table(res)
  colnames(res_dt) <- c('b', 'FP', 'TP')
  add_setup_info(res_dt)
}

run_MLE_sim <- function(proportions = 2/8, rhos = 0.7, n_sim = 100) {
  cost_types <- c('iid', 'cor_MLE', 'cor_aMLE')
  res <- lapply(proportions, function(prop) {
    rho_res <- lapply(rhos, function(rho) {
      cost_res <- lapply(cost_types, function(cost_type) {
        print(paste0('Prop = ', prop, ', rho = ', rho, ', cost = ', cost_type))
        change_setup <- init_setup_MLE(prop)
        start_seed <- sample(1:10^6, 1)
        root_seeds <- (start_seed + 1):(start_seed + 2 * n_sim)
        true_false_positives(cost_type, change_setup, rho, root_seeds, TRUE)
      })
      do.call('rbind', cost_res)
    })
    do.call('rbind', rho_res)
  })
  res <- do.call('rbind', res)
  save(res, file = get_tfp_file_name())
}

get_tfp_file_name <- function() {
  'roc_mle.RData'
}

# plot_roc <- function(prop = 0.3, correlation = 0.7) {
#   load(get_tfp_file_name())
#   sub_res <- res[proportions == prop & rho == correlation]
#   p <- sub_res$p[1]
#   title <- paste0('prop = ', round(prop * p), '/8',
#                   ', rho = ', correlation)
#   ggplot2::ggplot(sub_res, ggplot2::aes(FP, TP, colour = cost_type)) +
#     ggplot2::geom_line() +
#     ggplot2::ggtitle(title)
# }
#
# plot_all_roc <- function() {
#   props <- c(1/8, 2/8, 4/8)
#   cors <- c(0.7, 0.9, 0.99)
#   prop_cor_combs <- expand.grid(props, cors)
#   plots <- Map(plot_roc, prop_cor_combs[, 1], prop_cor_combs[, 2])
#   # gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = plots,
#   #                                                nrow = length(cors),
#   #                                                ncol = length(props)))
#   ggpubr::ggarrange(plotlist = plots, nrow = length(cors), ncol = length(props),
#                     common.legend = TRUE, legend = "right")
#
# }

aMLE_vs_iid <- function(setup = init_setup_MLE(durations = 98), rho = 0.5,
                        cor_mat_type = 'lattice_car',
                        b = 1, l = 5, M = 50) {
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
  sim_data <- simulate_cor(n = setup$n, p = setup$p,
                           vartheta = vartheta_from_mu(setup$mu), Sigma = Sigma$mat,
                           locations = setup$locations, durations = setup$durations,
                           proportions = setup$proportions, change_type = setup$change_type)
  A <- Sigma$inverse
  mean_x <- colMeans(sim_data)
  P_J <- power_set(setup$p)
  res_diff <- vapply(P_J, function(J) {
    mean_x_J <- rep(0, setup$p)
    mean_x_J[J] <- mean_x[J]
    (2 * mean_x - mean_x_J) %*% A %*% mean_x_J - mean_x_J %*% mean_x_J
  }, numeric(1))
  res_iid <- vapply(P_J, function(J) {
    mean_x_J <- rep(0, setup$p)
    mean_x_J[J] <- mean_x[J]
    98 * mean_x_J %*% mean_x_J - 3.4 * length(J)
  }, numeric(1))
  res_aMLE <- vapply(P_J, function(J) {
    mean_x_J <- rep(0, setup$p)
    mean_x_J[J] <- mean_x[J]
    98 * (2 * mean_x - mean_x_J) %*% A %*% mean_x_J - 3.4 * length(J)
  }, numeric(1))
  mu_MLE_func <- mu_MLE(A)
  res_MLE <- vapply(P_J, function(J) {
    mu_hat <- rep(0, setup$p)
    mu_hat[J] <- mu_MLE_func(mean_x, J)
    98 * (2 * mean_x - mu_hat) %*% A %*% mu_hat - 3.4 * length(J)
  }, numeric(1))
  print(res_diff)
  print(c(max(res_aMLE), max(res_MLE), max(res_iid)))
  sum(res_diff <= 0)
}
