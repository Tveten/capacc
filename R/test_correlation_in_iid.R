# DONE: Try other penalties by adjusting psi = a log(n).
# DONE: Create function for running simulation study with different a.
### DONE: Make sure within detectability boundaries.

# init_setup <- function(n = 10^3, p = 4, proportions = 0, mu = 1,
#                        locations = n - durations - 1, durations = 20,
#                        change_type = 'adjacent') {
#   return(list('n'           = n,
#               'p'           = p,
#               'proportions' = proportions,
#               'mu'          = mu,
#               'locations'   = locations,
#               'durations'   = durations,
#               'change_type' = change_type))
# }

adjusted_penalty <- function(a, n, p, C = 2) {
  a <- as.numeric(a)
  psi <- a * log(n)
  # k_star <- sqrt(p) * psi / log(p)
  k_star <- (p + C * sqrt(p * psi)) / (2 * log(p))
  sum_beta <- C * psi + C * 1:p * log(p)
  sum_beta[1:p > k_star] <- p + C * psi + C * sqrt(p * psi)
  diff(c(0, sum_beta))
}

penalty_func <- function(b, a, n, p, C = 2) {
  penalty <- linear_penalty(b, 2, n, p)
  beta_linear <- penalty$beta
  k_star <- floor(penalty$k_star)
  beta_const <- const_penalty(b, 2, n, p)
  beta <- rep(0, p)
  beta[1:k_star] <- 1:k_star * beta_linear
  beta[(k_star + 1):p] <- beta_const
  list('sum_beta' = beta, 'k_star' = k_star)
}

signal_strength2 <- function(mu, proportion, p, n, a) {
  psi <- a * log(n)
  k_star <- sqrt(p) * psi / log(p)
  size_J <- round(proportion * p)
  squared_norm_mu <- mu^2 * size_J

  if(size_J <= k_star) return(squared_norm_mu / (log(p) + psi / size_J))
  if(size_J > k_star) return(squared_norm_mu / (sqrt(p * psi) / size_J + psi / size_J))
}

min_duration <- function(a = 4, n = 10^4, p = 4, proportion = 0.5, mu = 1, C = 2) {
  delta_squared <- signal_strength2(mu, proportion, p, n, a)
  40 * C / delta_squared
}





simulate_mvcapa_iid <- function(setup = init_setup(10^3, 4), Sigma = diag(1, setup$p), a = 'default') {
  sim_data <- simulate_cor(n = setup$n, p = setup$p, vartheta = vartheta_from_mu(setup$mu),
                           Sigma = Sigma,
                           locations = setup$locations, durations = setup$durations,
                           proportions = setup$proportions)
  if (a != 'default') beta <- adjusted_penalty(a, n = setup$n, p = setup$p)
  else beta <- NULL

  anomaly::capa.mv(sim_data, type = 'mean', beta = beta,
                   max_seg_len = min(1500, max(100, 2 * setup$durations)))
}

rK_hat <- function(n_sim = 1, setup = init_setup(10^3, 4), Sigma = diag(1, setup$p), a = 'default') {
  vapply(1:n_sim, function(i) {
      res <- simulate_mvcapa(setup = setup, Sigma = Sigma, a = a)
      c_anomalies <- anomaly::collective_anomalies(res)
      if (is.na(c_anomalies$start[1])) return(0)
      else return(length(unique(c_anomalies$start)))
    },
    numeric(1))
}

rK_hat_for <- function(param_name, params, setup = init_setup(10^3, 4),
                       n_sim = 10^3, parallelise = TRUE) {
  valid_param_name <- function() {
    return(any(c('a', 'rho', 'phi', 'alphad') == param_name))
  }

  get_Sigma <- function(param_name, param) {
    if (param_name == 'a') return(diag(1, setup$p))
    if (param_name == 'rho') return(constant_cor_mat(setup$p, param))
    if (param_name == 'phi') return(ar_cor_mat(setup$p, param))
    if (param_name == 'alphad') return(tpca::rcor_mat(setup$p, alphad = param))
  }

  get_a <- function(param_name, param) {
    if (param_name == 'a') return(param)
    else {
      if (setup$proportions > 0) return(3)
      else return(2)
    }
  }

  valid_param_name()
  print(c(param_name, params))
  print(proc.time()[3])

  if (parallelise) {
    comp_cluster <- setup_parallel()
    `%dopar%` <- foreach::`%dopar%`
    K_hat_list_of_lists <- foreach::foreach(param = params,
                                            .export = c('ar_cor_mat',
                                                        'constant_cor_mat',
                                                        'adjusted_penalty',
                                                        'simulate_mvcapa',
                                                        'rK_hat')) %dopar% {
      list(param, rK_hat(n_sim, setup,
                         Sigma = get_Sigma(param_name, param),
                         a     = get_a(param_name, param))
      )
    }
    stop_parallel(comp_cluster)
  }
  else {
    K_hat_list_of_lists <- lapply(params, function(param) {
      list(param, rK_hat(n_sim, setup,
                         Sigma = get_Sigma(param_name, param),
                         a     = get_a(param_name, param))
      )
    })
  }

  K_hat_list <- extract_nested_element(2, K_hat_list_of_lists)
  names(K_hat_list) <- unlist(extract_nested_element(1, K_hat_list_of_lists))
  K_hat_dt <- data.table::as.data.table(reshape2::melt(K_hat_list, value.name = 'K_hat'))
  names(K_hat_dt)[2] <- 'param_value'
  K_hat_dt[, 'param' := rep(param_name, .N)]
  K_hat_dt
}

run_sim <- function(run_name, setup, n_sim = 10^3, parallelise = TRUE) {
  # WARNING: Assumes that 'full_K_hat_list' exists in the global environment.
  # If K > 0, remember to check min_duration().
  get_as <- function(prop) {
    if(prop > 0) return(c('default', 2, 2.5, 3))
    else return(c('default', 1, 1.5, 2))
  }

  print(c(run_name, proc.time()[3]))
  K_hat <- setup
  K_hat$n_sim <- n_sim

  as <- get_as(setup$proportions)
  correlation_params <- c(0.3, 0.45, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)
  param_names = c('a', 'rho', 'phi')
  params = list(as, correlation_params, correlation_params)
  K_hat_dt_list <- Map(rK_hat_for, param_names, params,
                       MoreArgs = list('setup' = setup, 'n_sim' = n_sim, 'parallelise' = parallelise))
  K_hat$dt <- do.call('rbind', K_hat_dt_list)

  all_K_hat_results[[run_name]] <<- K_hat
  save(all_K_hat_results, file = 'current_results.RData')
}

run_all <- function() {
  # p = 4
  run_sim('n1e3p4prop0',  init_setup(n = 10^3, p = 4, proportions = 0))
  run_sim('n1e4p4prop05', init_setup(n = 10^4, p = 4, proportions = 0.5, durations = 610))

  # p = 100
  run_sim('n1e3p100prop0',   init_setup(n = 10^3,     p = 100, proportions = 0))
  run_sim('n3e3p100prop003', init_setup(n = 3 * 10^3, p = 100, proportions = 0.03,
                                        durations = 280), n_sim = 100)
}

plot_K_hat <- function(K_hat_list, xlim = c(0, 7), ylim = c(0, n_sim)) {
  # K_hat_list: Output from run_sim
  get_title <- function(param_name) {
    # title <- paste0('P(K_hat = k) for different ', param_name)
    title <- 'P(K_hat = k)'
    if (param_name == 'a')
      title <- paste0(title, ', Sigma = I')
    title <- paste0(title,
                    '. n = ', K_hat_list$n,
                    ', p = ', K_hat_list$p,
                    ', prop = ', K_hat_list$proportions)
    if (K_hat_list$proportions > 0)
      title <- paste0(title, ', dur = ', K_hat_list$durations)
    title
  }

  get_color <- function(n, param_name) {
    if (param_name == 'a') return(RColorBrewer::brewer.pal(2 * n, 'OrRd')[(n + 1):(2 * n)])
    else return(rainbow(n))
  }

  plot_K_hat_pmf <- function(param_name) {
    K_hat_dt <- K_hat_list$dt[param == param_name]
    K_hat_dt[K_hat > xlim[2], K_hat := xlim[2]]
    ggplot2::ggplot(data = K_hat_dt, ggplot2::aes(K_hat, fill = param_value)) +
      ggplot2::geom_bar(alpha = 1, width = 0.7, position = ggplot2::position_dodge()) +
      ggplot2::ggtitle(get_title(param_name)) +
      ggplot2::labs(fill = param_name) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
      ggplot2::scale_fill_manual(values = get_color(length(unique(K_hat_dt$param_value)), param_name))
  }

  n_sim <- K_hat_list$n_sim
  param_names <- unique(K_hat_list$dt$param)
  pmf_plots <- lapply(param_names, plot_K_hat_pmf)
  gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs = pmf_plots, nrow = length(pmf_plots)))
}
