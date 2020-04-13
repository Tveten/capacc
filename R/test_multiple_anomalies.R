#' @export
init_data_setup <- function(n = 200, p = 10, proportions = sqrt(p)/p, mu = 1,
                            locations = n - durations - 1, durations = 10,
                            change_type = 'adjacent', point_locations = NA,
                            point_proportions = NA, point_mu = NA,
                            cor_mat_type = 'banded', rho = 0.8, band = 2,
                            block_size = round(p/2),
                            min_nbs = 1, max_nbs = 3) {
  get_Sigma <- function(cor_mat_type) {
    if (cor_mat_type == 'iid') {
      return(list('mat'     = diag(1, p),
                  'inverse' = diag(1, p)))
    } else if (cor_mat_type == 'ar1') {
      return(list('mat'     = ar_cor_mat(p, rho),
                  'inverse' = ar_precision_mat(p, rho)))
    } else if (cor_mat_type == 'lattice') {
      precision_mat <- car_precision_mat(lattice_neighbours(p), rho)
      return(list('mat'     = solve(precision_mat),
                  'inverse' = precision_mat))
    } else if (cor_mat_type == 'banded') {
      precision_mat <- car_precision_mat(banded_neighbours(band, p), rho)
      return(list('mat'     = solve(precision_mat),
                  'inverse' = precision_mat))
    } else if (cor_mat_type == 'random') {
      precision_mat <- car_precision_mat(list('random', p), rho,
                                         min_nbs = min_nbs, max_nbs = max_nbs)
      return(list('mat'     = solve(precision_mat),
                  'inverse' = precision_mat))
    } else if (cor_mat_type == "block_banded") {
      precision_mat <- block_precision_mat(p, block_size, within_block_type = "banded",
                                           rho, band = band)
      return(list('mat'     = solve(precision_mat),
                  'inverse' = precision_mat))
    }
  }

  Sigma_obj <- get_Sigma(cor_mat_type)
  list('n'                 = n,
       'p'                 = p,
       'mu'                = mu,
       'cor_mat_type'      = cor_mat_type,
       'rho'               = rho,
       'band'              = band,
       'block_size'        = block_size,
       'Sigma'             = Sigma_obj$mat,
       'Sigma_inv'         = Sigma_obj$inverse,
       'locations'         = locations,
       'durations'         = durations,
       'proportions'       = proportions,
       'change_type'       = change_type,
       'point_locations'   = point_locations,
       'point_proportions' = point_proportions,
       'point_mu'          = point_mu)
}

simulate_mvcapa <- function(setup = init_data_setup(), cost_type = "cor", b = 1,
                            seed = NULL, min_seg_len = 2,
                            max_seg_len = 100, prune = TRUE,
                            return_anom_only = FALSE, nb_struct = "correct", est_band = 2) {
  if (!is.null(seed)) set.seed(seed)
  x <- simulate_cor(n                 = setup$n,
                    p                 = setup$p,
                    mu                = setup$mu,
                    Sigma             = setup$Sigma,
                    locations         = setup$locations,
                    durations         = setup$durations,
                    proportions       = setup$proportions,
                    change_type       = setup$change_type,
                    point_locations   = setup$point_locations,
                    point_proportions = setup$point_proportions,
                    point_mu          = setup$point_mu)
  x <- anomaly::robustscale(x)
  if (cost_type == "cor") {
    if (nb_struct == "correct")
      Q_hat <- estimate_precision_mat(x, setup$Sigma_inv)
    else if (nb_struct == "banded")
      Q_hat <- estimate_precision_mat(x, adjacency_mat(banded_neighbours(est_band, setup$p)))
    res <- mvcapa_cor(x, Q_hat,
                      b                = b,
                      b_point          = max(0.05, b),
                      min_seg_len      = min_seg_len,
                      max_seg_len      = max_seg_len)
    if (!return_anom_only) return(res)
    else return(list("collective" = collective_anomalies(list("anoms" = res)),
                     "point"      = point_anomalies(list("anoms" = res))))
  } else if (cost_type == "iid") {
    beta <- iid_penalty(setup$n, setup$p, b)
    beta_tilde <- iid_point_penalty(setup$n, setup$p, max(0.05, b))
    res <- anomaly::capa.mv(x,
                            beta = beta,
                            beta_tilde = beta_tilde,
                            min_seg_len = min_seg_len,
                            max_seg_len = max_seg_len,
                            type        = "mean")
    if (!return_anom_only) return(res)
    else return(list("collective" = anomaly::collective_anomalies(res),
                     "point"      = anomaly::point_anomalies(res)))
  }
}

change_sd <- function(p, prop, duration) {
  xi <- - log(prop) / log(p)
  if (xi >= 0 && xi <= 1/2) mu_boundary <- p^(- 1/2 + xi)
  else if (xi > 1/2 && xi <= 3/4) mu_boundary <- sqrt((2 * xi - 1) * log(p))
  else if (xi > 3/4 && xi <= 1) mu_boundary <- sqrt(2 * log(p)) * (1 - sqrt(1 - xi))

  # > pnorm(1/8) - pnorm(-1/8)
  # [1] 0.09947645
  # I.e., 90% of drawn changes will be greater than the detection boundary
  # (when averaging over durations).
  mu_boundary / sqrt(duration)
}

draw_anomaly_setup <- function(data_setup, mean_duration = 10,
                               anomaly_rate = 3 / data_setup$n,
                               change_sd = log(data_setup$p), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  locations <- list(rnbinom(1, 1, anomaly_rate) + 1)
  durations <- list(rpois(1, mean_duration))
  if (locations[[1]] + durations[[1]] >= data_setup$n)
    locations[[1]] <- data_setup$n - durations[[1]] - 1
  mus <- list(rnorm(1, sd = change_sd(data_setup$p, data_setup$proportions, mean_duration)))
  i <- 1
  while (locations[[i]] + durations[[i]] < data_setup$n) {
    locations[[i + 1]] <- locations[[i]] + durations[[i]] + rnbinom(1, 1, anomaly_rate)
    durations[[i + 1]] <- rpois(1, mean_duration)
    mus[[i + 1]] <- rnorm(1, sd = change_sd(data_setup$p, data_setup$proportions, mean_duration))
    i <- i + 1
  }
  n_anoms <- length(locations) - 1
  data_setup$locations <- unlist(locations)[1:n_anoms]
  data_setup$durations <- unlist(durations)[1:n_anoms]
  data_setup$mu <- unlist(mus)[1:n_anoms]
  data_setup
}

sim_roc_curve <- function(cost_type, change_setup = init_data_setup(),
                          random_changes = FALSE, nb_struct = "correct", est_band = 2,
                          tol = 10, max_iter = 50, max_dist = 0.2,
                          root_seeds = 1:100) {
  add_setup_info <- function(res_dt) {
    res_dt$cost_type <- rep(cost_type, nrow(res_dt))
    res_dt$rho <- rep(change_setup$rho, nrow(res_dt))
    res_dt$band <- rep(change_setup$band, nrow(res_dt))
    res_dt$cor_mat_type <- rep(change_setup$cor_mat_type, nrow(res_dt))
    res_dt$n <- rep(change_setup$n, nrow(res_dt))
    res_dt$p <- rep(change_setup$p, nrow(res_dt))
    res_dt$proportion <- rep(change_setup$proportions, nrow(res_dt))
    res_dt$change_type <- rep(change_setup$change_type, nrow(res_dt))
    if (random_changes) {
      res_dt$mu <- rep(NA, nrow(res_dt))
      res_dt$duration <- rep(NA, nrow(res_dt))
      res_dt$location <- rep(NA, nrow(res_dt))
    } else {
      res_dt$mu <- rep(change_setup$mu[1], nrow(res_dt))
      res_dt$vartheta <- rep(res_dt$mu[1] * sqrt(round(change_setup$proportions * change_setup$p)), nrow(res_dt))
      res_dt$duration <- rep(change_setup$durations[1], nrow(res_dt))
      res_dt$location <- rep(change_setup$locations[1], nrow(res_dt))
    }
    res_dt$est_nb_struct <- rep(nb_struct, nrow(res_dt))
    if (nb_struct == "correct") res_dt$est_band <- rep(NA, nrow(res_dt))
    else if (nb_struct == "banded") res_dt$est_band <- rep(est_band, nrow(res_dt))
    res_dt
  }

  count_FP_TP2 <- function(anom_list) {
    true_starts <- change_setup$locations
    true_ends <- true_starts + change_setup$durations
    true_anom_inds <- inds_from_intervals(true_starts, true_ends, change_setup$n)
    true_points <- change_setup$point_locations
    true_anom_inds[true_points] <- TRUE

    est_starts <- unique(anom_list$collective$start)
    est_ends <- unique(anom_list$collective$end)
    est_anom_inds <- inds_from_intervals(est_starts, est_ends, change_setup$n)
    est_points <- unique(anom_list$point$location)
    est_anom_inds[est_points] <- TRUE

    n_tp <- sum(est_anom_inds & true_anom_inds)
    n_fp <- sum(est_anom_inds & !true_anom_inds)

    n_positives <- sum(true_anom_inds)
    n_negatives <- change_setup$n - n_positives
    return(data.frame("FP_rate" = n_fp / n_negatives, "TP_rate" = n_tp / n_positives))
  }

  count_FP_TP <- function(anom_list) {
    true_anoms <- data.frame("start" = change_setup$locations,
                             "end"   = change_setup$locations + change_setup$durations)
    est_anoms <- data.frame("start" = unique(anom_list$collective$start),
                            "end"   = unique(anom_list$collective$end))
    n_anom <- nrow(true_anoms)
    if (is.na(est_anoms$start[1])) {
      n_est_anom <- 0
      n_tp <- 0
      tp_rate <- 0
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
      tp_rate <- n_tp / n_anom
    }
    n_fp <- n_est_anom - n_tp
    min_seg_len <- 2
    n_negatives <- (change_setup$n  - n_anom * min_seg_len) / min_seg_len
    fp_rate <- n_fp / n_negatives
    return(data.frame("FP_rate" = fp_rate, "TP_rate" = tp_rate))
  }

  sim_FP_TP <- function(b, seeds) {
    fp_tp_df <- do.call("rbind", lapply(seeds, function(seed) {
      # cat('.')
      if (random_changes)
        change_setup <- draw_anomaly_setup(change_setup, seed = seed + 333)
      anom_list <- simulate_mvcapa(change_setup, cost_type, b, seed,
                                   return_anom_only = TRUE, nb_struct = nb_struct)
      count_FP_TP(anom_list)
    }))
    # cat('\n')
    data.table("b"       = b,
               "FP_rate" = mean(fp_tp_df$FP_rate),
               "TP_rate" = mean(fp_tp_df$TP_rate))
  }

  split_inds <- function(res) {
    exceeds_max_dist <- adjacent_dist(res[, .(FP_rate, TP_rate)]) > max_dist
    ind <- which(exceeds_max_dist) + 1
    ind <- ind[!(res[ind - 1]$TP_rate >= 0.99 & res[ind]$TP_rate >= 0.99)]
    ind <- ind[!(res[ind - 1]$FP_rate <= 0.01 & res[ind]$FP_rate <= 0.01)]
    ind
  }

  seeds <- lapply(seq.int(1, 10^6, length.out = 2 * max_iter),
                  function(i) i + root_seeds)
  bs <- sort(c(10^(-4), 1, 10, exp(c(-3, 0.5))))
  res <- do.call('rbind', Map(sim_FP_TP, bs, seeds[1:length(bs)]))
  ind <- split_inds(res)
  while (length(ind) > 0 && nrow(res) <= max_iter) {
    # print(res)
    new_bs <- exp((log(res$b[ind - 1]) + log(res$b[ind])) / 2)
    seed_ind <- (nrow(res) + 1):(nrow(res) + length(new_bs))
    res <- rbind(res, do.call('rbind', Map(sim_FP_TP, new_bs, seeds[seed_ind])))
    res <- res[order(b)]
    ind <- split_inds(res)
  }
  add_setup_info(res)
}

sim_fptp <- function(b, cost_type, change_setup = init_data_setup(),
                     random_changes = FALSE, nb_struct = "true", est_band = 2,
                     tol = 10, max_iter = 50, max_dist = 0.2,
                     root_seeds = 1:100) {
  count_FP_TP <- function(anom_list) {
    true_anoms <- data.frame("start" = change_setup$locations,
                             "end"   = change_setup$locations + change_setup$durations)
    est_anoms <- data.frame("start" = unique(anom_list$collective$start),
                            "end"   = unique(anom_list$collective$end))
    n_anom <- nrow(true_anoms)
    if (is.na(est_anoms$start[1])) {
      n_est_anom <- 0
      n_tp <- 0
      tp_rate <- 0
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
      tp_rate <- n_tp / n_anom
    }
    n_fp <- n_est_anom - n_tp
    min_seg_len <- 2
    n_negatives <- (change_setup$n  - n_anom * min_seg_len) / min_seg_len
    fp_rate <- n_fp / n_negatives
    return(data.frame("FP_rate" = fp_rate, "TP_rate" = tp_rate))
  }

  sim_FP_TP <- function(b, seeds) {
    fp_tp_df <- do.call("rbind", lapply(seeds, function(seed) {
      # cat('.')
      if (random_changes)
        change_setup <- draw_anomaly_setup(change_setup, seed = seed + 333)
      anom_list <- simulate_mvcapa(change_setup, cost_type, b, seed,
                                   return_anom_only = TRUE, nb_struct = nb_struct)
      count_FP_TP(anom_list)
    }))
    # cat('\n')
    data.table("b"       = b,
               "FP_rate" = fp_tp_df$FP_rate,
               "TP_rate" = fp_tp_df$TP_rate,
               "cost_type" = cost_type)
  }

  split_inds <- function(res) {
    exceeds_max_dist <- adjacent_dist(res[, .(FP_rate, TP_rate)]) > max_dist
    ind <- which(exceeds_max_dist) + 1
    ind <- ind[!(res[ind - 1]$TP_rate >= 0.99 & res[ind]$TP_rate >= 0.99)]
    ind <- ind[!(res[ind - 1]$FP_rate <= 0.01 & res[ind]$FP_rate <= 0.01)]
    ind
  }

  sim_FP_TP(b, root_seeds)
}

plot_tpfp <- function(res) {
  ggplot2::ggplot(res, ggplot2::aes(TP_rate, fill = cost_type)) +
    ggplot2::geom_histogram(alpha = 1, position = "dodge", binwidth = 0.01)
  # ggplot2::ggplot(res, ggplot2::aes(FP_rate, TP_rate, colour = cost_type)) +
  #   ggplot2::geom_point() +
  #   ggplot2::scale_x_continuous(name = "False positives", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  #   ggplot2::scale_y_continuous(name = "True positives", limits = c(0, 1), breaks = seq(0, 1, 0.2))
}


run_mvcapa_sim <- function(out_file, n = 500, p = 20, rhos = c(0.5, 0.7, 0.9, 0.99),
                           band = 2,
                           varthetas = 1, mus = 1, locations = n/2, durations = 10,
                           proportions = c(1/p, round(p^{-5/8}, 2), 1),
                           cor_mat_type = "banded", change_type = "adjacent",
                           random_changes = TRUE, nb_struct = "true", est_band = 2,
                           max_iter = 50, max_dist = 0.2, n_sim = 100, tol = 10,
                           init_seed = 4) {
  cost_types <- c('iid', 'cor')
  force(p)
  set.seed(init_seed)
  start_time <- proc.time()[3]
  res <- do.call("rbind", lapply(varthetas, function(vartheta) {
    do.call("rbind", lapply(proportions, function(prop) {
      do.call("rbind", lapply(rhos, function(rho) {
        start_seed <- sample(1:10^6, 1)
        root_seeds <- (start_seed + 1):(start_seed + n_sim)
        do.call("rbind", lapply(cost_types, function(cost_type) {
          print(paste0("vartheta = ", vartheta, ", prop = ", prop,  ', rho = ', rho, ', cost = ', cost_type))
          print(paste0("Time elapsed: ", round((proc.time()[3] - start_time) / 60, 1), " min."))
          change_setup <- init_data_setup(n = n,
                                          p = p,
                                          proportions = prop,
                                          mu = vartheta / sqrt(round(prop * p)),
                                          locations = locations,
                                          durations = durations,
                                          rho = rho,
                                          band = band,
                                          change_type = change_type,
                                          cor_mat_type = cor_mat_type)
          sim_roc_curve(cost_type, change_setup, random_changes,
                        nb_struct, est_band, tol, max_iter, max_dist,
                        root_seeds = root_seeds)
        }))
      }))
    }))
  }))
  fwrite(res, file = paste0("./results/", out_file), append = TRUE)
  res
}

runs <- function() {
  run_mvcapa_sim("mvcapa_roc_p10n100dense.csv", n = 100, p = 10, rhos = 0.99, varthetas = 1.5, locations = 50, durations = 10, proportions = 1, change_type = "adjacent", cor_mat_type = "banded", random_changes = FALSE, n_sim = 1000, init_seed = 86)
  run_mvcapa_sim("mvcapa_roc_p10n100dense_noscaling.csv", n = 100, p = 10, rhos = 0.99, varthetas = 3, locations = 50, durations = 10, proportions = 1, change_type = "adjacent", cor_mat_type = "banded", random_changes = FALSE, n_sim = 1000, init_seed = 86)
  run_mvcapa_sim("mvcapa_roc_p10n100dense.csv", n = 100, p = 10, rhos = 0.7, varthetas = 1.5, locations = 50, durations = 10, proportions = 1, change_type = "adjacent", cor_mat_type = "banded", random_changes = FALSE, n_sim = 1000, init_seed = 54)

  rhos <- c(0.5, 0.7, 0.9, 0.99)
  init_seeds <- c(19, 50, 22, 43)
  start_time <- proc.time()[3]
  for (i in 1:length(rhos)) {
    run_mvcapa_sim("mvcapa_roc_multiple_p20n300.csv", n = 300, p = 20, rhos = rhos[i], varthetas = seq(0.5, 1.5, 0.25), locations = c(50, 100, 200), durations = c(5, 10, 20), proportions = c(1, 4, 20)/20, change_type = "adjacent", cor_mat_type = "banded", random_changes = FALSE, n_sim = 100, init_seed = init_seeds[i])
  }

  rho <- 0.99
  start_time <- proc.time()[3]
  run_mvcapa_sim(n = 100, p = 20, rhos = rho, mus = seq(0.25, 1.3, 0.25), locations = 50, durations = 10, proportions = c(1, 4, 10)/20, change_type = "adjacent", cor_mat_type = "lattice", random_changes = FALSE, n_sim = 100)
  rhos <- c(0.5, 0.7, 0.9, 0.99)
  start_time <- proc.time()[3]
  for (rho in rhos) {
    run_mvcapa_sim(n = 100, p = 20, rhos = rho, mus = seq(0.25, 1.3, 0.25), locations = 50, durations = 10, proportions = c(1, 4, 10)/20, change_type = "adjacent", cor_mat_type = "lattice", nb_struct = "banded", est_band = 2, random_changes = FALSE, n_sim = 100)
    # run_mvcapa_sim(n = 100, p = 10, rhos = rho, mus = seq(0.5, 1.3, 0.25), locations = 50, durations = 10, proportions = c(1, 2, 10)/10, change_type = "adjacent", random_changes = FALSE, n_sim = 100)
  }

  rhos <- c(0.5, 0.7, 0.9, 0.99)
  start_time <- proc.time()[3]
  for (rho in rhos) {
    run_mvcapa_sim(n = 100, p = 20, rhos = rho, mus = seq(0.25, 1.3, 0.25), locations = 50, durations = 10, proportions = c(1, 4, 10)/20, change_type = "adjacent", cor_mat_type = "lattice", random_changes = FALSE, n_sim = 100)
    # run_mvcapa_sim(n = 100, p = 10, rhos = rho, mus = seq(0.5, 1.3, 0.25), locations = 50, durations = 10, proportions = c(1, 2, 10)/10, change_type = "adjacent", random_changes = FALSE, n_sim = 100)
  }
  run_mvcapa_sim(n = 100, p = 10, rhos = 0.99, mus = seq(0.5, 1.3, 0.1), locations = 50, durations = 10, proportions = c(1, 2)/10, change_type = "scattered", random_changes = FALSE, n_sim = 100)
  run_mvcapa_sim(n = 100, p = 10, rhos = 0.99, mus = seq(0.5, 1.3, 0.1), locations = 50, durations = 10, proportions = c(1, 2)/10, change_type = "adjacent", random_changes = FALSE, n_sim = 100)
  run_mvcapa_sim(n = 100, p = 50, rhos = 0.99, mus = seq(0.5, 1.3, 0.1), locations = 50, durations = 10, change_type = "adjacent", random_changes = FALSE, n_sim = 100)
}

roc <- function(res, prop = 0.2, m = 1, ss = 1) {
  if (is.na(m)) sub_res <- res[is.na(mu)][proportion == prop]
  else sub_res <- res[mu == m][proportion == prop]
  # title <- paste0('prop = ', round(prop * p), '/10', ', rho = ', r)
  title <- paste0("pr=", round(prop, 2), ", mu=", m)
  ggplot2::ggplot(sub_res, ggplot2::aes(FP_rate, TP_rate, colour = cost_type)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = sub_res[b == 1],
                        ggplot2::aes(FP_rate, TP_rate, colour = cost_type),
                        shape = 4, size = 3) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous(name = "False positives", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(name = "True positives", limits = c(0, 1), breaks = seq(0, 1, 0.2))
}

roc_ss <- function(res, prop = 0.2, ss = 1) {
  sub_res <- res[vartheta == ss][proportion == prop]
  # title <- paste0('prop = ', round(prop * p), '/10', ', rho = ', r)
  title <- paste0("pr=", round(prop, 2), ", ss=", ss)
  ggplot2::ggplot(sub_res, ggplot2::aes(FP_rate, TP_rate, colour = cost_type)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(data = sub_res[b == 1],
                        ggplot2::aes(FP_rate, TP_rate, colour = cost_type),
                        shape = 4, size = 3) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous(name = "False positives", limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(name = "True positives", limits = c(0, 1), breaks = seq(0, 1, 0.2))
}

adjust_res_lattice <- function(res) {
  res <- res[!(cost_type == "iid" & est_nb_struct == "banded")]
  res[, "type" := "iid"]
  res[cost_type == "cor" & est_nb_struct == "true", "type" := "cor_true"]
  res[cost_type == "cor" & est_nb_struct == "banded", "type" := "cor_banded"]
  res[, "cost_type" := type]
  res
}

#' @export
prop_mean_roc <- function(file_name, rho_ = 0.9, mus = NULL, props = NULL,
                          n_obs = 100, n_var = 10, cmt = "banded",
                          band = 2, duration = NA, location = NA,
                          ct = "adjacent", nb_struct = "banded") {
  res <- fread(paste0("./results/", file_name))
  res <- res[rho == rho_]
  if (cmt == "lattice") res <- adjust_res_lattice(res)

  if (is.null(mus)) mus <- unique(res$mu)
  if (is.null(props)) props <- unique(res$proportion)

  mean_prop_combs <- expand.grid(mus, props)
  plots <- Map(roc, prop = mean_prop_combs[, 2], m = mean_prop_combs[, 1],
               MoreArgs = list("res" = res))
  ggpubr::ggarrange(plotlist = plots, nrow = length(props), ncol = length(mus),
                    common.legend = TRUE, legend = "right")
}

#' @export
prop_ss_roc <- function(file_name, rho_ = 0.9, varthetas = NULL, props = NULL,
                        n_obs = 100, n_var = 10, cmt = "banded",
                        band = 2, duration = NA, location = NA,
                        ct = "adjacent", nb_struct = "banded") {
  res <- fread(paste0("./results/", file_name))
  res <- res[rho == rho_]
  if (cmt == "lattice") res <- adjust_res_lattice(res)

  if (is.null(varthetas)) varthetas <- unique(res$vartheta)
  if (is.null(props)) props <- unique(res$proportion)

  mean_prop_combs <- expand.grid(varthetas, props)
  plots <- Map(roc_ss, prop = mean_prop_combs[, 2], ss = mean_prop_combs[, 1],
               MoreArgs = list("res" = res))
  ggpubr::ggarrange(plotlist = plots, nrow = length(props), ncol = length(varthetas),
                    common.legend = TRUE, legend = "right")
}

saving_iid <- function(J, x) {
  means <- colMeans(x[, J])
  nrow(x) * sum(means^2)
}

compare_signal_strength <- function(mu, Q) {
  data.frame("type"            = c("cor", "iid"),
             "signal_strength" = c(sqrt(as.numeric(mu %*% Q %*% mu)),
                                   sqrt(as.numeric(mu %*% mu))))
}

find_threshold <- function(setup = init_data_setup(), alpha = 0.99, n_sim = 10^4) {
  current_thresholds <- fread("./results/alpha_dense.csv")
  current_thresholds <- current_thresholds[n == setup$n & p == setup$p & alpha == a & cor_mat_type == setup$cor_mat_type & band == band & rho == rho]
  if (nrow(current_thresholds) > 0) stop("A threshold for this setup already exists.")

  res_dt <- data.table("value" = rep(0, 2 * n_sim),
                       "type" = c(rep("cor_saving", n_sim),
                                  rep("iid_saving", n_sim)))
  for (i in 1:n_sim) {
    x <- simulate_cor(n = setup$n, p = setup$p, mu = 0, Sigma = setup$Sigma)
    res_dt$value[i] <- dense_mvnormal_savings(matrix(colMeans(x[1:2, ]), nrow = setup$p),
                                              setup$Sigma_inv, 2)
    res_dt$value[n_sim + i] <- saving_iid(1:setup$p, x[1:2, ])
  }
  threshold_dt <- res_dt[, .(threshold    = quantile(value, alpha),
                             alpha        = alpha,
                             n            = setup$n,
                             p            = setup$p,
                             cor_mat_type = setup$cor_mat_type,
                             band         = setup$band,
                             rho          = setup$rho),
                         by = type]
  fwrite(threshold_dt, "./results/alpha_dense.csv", append = TRUE)
}


#' @export
saving_distr <- function(s, e, setup = init_data_setup(n = 100, p = 10, mu = 0, rho = 0.9),
                         b = 1, n_sim = 10^4, a = 0.99, seed = NULL) {
  get_title <- function(s, e) {
    if (is_in_interval(setup$locations, c(s, e)) | is_in_interval(setup$locations + setup$durations, c(s, e))) {
      mu <- setup$mu
      prop <- setup$proportions
    } else {
      mu <- 0
      prop <- 0
    }
    paste0("2-banded, rho=", setup$rho, ", p=", setup$p, ", n=", setup$n,
           ", s=", s, ", e=", e, ", mu=", setup$mu, ", prop=", setup$proportions)
  }

  get_threshold <- function(setup, a) {
    alpha <- fread("./results/alpha_dense.csv")
    alpha <- alpha[n == setup$n & p == setup$p & alpha == a & cor_mat_type == setup$cor_mat_type & band == setup$band & rho == setup$rho]
    if (nrow(alpha) == 0) {
      message("Finding alpha_dense. May take time.")
      find_threshold(setup, alpha = a)
      return(get_threshold(setup, a))
    } else return(alpha)
  }

  if (!is.null(seed)) set.seed(seed)
  res_dt <- data.table("value" = rep(0, 3 * n_sim),
                       "type" = c(rep("cor_saving", n_sim),
                                  rep("iid_saving", n_sim),
                                  rep("chisq_p", n_sim)))
  for (i in 1:n_sim) {
    x <- simulate_cor(n                 = setup$n,
                      p                 = setup$p,
                      mu                = setup$mu,
                      Sigma             = setup$Sigma,
                      locations         = setup$locations,
                      durations         = setup$durations,
                      proportions       = setup$proportions,
                      change_type       = setup$change_type,
                      point_locations   = setup$point_locations,
                      point_proportions = setup$point_proportions,
                      point_mu          = setup$point_mu)
    res_dt$value[i] <- dense_mvnormal_savings(matrix(colMeans(x[s:e, ]), nrow = setup$p),
                                              setup$Sigma_inv, e - s + 1)
    res_dt$value[n_sim + i] <- saving_iid(1:setup$p, x[s:e, ])
  }
  res_dt$value[(2 * n_sim + 1):(3 * n_sim)] <- rchisq(n_sim, setup$p)
  alpha <- get_threshold(setup, a)
  print(paste0("cor detection rate: ", res_dt[type == "cor_saving", sum(value > alpha[type == "cor_saving", threshold]) / n_sim]))
  print(paste0("iid detection rate: ", res_dt[type == "iid_saving", sum(value > alpha[type == "iid_saving", threshold]) / n_sim]))

  ggplot2::ggplot(data = res_dt, ggplot2::aes(value, colour = type)) +
    ggplot2::geom_density() +
    ggplot2::scale_x_continuous("Saving") +
    ggplot2::ggtitle(get_title(s, e)) +
    ggplot2::geom_vline(xintercept = alpha[, threshold],
                        color = c("green", "blue"), size=1.5)
}

#' @export
plot_signal_strength <- function(p, mu, rho, band = 2) {
  Q <- car_precision_mat(banded_neighbours(band, p), rho = rho)
  ss <- do.call("rbind", lapply(1:p, function(i) {
    mu_vec <- rep(0, p)
    mu_vec[1:i] <- mu
    cbind(compare_signal_strength(mu_vec, Q), data.frame("prop" = round(i/p, 4)))
  }))
  ggplot2::ggplot(data = ss, ggplot2::aes(prop, signal_strength, colour = type)) +
    ggplot2::geom_line() +
    ggplot2::scale_y_continuous("Signal strength", limits = c(0, max(ss$signal_strength))) +
    ggplot2::scale_x_continuous("Proportion") +
    ggplot2::ggtitle(paste0("2-banded, p = ", p, ", rho = ", rho, ", mu = ", mu))
}



find_penalty_scale <- function(setup = init_data_setup(), alpha = 0.01, n_sim = 10^4) {
  current_scales <- fread("./results/penalty_scales.csv")
  current_scales <- current_scales[n == setup$n & p == setup$p & alpha == a & cor_mat_type == setup$cor_mat_type & band == setup$band & rho == setup$rho]
  if (nrow(current_scales) > 0) stop("A penalty scale for this setup already exists.")

  res_dt <- data.table("value" = rep(0, 2 * n_sim),
                       "type" = c(rep("cor", n_sim), rep("iid", n_sim)))
  for (i in 1:n_sim) {
    x <- simulate_cor(n = 2, p = setup$p, mu = 0, Sigma = setup$Sigma, locations = 1, durations = 1)
    res_dt$value[i] <- dense_mvnormal_savings(matrix(colMeans(x), nrow = setup$p),
                                              setup$Sigma_inv, 2)
    res_dt$value[n_sim + i] <- saving_iid(1:setup$p, x)
  }
  alpha_dense <- get_penalty("dense", setup$n, setup$p)$alpha
  q <- 1 - alpha / setup$n
  scale_dt <- res_dt[, .(b            = quantile(value, q) / alpha_dense,
                         alpha        = alpha,
                         n            = setup$n,
                         p            = setup$p,
                         cor_mat_type = setup$cor_mat_type,
                         band         = setup$band,
                         rho          = setup$rho),
                     by = type]
  fwrite(scale_dt, "./results/penalty_scales.csv", append = TRUE)
}

find_penalty_scale_sparse <- function(setup = init_data_setup(), alpha = 0.01, n_sim = 10^4) {
  current_scales <- fread("./results/penalty_scales_sparse.csv")
  current_scales <- current_scales[n == setup$n & p == setup$p & alpha == a & cor_mat_type == setup$cor_mat_type & band == setup$band & rho == setup$rho]
  if (nrow(current_scales) > 0) stop("A penalty scale for this setup already exists.")

  res_dt <- data.table("value" = rep(0, 2 * n_sim),
                       "type" = c(rep("cor", n_sim), rep("iid", n_sim)))
  Q <- setup$Sigma_inv
  for (i in 1:n_sim) {
    x <- simulate_cor(n = 2, p = setup$p, mu = 0, Sigma = setup$Sigma, locations = 1, durations = 1)
    mean_x <- colMeans(x)
    res_dt$value[i] <- max(2 * ((2 * mean_x * Q %*% mean_x) - diag(as.matrix(Q)) * mean_x^2))
    res_dt$value[n_sim + i] <- max(2 * mean_x^2)
  }
  penalty <- get_penalty("sparse", setup$n, setup$p)
  q <- 1 - alpha / setup$n
  scale_dt <- res_dt[, .(b            = quantile(value, q) / (penalty$alpha + penalty$beta),
                         alpha        = alpha,
                         n            = setup$n,
                         p            = setup$p,
                         cor_mat_type = setup$cor_mat_type,
                         band         = setup$band,
                         rho          = setup$rho),
                     by = type]
  fwrite(scale_dt, "./results/penalty_scales_sparse.csv", append = TRUE)
}


check_fpr <- function(setup = init_data_setup(), scale_by = "sparse", b = NULL,
                      n_sim = 10^2, a = 0.01, seed = NULL) {
  get_b <- function(setup, a) {
    if (scale_by == "dense") file_name <- "./results/penalty_scales.csv"
    if (scale_by == "sparse") file_name <- "./results/penalty_scales_sparse.csv"
    bs <- fread(file_name)
    bs <- bs[n == setup$n & p == setup$p & alpha == a & cor_mat_type == setup$cor_mat_type & band == setup$band & rho == setup$rho]
    if (nrow(bs) == 0) {
      message("Finding penalty scale. May take time.")
      if (scale_by == "dense") find_penalty_scale(setup, alpha = a)
      if (scale_by == "sparse") find_penalty_scale_sparse(setup, alpha = a)
      return(get_b(setup, a))
    } else return(bs)
  }

  costs <- c("cor", "iid")
  if (is.null(b)) bs <- get_b(setup, a)
  else bs <- data.table("type" = costs, "b" = b)
  print(bs)
  seeds <- sample(1:10^7, n_sim)
  res <- do.call("rbind", lapply(seeds, function(seed) {
    do.call("rbind", lapply(costs, function(cost) {
      sim_res <- simulate_mvcapa(setup,
                                 cost_type        = cost,
                                 b                = bs[type == cost, b],
                                 seed             = seed,
                                 return_anom_only = TRUE)
      data.table("type" = cost, "fp" = !is.na(sim_res$collective$start[1]))
    }))
  }))
  res[, .(fp = mean(fp)), by = type]
}

tune_penalty2 <- function(setup = init_data_setup(), cost = "cor",
                         Q_est_struct = "true", est_band = 2,
                         minsl = 2, maxsl = 100, alpha = 0.05,
                         init_b = c(0.05, 1, 10), tol = 0.02, max_iter = 50,
                         n_sim = 2 * 10^2, seed = NA) {

  add_setup_info <- function(res) {
    res$alpha <- alpha
    res$n <- setup$n
    res$p <- setup$p
    res$cost <- cost
    res$rho <- setup$rho
    res$cor_mat_type <- setup$cor_mat_type
    res$band <- setup$band
    if (Q_est_struct == "true") res$Q_est_struct <- setup$cor_mat_type
    else if (Q_est_struct == "banded") res$Q_est_struct <- Q_est_struct
    res$est_band <- est_band
    res$minsl <- minsl
    res$maxsl <- maxsl
    res$tol <- tol
    res$max_iter <- 50
    res$n_sim <- n_sim
    res$seed <-  seed
    res
  }

  # fp_rate <- function(b, seeds) {
  fp_rate <- function(b) {
    fps <- unlist(lapply(1:n_sim, function(i) {
      # !is.na(simulate_mvcapa(setup, cost = cost, b = b, seed = seed,
      !is.na(simulate_mvcapa(setup, cost = cost, b = b,
                             min_seg_len = minsl, max_seg_len = maxsl,
                             nb_struct = Q_est_struct, est_band = est_band,
                             return_anom_only = TRUE)$collective$start[1])
    }))
    data.table("b" = b, "fp" = mean(fps), "diff" = mean(fps) - alpha)
  }

  new_b <- function(res) {
    res <- res[order(diff)]
    b_lower <- res[diff <= 0, min(b)]
    b_upper <- res[diff > 0][which.min(diff), b]
    exp((log(b_lower) + log(b_upper)) / 2)
  }

  if (!is.na(seed)) set.seed(seed)
  # seeds <- lapply(seq.int(1, 10^6, length.out = 2 * max_iter),
  #                 function(i) i + 1:n_sim)
  # res <- do.call('rbind', Map(fp_rate, bs, seeds[1:length(bs)]))
  res <- do.call('rbind', Map(fp_rate, init_b))
  # ind <- split_inds(res)
  while (all(abs(res$diff) > tol) && nrow(res) <= max_iter) {
    # ind <- split_inds(res)
    # new_bs <- exp((log(res$b[ind - 1]) + log(res$b[ind])) / 2)
    # seed_ind <- (nrow(res) + 1):(nrow(res) + length(new_bs))
    # res <- rbind(res, do.call('rbind', Map(fp_rate, new_bs, seeds[seed_ind])))
    res <- rbind(res, fp_rate(new_b(res)))[order(b)]
    print(abs(res$diff))
    print(tol)
    print(abs(res$diff) > tol)
    print(res)
  }
  print(res)
  add_setup_info(res[which.min(abs(diff))])
}
