sim_roc <- function(data = init_data(), params = method_params(),
                    tuning = tuning_params(),
                    curve = curve_params(init_values = sort(c(10^(-4), 1, 10, exp(c(-3, 0.5))))),
                    loc_tol = 10, seed = NA) {
  add_setup_info <- function(res) {
    res <- cbind(res,
                 as.data.table(data[!(grepl("Sigma", names(data)) |
                                        names(data) == "changing_vars")]),
                 as.data.table(params[names(params) != "b"]),
                 as.data.table(tuning[names(tuning) != "init_b"]),
                 as.data.table(curve[names(curve) != "init_values"]))
    res$vartheta <- res$mu / sqrt(round(data$proportions * data$p))
    res$loc_tol <- loc_tol
    res$seed <-  seed
    res
  }

  sim_tpfp <- function(b) {
    cat('.')
    tpfp_df <- do.call("rbind", lapply(1:curve$curve_n_sim, function(i) {
      params$b <- b
      anom_list <- simulate_mvcapa(data, params, return_anom_only = TRUE)
      tpfp_rate(anom_list, loc_tol, data)
    }))
    data.table("b" = b, "tp" = mean(tpfp_df$tp), "fp" = mean(tpfp_df$fp))
  }

  split_inds <- function(res) {
    exceeds_max_dist <- adjacent_dist(res[, .(fp, tp)]) > curve$curve_max_dist
    ind <- which(exceeds_max_dist) + 1
    ind <- ind[!(res[ind - 1]$tp >= 0.99 & res[ind]$fp >= 0.99)]
    ind <- ind[!(res[ind - 1]$tp <= 0.01 & res[ind]$fp <= 0.01)]
    ind
  }

  message(paste0("Estimating roc for n=", data$n,
                 ", p=", data$p,
                 ", cost=", params$cost,
                 ", precision=", data$precision_type,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", prop=", data$proportions,
                 ", mu=", data$mu,
                 ", precision_est_struct=", params$precision_est_struct,
                 ", est_band=", params$est_band,
                 "."))

  res <- do.call('rbind', Map(sim_tpfp, curve$init_values))
  ind <- split_inds(res)
  while (length(ind) > 0 && nrow(res) <= curve$curve_max_iter) {
    new_bs <- exp((log(res$b[ind - 1]) + log(res$b[ind])) / 2)
    res <- rbind(res, do.call('rbind', Map(sim_tpfp, new_bs)))
    res <- res[order(b)]
    ind <- split_inds(res)
  }
  cat('\n')
  print(res)
  fwrite(add_setup_info(res), "./results/roc.csv", append = TRUE)
}

many_rocs <- function(data = init_data(), params = method_params(),
                      costs = c("iid", "cor"), bands = 2, rhos = 0.9,
                      varthetas = 1, props = 0.1, precision_est_structs = "correct",
                      est_bands = NA, tuning = tuning_params(),
                      curve = curve_params(), loc_tol = 10) {

}

roc_runs <- function() {

}

# run_mvcapa_sim <- function(out_file, n = 500, p = 20, rhos = c(0.5, 0.7, 0.9, 0.99),
#                            band = 2,
#                            varthetas = 1, mus = 1, locations = n/2, durations = 10,
#                            proportions = c(1/p, round(p^{-5/8}, 2), 1),
#                            cor_mat_type = "banded", change_type = "adjacent",
#                            random_changes = TRUE, nb_struct = "true", est_band = 2,
#                            max_iter = 50, max_dist = 0.2, n_sim = 100, tol = 10,
#                            init_seed = 4) {
#   cost_types <- c('iid', 'cor')
#   force(p)
#   set.seed(init_seed)
#   start_time <- proc.time()[3]
#   res <- do.call("rbind", lapply(varthetas, function(vartheta) {
#     do.call("rbind", lapply(proportions, function(prop) {
#       do.call("rbind", lapply(rhos, function(rho) {
#         start_seed <- sample(1:10^6, 1)
#         root_seeds <- (start_seed + 1):(start_seed + n_sim)
#         do.call("rbind", lapply(cost_types, function(cost_type) {
#           print(paste0("vartheta = ", vartheta, ", prop = ", prop,  ', rho = ', rho, ', cost = ', cost_type))
#           print(paste0("Time elapsed: ", round((proc.time()[3] - start_time) / 60, 1), " min."))
#           change_setup <- init_data(n = n,
#                                     p = p,
#                                     proportions = prop,
#                                     mu = vartheta / sqrt(round(prop * p)),
#                                     locations = locations,
#                                     durations = durations,
#                                     rho = rho,
#                                     band = band,
#                                     change_type = change_type,
#                                     cor_mat_type = cor_mat_type)
#           sim_roc_curve(cost_type, change_setup, random_changes,
#                         nb_struct, est_band, tol, max_iter, max_dist,
#                         root_seeds = root_seeds)
#         }))
#       }))
#     }))
#   }))
#   fwrite(res, file = paste0("./results/", out_file), append = TRUE)
#   res
# }
#
# runs <- function() {
#   run_mvcapa_sim("mvcapa_roc_p10n100dense.csv", n = 100, p = 10, rhos = 0.99, varthetas = 1.5, locations = 50, durations = 10, proportions = 1, change_type = "adjacent", cor_mat_type = "banded", random_changes = FALSE, n_sim = 1000, init_seed = 86)
#   run_mvcapa_sim("mvcapa_roc_p10n100dense_noscaling.csv", n = 100, p = 10, rhos = 0.99, varthetas = 3, locations = 50, durations = 10, proportions = 1, change_type = "adjacent", cor_mat_type = "banded", random_changes = FALSE, n_sim = 1000, init_seed = 86)
#   run_mvcapa_sim("mvcapa_roc_p10n100dense.csv", n = 100, p = 10, rhos = 0.7, varthetas = 1.5, locations = 50, durations = 10, proportions = 1, change_type = "adjacent", cor_mat_type = "banded", random_changes = FALSE, n_sim = 1000, init_seed = 54)
#
#   rhos <- c(0.5, 0.7, 0.9, 0.99)
#   init_seeds <- c(19, 50, 22, 43)
#   start_time <- proc.time()[3]
#   for (i in 1:length(rhos)) {
#     run_mvcapa_sim("mvcapa_roc_multiple_p20n300.csv", n = 300, p = 20, rhos = rhos[i], varthetas = seq(0.5, 1.5, 0.25), locations = c(50, 100, 200), durations = c(5, 10, 20), proportions = c(1, 4, 20)/20, change_type = "adjacent", cor_mat_type = "banded", random_changes = FALSE, n_sim = 100, init_seed = init_seeds[i])
#   }
#
#   rho <- 0.99
#   start_time <- proc.time()[3]
#   run_mvcapa_sim(n = 100, p = 20, rhos = rho, mus = seq(0.25, 1.3, 0.25), locations = 50, durations = 10, proportions = c(1, 4, 10)/20, change_type = "adjacent", cor_mat_type = "lattice", random_changes = FALSE, n_sim = 100)
#   rhos <- c(0.5, 0.7, 0.9, 0.99)
#   start_time <- proc.time()[3]
#   for (rho in rhos) {
#     run_mvcapa_sim(n = 100, p = 20, rhos = rho, mus = seq(0.25, 1.3, 0.25), locations = 50, durations = 10, proportions = c(1, 4, 10)/20, change_type = "adjacent", cor_mat_type = "lattice", nb_struct = "banded", est_band = 2, random_changes = FALSE, n_sim = 100)
#     # run_mvcapa_sim(n = 100, p = 10, rhos = rho, mus = seq(0.5, 1.3, 0.25), locations = 50, durations = 10, proportions = c(1, 2, 10)/10, change_type = "adjacent", random_changes = FALSE, n_sim = 100)
#   }
#
#   rhos <- c(0.5, 0.7, 0.9, 0.99)
#   start_time <- proc.time()[3]
#   for (rho in rhos) {
#     run_mvcapa_sim(n = 100, p = 20, rhos = rho, mus = seq(0.25, 1.3, 0.25), locations = 50, durations = 10, proportions = c(1, 4, 10)/20, change_type = "adjacent", cor_mat_type = "lattice", random_changes = FALSE, n_sim = 100)
#     # run_mvcapa_sim(n = 100, p = 10, rhos = rho, mus = seq(0.5, 1.3, 0.25), locations = 50, durations = 10, proportions = c(1, 2, 10)/10, change_type = "adjacent", random_changes = FALSE, n_sim = 100)
#   }
#   run_mvcapa_sim(n = 100, p = 10, rhos = 0.99, mus = seq(0.5, 1.3, 0.1), locations = 50, durations = 10, proportions = c(1, 2)/10, change_type = "scattered", random_changes = FALSE, n_sim = 100)
#   run_mvcapa_sim(n = 100, p = 10, rhos = 0.99, mus = seq(0.5, 1.3, 0.1), locations = 50, durations = 10, proportions = c(1, 2)/10, change_type = "adjacent", random_changes = FALSE, n_sim = 100)
#   run_mvcapa_sim(n = 100, p = 50, rhos = 0.99, mus = seq(0.5, 1.3, 0.1), locations = 50, durations = 10, change_type = "adjacent", random_changes = FALSE, n_sim = 100)
# }

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



sim_fptp <- function(b, cost_type, change_setup = init_data(),
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
