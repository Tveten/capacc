
tpfp_rate <- function(anom_list, tol, data) {
  true_anoms <- data.frame("start" = data$locations,
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
  return(data.frame("fp" = fp, "tp" = tp))
}

#' @export
curve_params <- function(max_dist = 0.2, max_iter = 50, n_sim = 100,
                         init_values = c(0.1, 5)) {
  list("curve_max_dist" = max_dist,
       "curve_max_iter" = max_iter,
       "curve_n_sim"    = n_sim,
       "init_values"    = init_values)
}

read_single_curve <- function(data = init_data(), params = mvcapa_params(),
                              tuning = tuning_params(), curve = curve_params(),
                              .loc_tol = 10) {
  file_name <- "./results/power.csv"
  res <- fread(file_name)
  res[n == data$n &
        p == data$p &
        rho == data$rho &
        precision_type == data$precision_type &
        band == data$band &
        block_size == data$block_size &
        is_equal(proportions, data$proportions) &
        locations == data$locations &
        durations == data$durations &
        change_type == data$change_type &
        precision_est_struct == params$precision_est_struct &
        est_band == params$est_band &
        minsl == params$minsl &
        maxsl == params$maxsl &
        is_equal(alpha, tuning$alpha) &
        is_equal(alpha_tol, tuning$alpha_tol) &
        tuning_n_sim == tuning$tuning_n_sim &
        curve_n_sim == curve$curve_n_sim &
        curve_max_dist == curve$curve_max_dist &
        loc_tol == .loc_tol]
}

already_estimated <- function(data, params, tuning, curve, loc_tol) {
  res <- read_single_curve(data, params, tuning, curve, loc_tol)
  res <- res[cost == params$cost]
  if (nrow(res) > 0) {
    message("The power for this setup has already been estimated. Exiting.")
    return(TRUE)
  } else return(FALSE)
}

#' @export
power_curve <- function(data = init_data(), params = mvcapa_params(),
                        tuning = tuning_params(), curve = curve_params(),
                        loc_tol = 10, seed = NA) {
  add_setup_info <- function(res_dt) {
    res <- cbind(res,
                 as.data.table(data[!(grepl("Sigma", names(data)) |
                                        names(data) == "changing_vars")]),
                 as.data.table(params[names(params) != "b"]),
                 as.data.table(tuning[names(tuning) != "init_b"]),
                 as.data.table(curve[names(curve) != "init_values"]))
    res$loc_tol <- loc_tol
    res$seed <-  seed
    res
  }

  est_power <- function(vartheta) {
    power <- mean(unlist(lapply(1:curve$curve_n_sim, function(i) {
      # cat('.')
      data$mu = vartheta / sqrt(round(data$proportions * data$p))
      tpfp_rate(simulate_mvcapa(data, params, return_anom_only = TRUE),
                loc_tol, data)$tp
    })))
      # cat('\n')
    data.table("vartheta" = vartheta, "power" = power)
  }

  split_inds <- function(res) {
    # exceeds_max_dist <- adjacent_dist(res[, .(vartheta, power)]) > curve$curve_max_dist
    exceeds_max_dist <- diff(res$power) > curve$curve_max_dist
    ind <- which(exceeds_max_dist) + 1
    ind <- ind[!(res[ind - 1]$power >= 0.98 & res[ind]$power >= 0.98)]
    ind <- ind[!(res[ind - 1]$power <= 0.02 & res[ind]$power <= 0.02)]
    ind
  }

  message(paste0("Estimating power curve for n=", data$n,
                 ", p=", data$p,
                 ", cost=", params$cost,
                 ", precision=", data$precision_type,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", prop=", data$proportions,
                 ", precision_est_struct=", params$precision_est_struct,
                 ", est_band=", params$est_band,
                 "."))
  if (already_estimated(data, params, tuning, curve, loc_tol)) return(NULL)

  params$b <- get_tuned_penalty(data, params, tuning, seed = sample(1:1000))$b
  if (!is.na(seed)) set.seed(seed)
  res <- do.call('rbind', Map(est_power, curve$init_values))
  ind <- split_inds(res)
  while (length(ind) > 0 && nrow(res) <= curve$curve_max_iter) {
    new_varthetas <- (res$vartheta[ind - 1] + res$vartheta[ind]) / 2
    res <- rbind(res, do.call('rbind', Map(est_power, new_varthetas)))
    res <- res[order(vartheta)]
    ind <- split_inds(res)
  }
  print(res)
  res <- add_setup_info(res)
  fwrite(res, "./results/power.csv", append = TRUE)
}

#' @export
many_power_curves <- function(ns = 100, ps = 10, costs = c("iid", "cor"),
                              precision_types = "banded", bands = 2,
                              rhos = 0.9, props = 0.1, loc = 50,
                              precision_est_structs = "correct", est_bands = NA,
                              tuning = tuning_params(), curve = curve_params(),
                              loc_tol = 10) {
  lapply(ns, function(n) {
    lapply(ps, function(p) {
      lapply(props, function(prop) {
        lapply(precision_types, function(precision_type) {
          lapply(rhos, function(rho) {
            lapply(bands, function(band) {
              lapply(precision_est_structs, function(precision_est_struct) {
                lapply(est_bands, function(est_band) {
                  seed <- sample(1:10^6, 1)
                  lapply(costs, function(cost) {
                    if (precision_type == "block_banded") {
                      m <- p - 1
                      change_type <- "custom"
                      changing_vars <- p
                    } else {
                      m <- p
                      change_type <- "adjacent"
                      changing_vars <- NA
                    }
                    data <- init_data(n = n,
                                      p = p,
                                      precision_type = precision_type,
                                      band = band,
                                      rho = rho,
                                      block_size = m,
                                      proportions = prop,
                                      locations = loc,
                                      change_type = change_type,
                                      changing_vars = changing_vars)
                    params <- mvcapa_params(cost,
                                            precision_est_struct = precision_est_struct,
                                            est_band = est_band)
                    power_curve(data, params, tuning, curve, loc_tol, seed)
                  })
                })
              })
            })
          })
        })
      })
    })
  })
  NULL
}

#' @export
power_runs <- function() {
  curve <- curve_params(max_dist = 0.1, n_sim = 300)
  many_power_curves(precision_types = c("banded"),
                    rhos = c(-0.3, 0.7, 0.9, 0.99), props = c(0.1, 0.3, 1))
  many_power_curves(precision_types = c("block_banded"),
                    rhos = c(0.5, 0.7, 0.9, 0.99), props = 0.1)
  many_power_curves(precision_types = c("lattice"), rhos = c(0.7, 0.9, 0.99),
                    ns = 200, p = 20, props = c(0.1, 0.3, 1), loc = 100,
                    precision_est_structs = c("correct", "banded"), est_bands = 2)
}

plot_power_curve <- function(data = init_data(), params = mvcapa_params(),
                             tuning = tuning_params(), curve = curve_params(),
                             loc_tol = 10, seed = NA) {
    get_title <- function() {
      if (data$precision_type == "lattice")
        title <- "lattice"
      if (data$precision_type == "banded")
        title <- paste0(data$band, "-banded")
      if (data$precision_type == "block_banded")
        title <- paste0(data$band, "-banded, m=", data$block_size)

      title <- paste0(title, " rho=", data$rho, ", p=", data$p, ", n=", data$n,
             ", s=", data$locations, ", e=", data$locations + data$durations,
             ", prop=", round(data$proportions, 2))
  }

  res <- read_single_curve(data, params, tuning, curve, loc_tol)
  ggplot2::ggplot(data = res, ggplot2::aes(x = vartheta, y = power, colour = cost)) +
    ggplot2::geom_line() +
    ggplot2::ggtitle(get_title()) +
    ggplot2::scale_x_continuous("Signal strength") +
    ggplot2::scale_y_continuous("Power")
}
