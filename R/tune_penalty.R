
tune_penalty <- function(setup = init_data_setup(mu = 0), cost = "cor",
                         Q_est_struct = "correct", est_band = NA,
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
    res$block_size <- setup$block_size
    res$Q_est_struct <- Q_est_struct
    if (Q_est_struct == "correct") res$est_band <- NA
    else res$est_band <- est_band
    res$minsl <- minsl
    res$maxsl <- maxsl
    res$tol <- tol
    res$max_iter <- 50
    res$n_sim <- n_sim
    res$seed <-  seed
    res
  }

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
    b_lower <- as.numeric(res[diff <= 0, min(b)])
    b_upper <- as.numeric(res[diff > 0][which.min(diff), b])
    exp((log(b_lower) + log(b_upper)) / 2)
  }

  if (!is.na(seed)) set.seed(seed)
  setup$mu <- 0
  message(paste0("Tuning ", cost,  " penalty for n = ", setup$n, ", p = ", setup$p,
                 ", cor_mat_type = ", setup$cor_mat_type,
                 ", band = ", setup$band,
                 " and rho = ", setup$rho))
  res <- do.call('rbind', Map(fp_rate, init_b))
  print(res)
  while (!any(abs(res$diff) <= tol + 1e-6) && nrow(res) <= max_iter) {
    res <- rbind(res, fp_rate(new_b(res)))
    print(res[.N])
  }
  res <- add_setup_info(res[which.min(abs(diff))])
  fwrite(res, "./results/penalties.csv", append = TRUE)
}

#' @export
get_tuned_penalty <- function(setup = init_data_setup(mu = 0), .cost = "cor",
                              .Q_est_struct = "correct", .est_band = NA,
                              .minsl = 2, .maxsl = 100, .alpha = 0.05,
                              init_b = c(0.05, 1, 10), .tol = 0.02, .max_iter = 50,
                              .n_sim = 2 * 10^2, seed = NA) {

  pen <- fread("./results/penalties.csv")
  pen <- pen[alpha == .alpha &
               n == setup$n &
               p == setup$p &
               cost == .cost &
               rho == setup$rho &
               cor_mat_type == setup$cor_mat_type &
               band == setup$band &
               block_size == setup$block_size &
               Q_est_struct == .Q_est_struct &
               minsl == .minsl &
               maxsl == .maxsl &
               tol == .tol &
               n_sim == .n_sim]
  if (is.na(.est_band)) pen <- pen[is.na(est_band)]
  else pen <- pen[est_band == .est_band]
  if (nrow(pen) == 0) {
    tune_penalty(setup, .cost, .Q_est_struct, .est_band, .minsl, .maxsl, .alpha,
                 init_b, .tol, .max_iter, .n_sim, seed)
    return(get_tuned_penalty(setup, .cost, .Q_est_struct, .est_band, .minsl,
                             .maxsl, .alpha, init_b, .tol, .max_iter, .n_sim, seed))
  } else return(pen)
}

#' @export
tune_many <- function(n = 100, p = 10, band = 2) {
  init_b <- c(0.5, 1, 4)
  cor_mat_types <- c("banded", "block_banded", "lattice")
  costs <- c("iid", "cor")
  rhos = c(-0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99)
  lapply(cor_mat_types, function(cor_mat_type) {
    lapply(rhos, function(rho) {
      lapply(costs, function(cost) {
        if (cor_mat_type == "block_banded") m <- p - 1
        else m <- p
        setup <- init_data_setup(n = n, p = p, cor_mat_type = cor_mat_type,
                                 band = band, rho = rho, block_size = m)
        get_tuned_penalty(setup, cost, init_b = init_b)
      })
    })
  })
}
