#### C++ helpers #####

collective_anomalies <- function(mvcapa_cor_res) {

  if (nrow(mvcapa_cor_res$anoms) == 0)
    return(data.frame('start' = NA, 'end' = NA, 'variate' = NA, 'mean_change' = NA))
  else {
    res_dt <- as.data.table(mvcapa_cor_res$anoms)
    res_dt <- res_dt[start != end]
    names(res_dt)[4] <- "mean_change"
    return(res_dt)
  }
}

point_anomalies <- function(mvcapa_cor_res) {
  if (nrow(mvcapa_cor_res$anoms) == 0)
    return(data.frame('location' = NA, 'variate' = NA, 'strength' = NA))
  else {
    res_dt <- as.data.table(mvcapa_cor_res$anoms)
    res_dt <- res_dt[start == end][, 2:4]
    names(res_dt)[c(1, 3)] <- c("location", "strength")
    return(res_dt)
  }
}



#### R functions #####

init_precision_mat <- function(precision_mat, tol = sqrt(.Machine$double.eps)) {
  p <- ncol(precision_mat)
  precision_mat[abs(precision_mat) <= tol] <- 0
  lower_nbs <- lower_nbs(precision_mat)
  extended_nbs <- extended_lower_nbs(lower_nbs)
  # lower_nbs <- lapply(1:p, get_neighbours_less_than, sparse_mat = precision_mat)
  # extended_nbs <- lapply(1:p, function(d) remaining_neighbours_below(d, lower_nbs, p))
  list('Q' = precision_mat, 'nbs' = lower_nbs, 'extended_nbs' = extended_nbs)
}

mvcapa_corR <- function(x, Q, b = 1,
                        min_seg_len = 3, max_seg_len = round(nrow(x) / 3),
                        prune = TRUE, print_progress = FALSE) {
  # TODO: Restrict min_seg_len and max_seg_len
  # Robustly normalise x
  n <- nrow(x)
  p <- ncol(x)
  # x <- anomaly::robustscale(x)
  # TODO: Implement standardisation based on Q.
  # Now: Assumes x and Q are standardised before input.
  # precision_mat <- estimate_precision_mat(x, adj_model)
  if (!any(class(Q) == "dsCMatrix")) Q <- Matrix::Matrix(Q, sparse = TRUE)
  Q_obj <- init_precision_mat(Q)
  # Q_order <- 1:p

  x <- rbind(rep(0, p), x)  # So the index m always means m - 1.

  # Setting up penalties
  penalty_names <- c("sparse", "dense", "point")
  penalty <- lapply(penalty_names, get_penalty, n = n, p = p, b = b)
  names(penalty) <- penalty_names

  # Initialising the DP for the optimal cost.
  C <- rep(0, n + 1)
  S <- matrix(0, nrow = n + 1, ncol = n + 1)
  J_max <- list()
  S_point <- rep(0, n + 1)
  J_point <- list()
  anom <- matrix(0, nrow = n + 1, ncol = 2)
  prune_t <- rep(0, n + 1)

  # Running OP
  for (m in (min_seg_len + 1):(n + 1)) {
    if (print_progress)
      if (m %% 10 == 0) print(paste0(round(100 * m/(n + 1), 2), '% complete.'))
    ## Collective anomalies:
    ts <- max(m - min_seg_len - max_seg_len + 1, 1):(m - min_seg_len)
    ts <- setdiff(ts, unique(prune_t))
    J_m <- list()
    for (t in ts) {
      saving <- optimise_mvnormal_saving(x[(t + 1):m, , drop = FALSE], Q,
                                         Q_obj$nbs,
                                         Q_obj$extended_nbs,
                                         penalty$dense$alpha,
                                         penalty$sparse$beta,
                                         penalty$sparse$alpha)

      # J_m[[t]] <- Q_order[saving$J_max]  # To refer to the original variates.
      J_m[[t]] <- saving$J_max
      S[m, t] <- saving$S_max
    }
    C1s <- C[ts] + S[m, ts]
    C1_max <- max(C1s)
    t_max <- ts[which.max(C1s)]
    J_max[[m]] <- J_m[[t_max]]

    ## Point anomalies:
    point_saving <- optimise_mvnormal_saving(x[m, , drop = FALSE], Q,
                                             Q_obj$nbs,
                                             Q_obj$extended_nbs,
                                             Inf,
                                             penalty$point$beta,
                                             penalty$point$alpha)
    # J_point[[m]] <- Q_order[point_saving$J_max]  # To refer to the original variates.
    J_point[[m]] <- point_saving$J_max
    S_point[m] <- point_saving$S_max
    C2 <- C[m - 1] + S_point[m]

    C0 <- C[m - 1]
    C[m] <- max(C0, C1_max, C2)

    if (prune) {
      t_to_prune <- ts[C0 >= C1s + penalty$dense$alpha]
      # print("Prune:")
      # print(t_to_prune)
      prune_t[t_to_prune] <- t_to_prune
    }

    # anomaly_type == 0: Previous cost is maximal. Points one time-step back.
    # anomaly_type == 1: Current cost for a collective anomaly is maximal. Points to the maximising time-step.
    # anomaly_type == 2: Current cost for a point anomaly is maximal. Points to the current time-step.
    anomaly_type <- which.max(c(C0, C1_max, C2)) - 1
    if (anomaly_type == 0) {
      anom[m, ] <- c(m - 1, anomaly_type)
    } else if (anomaly_type == 1) {
      anom[m, ] <- c(t_max, anomaly_type)
    } else if (anomaly_type == 2) {
      anom[m, ] <- c(m - 1, anomaly_type)
    }
  }
  # return(list('x'       = x[-1, order(Q_order)],
  return(list('x'       = x[-1, ],
              'S'       = S,
              'J'       = J_max,
              'S_point' = S_point,
              'J_point' = J_point,
              'C'       = C,
              'Q'       = Q,
              'anom'    = anom))
}

collective_anomaliesR <- function(mvcapa_cor_res) {
  anom_list <- list()
  m <- nrow(mvcapa_cor_res$anom)
  # Indices m correspond to m - 1 in the data.
  while (m >= 1) {
    if (mvcapa_cor_res$anom[m, 2] == 0) m <- mvcapa_cor_res$anom[m, 1]
    else if (mvcapa_cor_res$anom[m, 2] == 1) {
      end <- m - 1
      start <- mvcapa_cor_res$anom[m, 1]
      J <- mvcapa_cor_res$J[[m]]
      means <- colMeans(mvcapa_cor_res$x[start:end, J, drop = FALSE])
      anom_df <- data.frame('start'       = rep(start, length(J)),
                            'end'         = rep(end, length(J)),
                            'variate'     = J,
                            'mean_change' = means)
      anom_list[[length(anom_list) + 1]] <- anom_df
      m <- start
    } else if (mvcapa_cor_res$anom[m, 2] == 2) {
      m <- mvcapa_cor_res$anom[m, 1]
    }
  }
  if (length(anom_list) == 0)
    return(data.frame('start' = NA, 'end' = NA, 'variate' = NA, 'mean_change' = NA))
  else {
    anom_df <- do.call('rbind', anom_list)
    anom_dt <- data.table::as.data.table(anom_df)
    anom_dt <- anom_dt[order(start)][, .SD[order(variate)], by = start]
    return(as.data.frame(anom_dt))
  }
}

point_anomaliesR <- function(mvcapa_cor_res) {
  anom_list <- list()
  m <- nrow(mvcapa_cor_res$anom)
  while (m >= 1) {
    # if (mvcapa_cor_res$anom[m, 2] == 0) m <- mvcapa_cor_res$anom[m, 1]
    # else if (mvcapa_cor_res$anom[m, 2] == 1) {
    #   m <- mvcapa_cor_res$anom[m, 1]
    # } else if (mvcapa_cor_res$anom[m, 2] == 2) {
    if (mvcapa_cor_res$anom[m, 2] == 2) {
      J_point <- mvcapa_cor_res$J_point[[m]]
      anom_df <- data.frame('location' = rep(m - 1, length(J_point)),
                            'variate'  = J_point,
                            'strength' = mvcapa_cor_res$x[m - 1, J_point])
      anom_list[[length(anom_list) + 1]] <- anom_df
      # m <- mvcapa_cor_res$anom[m, 1]
    }
    m <- mvcapa_cor_res$anom[m, 1]
  }
  if (length(anom_list) == 0)
      return(data.frame('location' = NA, 'variate' = NA, 'strength' = NA))
  else {
    anom_df <- do.call('rbind', anom_list)
    anom_dt <- data.table::as.data.table(anom_df)
    anom_dt <- anom_dt[order(location)][, .SD[order(variate)], by = location]
    return(as.data.frame(anom_dt))
  }
}
