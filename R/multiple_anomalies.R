BIC_precision_mat <- function(A, S, n, rho) {
  nnz <- sum(A[lower.tri(A, diag = TRUE)] != 0)
  - log(det(A)) + sum(diag(A %*% S)) + log(n) / n * nnz
}

select_precision_mat <- function(glasso_res, S, n) {
  BIC_per_lambda <- unlist(lapply(1:length(glasso_res$lambda), function(i) {
    BIC_precision_mat(glasso_res$icov[[i]], S, n, glasso_res$lambda[i])
  }))
  glasso_res$icov[[which.min(BIC_per_lambda)]]
}

estimate_sparse_precision_mat <- function(x, lambda_min_ratio = 0.1, nlambda = 10,
                                          robust = TRUE, standardise = TRUE) {
  if (standardise) x <- anomaly::robustscale(x)
  if (robust) S <- robust_cov_mat(x)
  else        S <- var(x)
  glasso_res <- huge::huge(S, lambda.min.ratio = lambda_min_ratio,
                           nlambda = nlambda, method = 'glasso')
  select_precision_mat(glasso_res, S, nrow(x))
}

robust_cov_mat <- function(x) {
  p <- ncol(x)
  n <- nrow(x)
  scores_x <- apply(x, 2, function(x_i) qnorm(rank(x_i)/(n + 1)))
  mad_x <- diag(apply(x, 2, mad))
  mad_x %*% cor(scores_x) %*% mad_x
}


mvcapa_cor <- function(x, precision_mat = 'automatic', lambda_min_ratio = 0.1, b = 1,
                       min_seg_len = 1, max_seg_len = nrow(x) / 3) {
  # Robustly normalise x
  n <- nrow(x)
  p <- ncol(x)
  x <- anomaly::robustscale(x)
  if (precision_mat == 'automatic')
    precision_mat <- estimate_sparse_precision_mat(x, lambda_min_ratio, standardise = FALSE)
  precision_mat_order <- 1:p
  band_precision_mat <- band(precision_mat)
  reordered_precision_mat <- cuthill_mckee(precision_mat, return_all = TRUE)
  if (reordered_precision_mat$band < band_precision_mat) {
    # print('Reordering')
    precision_mat <- reordered_precision_mat$x
    precision_mat_order <- reordered_precision_mat$order
    band_precision_mat <- reordered_precision_mat$band
    x <- x[, precision_mat_order]
  }
  print('Estimated precision matrix:')
  print(precision_mat)
  print(paste0('Band = ', band(precision_mat)))
  if (band_precision_mat > 15)
    warning(paste0('Computation will be slow because the band of the precision matrix is ', band_precision_mat))
  if (band_precision_mat >= 20)
    stop(paste0('The band of the precision matrix is too high for computation: ', band_precision_mat))
  precision_mat_obj <- init_precision_mat(precision_mat)

  x <- rbind(rep(0, p), x)  # So the index m is always means m - 1.

  # Initialising the DP for the optimal cost.
  C <- rep(0, n + 1)
  S <- matrix(0, nrow = n + 1, ncol = n + 1)
  J <- list()
  anom <- matrix(0, nrow = n + 1, ncol = 2)

  # Running OP
  for (m in (min_seg_len + 1):(n + 1)) {
    if (m %% 10 == 0) print(paste0(round(100 * m/(n + 1), 2), '% complete.'))
    # print(m)
    J[[m]] <- list()
    # TODO: Restrict min_seg_len and max_seg_len
    max_ind <- max(m - min_seg_len - max_seg_len + 1, 1):(m - min_seg_len)
    for (t in max_ind) {
      optim_savings_obj <- optim_penalised_savings_c(x[(t + 1):m, , drop = FALSE],
                                                     n, precision_mat_obj)
      # optim_savings_obj$J_max
      # J[[m]][[t]] <- optim_savings_obj$J_max
      J[[m]][[t]] <- precision_mat_order[optim_savings_obj$J_max]  # To refer to the original variates.
      S[m, t] <- optim_savings_obj$B_max
    }
    C0 <- C[m - 1]
    C1s <- C[max_ind] + S[m, max_ind]
    C1_max <- max(C1s)
    t_max <- max_ind[which.max(C1s)]
    C[m] <- max(C0, C1_max)

    # anomaly_type == 0: Previous cost is maximal. Points one time-step back.
    # anomaly_type == 1: Current cost for a collective anomaly is maximal. Points to the maximising time-step.
    # anomaly_type == 2: Current cost for a point anomaly is maximal. Points to the current time-step.
    anomaly_type <- which.max(c(C0, C1_max)) - 1
    if (anomaly_type == 0) {
      anom[m, ] <- c(m - 1, anomaly_type)
    } else if (anomaly_type == 1) {
      anom[m, ] <- c(t_max + 1, anomaly_type)
    }
  }
  return(list('x' = x[, order(precision_mat_order)],
              'S' = S,
              'J' = J,
              'C' = C,
              'A' = precision_mat,
              'anom' = anom))
}

init_precision_mat <- function(precision_mat) {
  p <- ncol(precision_mat)
  lower_nbs <- lapply(1:p, get_neighbours_less_than, sparse_mat = precision_mat)
  extended_nbs <- lapply(1:p, function(d) remaining_neighbours_below(d, lower_nbs, p))
  list('A' = precision_mat, 'extended_nbs' = extended_nbs)
}

init_data_setup <- function(n = 200, p = 4, proportions = sqrt(p)/p, mu = 1,
                            locations = n - durations - 1, durations = 10,
                            change_type = 'adjacent', cor_mat_type = 'lattice',
                            rho = 0.8, band = 2, min_nbs = 1, max_nbs = 3) {
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
    }
  }

  Sigma_obj <- get_Sigma(cor_mat_type)
  list('n'                 = n,
       'p'                 = p,
       'mu'                = mu,
       'Sigma'             = Sigma_obj$mat,
       'Sigma_inv'         = Sigma_obj$inverse,
       'locations'         = locations,
       'durations'         = durations,
       'proportions'       = proportions,
       'change_type'       = change_type,
       'precision_mat_obj' = init_precision_mat(Sigma_obj$inverse))
}

simulate_mvcapa_cor <- function(setup = init_data_setup(), b = 1,
                                min_seg_len = 5, max_seg_len = round(nrow(sim_data) / 3)) {
  sim_data <- simulate_cor(n           = setup$n,
                           p           = setup$p,
                           mu          = setup$mu,
                           Sigma       = setup$Sigma,
                           locations   = setup$locations,
                           durations   = setup$durations,
                           proportions = setup$proportions,
                           change_type = setup$change_type)
  print(setup$Sigma_inv)
  res <- mvcapa_cor(sim_data,
                    # setup$precision_mat_obj,
                    b           = b,
                    min_seg_len = min_seg_len,
                    max_seg_len = max_seg_len)
  res
}

collective_anomalies <- function(mvcapa_cor_res) {

  anom_list <- list()
  m <- nrow(mvcapa_cor_res$anom)
  while (m >= 1) {
    if (mvcapa_cor_res$anom[m, 2] == 0) m <- mvcapa_cor_res$anom[m, 1]
    else if (mvcapa_cor_res$anom[m, 2] == 1) {
      end <- m
      start <- mvcapa_cor_res$anom[m, 1]
      J <- mvcapa_cor_res$J[[end]][[start - 1]]
      means <- colMeans(mvcapa_cor_res$x[start:end, J, drop = FALSE])
      m <- start - 1
      anom_df <- data.frame('start'       = rep(start, length(J)),
                            'end'         = rep(end, length(J)),
                            'variate'     = J,
                            'mean_change' = means)
      anom_list[[length(anom_list) + 1]] <- anom_df
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

capa_line_plot <- function(object, epoch = dim(object$x)[1],
                           subset = 1:ncol(object$x), variate_names = TRUE) {
    # creating null entries for ggplot global variables so as to pass CRAN checks
    x <- value <- ymin <- ymax <- x1 <- x2 <- y1 <- y2 <- x1 <- x2 <- y1 <- y2 <- NULL
    data_df <- as.data.frame(object$x)
    names <- paste("y", 1:ncol(object$x), sep = "")
    colnames(data_df) <- names
    data_df <- as.data.frame(data_df[, subset, drop = FALSE])
    n <- nrow(data_df)
    p <- ncol(data_df)
    data_df <- cbind(data.frame("x" = 1:n), data_df)
    data_df <- reshape2::melt(data_df, id = "x")
    out <- ggplot2::ggplot(data = data_df)
    out <- out + ggplot2::aes(x = x, y = value)
    out <- out + ggplot2::geom_point()
    # highlight the collective anomalies
    c_anoms <- collective_anomalies(object)
    c_anoms <- c_anoms[c_anoms$variate %in% subset, ]
    if(!any(is.na(c_anoms)) & nrow(c_anoms) > 0)
    {
        c_anoms_data_df <- c_anoms[, 1:3]
        names(c_anoms_data_df) <- c(names(c_anoms_data_df)[1:2], "variable")
        c_anoms_data_df$variable <- names[c_anoms_data_df$variable]
        c_anoms_data_df$ymin <- -Inf
        c_anoms_data_df$ymax <- Inf
        out <- out + ggplot2::geom_rect(data = c_anoms_data_df,
                                       inherit.aes = FALSE,
                                       mapping = ggplot2::aes(xmin = start,
                                                              xmax = end,
                                                              ymin = ymin,
                                                              ymax = ymax),
                                       fill = "blue", alpha=0.4)
        # c_anoms_data_df$start <- c_anoms_data_df$start + c_anoms$start.lag
        # c_anoms_data_df$end<-c_anoms_data_df$end-c_anoms$end.lag
        # out <- out + ggplot2::geom_rect(data = c_anoms_data_df,
        #                                 inherit.aes = FALSE,
        #                                 mapping = ggplot2::aes(xmin = start,
        #                                                        xmax = end,
        #                                                        ymin = ymin,
        #                                                        ymax = ymax),
        #                                 fill = "blue", alpha = 0.5)
    }
    # out<-out+facet_grid(variable~.,scales="free_y")
    # highlight the point anomalies
    # TODO: Add point anomalies.
    # p_anoms<-point_anomalies(object)
    # p_anoms<-p_anoms[p_anoms$variate %in% subset,]
    # if(!any(is.na(p_anoms)) & nrow(p_anoms) > 0)
    #     {
    #         p_anoms_data_df<-Reduce(rbind,Map(function(a,b) data_df[data_df$variable==names[a] & data_df$x==b,],p_anoms$variate,p_anoms$location))
    #         out<-out+geom_point(data=p_anoms_data_df,colour="red", size=1.5)
    # }

    out <- out + ggplot2::facet_grid(factor(variable, levels = rev(names)) ~ .,
                                     scales = "free_y")
    # grey out the data after epoch
    if(epoch != nrow(object$x))
    {
	      d <- data.frame(variable = names[subset], x1 = epoch, x2 = n,
	                    y1 = -Inf, y2 = Inf)
        out <- out + geom_rect(data = d, inherit.aes = FALSE,
                               mapping = aes(xmin = x1, xmax = x2,
                                             ymin = y1, ymax = y2),
                               fill = "yellow", alpha = 0.2)
    }
    if(variate_names == FALSE)
    {
        out <- out + theme(strip.text.y = element_blank())
    }
    # change background
    out <- out + ggplot2::theme(
                                # Hide panel borders and remove grid lines
                                panel.border = ggplot2::element_blank(),
                                panel.grid.major = ggplot2::element_blank(),
                                panel.grid.minor = ggplot2::element_blank(),
                                # Change axis line
                                axis.line = ggplot2::element_line(colour = "black")
                              )
    return(out)
}

capa_tile_plot <- function(object, variate_names = FALSE,
                           epoch = dim(object$x)[1], subset = 1:ncol(object$x)) {
    # nulling out variables used in ggplot to get the package past CRAN checks
    x1 <- y1 <- x2 <- y2 <- variable <- value <- NULL
    df <- as.data.frame(object$x)
    df <- as.data.frame(df[, subset, drop = FALSE])
    # normalise data
    for(i in 1:ncol(df))
    {
        df[, i] <- (df[, i] - min(df[, i])) / (max(df[, i]) - min(df[, i]))
    }
    n <- data.frame("n" = seq(1, nrow(df)))
    molten.data <- reshape2::melt(cbind(n, df), id = "n")
    out <- ggplot2::ggplot(molten.data, ggplot2::aes(n, variable))
    out <- out + ggplot2::geom_tile(ggplot2::aes(fill = value))
    # get any collective anomalies
    c_anoms <- collective_anomalies(object)
    c_anoms <- c_anoms[c_anoms$variate %in% subset, ]
    c_anoms <- unique(c_anoms[, 1:2])
    if(!any(is.na(c_anoms)) & nrow(c_anoms) > 0)
    {
        ymin <- 0
        ymax <- ncol(df)
        out <- out + ggplot2::annotate("rect", xmin = c_anoms$start,
                                       xmax = c_anoms$end, ymin = ymin,
                                       ymax = ymax + 1, alpha = 0.0,
                                       color = "red", fill = "yellow")
    }
    if(epoch != nrow(object$x))
    {
        d <- data.frame(x1 = epoch, x2 = nrow(object$x), y1 = -Inf, y2 = Inf)
        out <- out + ggplot2:::geom_rect(data = d, inherit.aes = FALSE,
                                         mapping = ggplot2::aes(xmin = x1,
                                                                xmax = x2,
                                                                ymin = y1,
                                                                ymax = y2),
                                         fill = "yellow", alpha=0.2)
    }
    if(variate_names == FALSE)
    {
        out <- out + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                    axis.title = ggplot2::element_blank())
    }
    out <- out + ggplot2::theme(
                                # Hide panel borders and remove grid lines
                                panel.border = ggplot2::element_blank(),
                                panel.grid.major = ggplot2::element_blank(),
                                panel.grid.minor = ggplot2::element_blank(),
                                # Change axis line
                                axis.line = ggplot2::element_line(colour = "black")
                                )
    return(out)
}

plot_capa <- function(object, epoch = dim(object$x)[1],
                      subset = 1:ncol(object$x), variate_names = TRUE) {
  if (length(subset) <= 20)
    return(capa_line_plot(object, epoch, subset, variate_names))
  else
    return(capa_tile_plot(object, variate_names, epoch, subset))
}








