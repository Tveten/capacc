#' @export
tune_penalty <- function(data = init_data(mu = 0), method = method_params(),
                         tuning = tuning_params(), known = FALSE, seed = NA) {

  add_setup_info <- function(res) {
    which_data <- !(grepl("Sigma", names(data)) | names(data) == c("changing_vars"))
    if (known) excluded_method_names <- "b"
    else excluded_method_names <- c("b", "size_mu")
    res <- cbind(res,
                 as.data.table(data[which_data]),
                 as.data.table(method[!names(method) %in% excluded_method_names]),
                 as.data.table(tuning[names(tuning) != "init_b"]))
    res$seed <-  seed
    res[1]  # There might be more than one anomaly => > 1 location => nrow(res) > 1.
  }

  fp_rate <- function(b) {
    fps <- unlist(lapply(1:tuning$tuning_n_sim, function(i) {
      method$b <- b
      if (known) {
        return(simulate_detection_known(data, method)$S_max > 0)
      } else {
        sim <- simulate_detection(data, method, standardise_output = TRUE)
        if (is.null(sim$collective)) return(nrow(sim$cpt) > 0)
        else return(nrow(sim$collective) > 0)
      }
    }))
    data.table("b" = b, "fp" = mean(fps), "diff" = mean(fps) - tuning$alpha)
  }

  new_b <- function(res) {
    res <- res[order(diff)]
    b_lower <- as.numeric(res[diff <= 0, min(b)])
    b_upper <- as.numeric(res[diff > 0][diff == min(diff), max(b)])
    if (is.na(b_lower) || is.na(b_upper) || is.infinite(b_lower) || is.infinite(b_upper))
      stop("Will not converge. init_b does not cover the selected alpha level.")
    exp((log(b_lower) + log(b_upper)) / 2)
  }

  if (known) {
    file_name <- "penalties_known.csv"
    read_func <- read_penalties_known
  } else {
    file_name <- "penalties.csv"
    read_func <- read_penalties
  }

  seed <- get_previous_seed(seed, c(data, method, tuning), read_func, file_name)
  if (!is.na(seed)) set.seed(seed)
  data$mu <- 0
  data$vartheta <- 0
  data$point_mu <- data$point_locations <- data$point_proportions <- NA
  message(paste0("Tuning ", method$cost,
                 " penalty for n=", data$n,
                 ", p=", data$p,
                 ", loc=", data$locations[1],
                 ", dur=", data$durations[1],
                 ", precision_type=", data$precision_type,
                 ", block_size=", data$block_size[1],
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", precision_est_struct=", method$precision_est_struct,
                 ", est_band=", method$est_band))
  res <- do.call('rbind', Map(fp_rate, tuning$init_b))
  print(res)
  while (!any(abs(res$diff) <= tuning$alpha_tol + 1e-6) && nrow(res) <= tuning$tuning_max_iter) {
    res <- rbind(res, fp_rate(new_b(res)))
    print(res[.N])
  }
  res <- add_setup_info(res[which.min(abs(diff))])
  fwrite(res, paste0("./results/", file_name), append = TRUE)
}

#' @export
get_tuned_penalty <- function(data = init_data(mu = 0), method = method_params(),
                              tuning = tuning_params(), known = FALSE, seed = NA) {
  if (method$cost == "cor_exact") method$cost <- "cor"
  else if (method$cost %in% c("gflars", "var_pgl")) return(list(b = NA))
  if (known) {
    file_name <- "penalties_known.csv"
    read_func <- read_penalties_known
  } else {
    file_name <- "penalties.csv"
    read_func <- read_penalties
  }
  if (!known && data$n > 200) {
    data$n <- 200
    if (data$p <= 100) data$n <- 200
    else if (data$p > 100) data$n <- 2 * data$p
    data$locations <- data$durations <- 2
    data$change_type <- "adjacent"
    data$proportions <- data$vartheta <- data$mu <- 0
    data$point_mu <- data$point_locations <- data$point_proportions <- NA
  }
  if (!file.exists(paste0("./results/", file_name))) {
    message(paste0("File does not exist. Making file ", file_name, " in ./results/"))
    tune_penalty(data, method, tuning, known, seed)
    return(get_tuned_penalty(data, method, tuning, known, seed))
  } else {
    all_pen <- read_results(file_name)
    res <- read_func(all_pen, c(data, method, tuning))
    if (nrow(res) == 0) {
      tune_penalty(data, method, tuning, known, seed)
      return(get_tuned_penalty(data, method, tuning, known, seed))
    } else
      return(res[1])
  }
}

#' @export
tune_many <- function(n = 100, p = 10, band = 2) {
  init_b <- c(0.5, 1, 4)
  precision_types <- c("banded", "block_banded", "lattice")
  costs <- c("iid", "cor")
  rhos = c(0.5, 0.7, 0.8, 0.9, 0.95, 0.99)
  lapply(precision_types, function(precision_type) {
    lapply(rhos, function(rho) {
      lapply(costs, function(cost) {
        if (precision_type == "block_banded") m <- p - 1
        else m <- p
        data <- init_data(n = n, p = p, precision_type = precision_type,
                          band = band, rho = rho, block_size = m)
        get_tuned_penalty(data, method_params(cost), tuning_params(init_b = init_b))
      })
    })
  })
}

#' @export
plot_penalty <- function(data, method, tuning, known = FALSE, vars_in_title = NA) {
  all_params <- c(data, method, tuning)
  params_list <- combine_lists(split_params(
    expand_list(all_params,
                list("cost" = c("iid", "cor"),
                     "precision_est_struct" = c(NA, "correct"),
                     "rho" = c(0.01, 0.2, 0.5, 0.7, 0.9, 0.99))),
    list("data"    = names(data),
         "method"  = names(method),
         "tuning"  = names(tuning))
  ))
  if (known) {
    all_pen <- read_results("penalties_known.csv")
    read_func <- read_penalties_known
  } else {
    all_pen <- read_results("penalties.csv")
    read_func <- read_penalties
  }
  res <- do.call("rbind",
                 Map(read_func,
                     params_list,
                     MoreArgs = list(res = all_pen)))
  res <- rename_precision_est_struct(res)
  print(res[cost == "cor", .(rho, b, cost, precision_est_struct)])
  print(res[cost == "iid", .(rho, b, cost, precision_est_struct)])
  title <- make_title(all_params, penalties_title_parts(vars_in_title))
  ggplot2::ggplot(data = res, ggplot2::aes(x = rho, y = b,
                                           colour = cost,
                                           linetype = precision_est_struct)) +
    ggplot2::geom_line() +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous(latex2exp::TeX("$\\rho$")) +
    ggplot2::scale_y_continuous("Penalty scale") +
    ggplot2::scale_colour_discrete(name = "Cost") +
    ggplot2::scale_linetype_discrete(name = "Precision estimation")
}

penalty_plots <- function() {
  banded_data <- init_data(n = 100, p = 8, precision_type = "banded",
                           band = 2, locations = 50, durations = 10,
                           change_type = "adjacent", shape = 0)
  plot_penalty(banded_data, method_params(), tuning_params(), known = TRUE)
}
