
#' @export
sim_cpt_distr <- function(out_file,
                          data = init_data(n = 200, p = 10, rho = 0.8, band = 2,
                                           locations = 80, durations = 120),
                          method = method_params(cost = "mvlrt"),
                          n_sim = 100, seed = NA) {
  add_setup_info <- function(res) {
    res <- cbind(res,
                 as.data.table(data[!(grepl("Sigma", names(data)) |
                                        names(data) == "changing_vars")]),
                 as.data.table(method[names(data) != "maxsl"]))
    res$n_sim <- n_sim
    res$seed <-  seed
    res
  }

  message(paste0("Estimating changepoint distribution for n=", data$n,
                 ", p=", data$p,
                 ", cost=", method$cost,
                 ", precision=", data$precision_type,
                 ", band=", data$band,
                 ", rho=", data$rho,
                 ", prop=", data$proportions,
                 ", vartheta=", data$vartheta,
                 ", shape=", data$shape,
                 ", location=", data$locations,
                 ", precision_est_struct=", method$precision_est_struct,
                 ", est_band=", method$est_band,
                 "."))
  all_params <- c(data, method, list("n_sim" = n_sim))
  if (already_estimated(out_file, all_params, read_cpt_distr)) return(NULL)

  if (!is.na(seed)) set.seed(seed)
  res <- data.table(cpt = replicate(n_sim, simulate_detection(data, method, standardise_output = FALSE)$cpt))
  fwrite(add_setup_info(res), paste0("./results/", out_file), append = TRUE)
}

#' #' @export
many_cpt_distr <- function(out_file = "cpt_distr.csv",
                           variables = list("cost" = c("inspect", "inspect_iid",
                                                       "mvlrt", "mvlrt_iid"),
                                            "rho" = c(0.7, 0.9, 0.99),
                                            "proportions" = c(0.1, 0.3, 1)),
                           data   = init_data(n = 100, locations = 50, durations = 50),
                           method = method_params(),
                           loc_tol = 10) {
  params_list <- split_params(
    expand_list(c(data, method), variables),
    list("data"   = names(data),
         "method" = names(method))
  )
  Map(sim_cpt_distr,
      data = params_list$data,
      params = params$list$method,
      seed = get_sim_seeds(params_list, variables),
      MoreArgs = list("out_file" = out_file,
                      "n_sim"    = n_sim))
  NULL
}

#' @export
cpt_runs <- function(n_sim = 100) {
  out_file <- "cpt_distr.csv"

  #### FIX DURATIONS
  # many_cpt_distr(out_file,
  #                init_data(n = 100, p = 10),
  #                method_params(b = 1),
  #                n_sim = n_sim,
  #                props = c(0.1, 1),
  #                varthetas = c(0.5, 1, 2),
  #                rhos = c(0.7, 0.9, 0.99),
  #                locations = c(10, 50),
  #                precision_est_structs = c("correct", "banded"),
  #                est_bands = 0)
  # many_cpt_distr(out_file,
  #                init_data(n = 100, p = 10),
  #                method_params(b = 0.3),
  #                n_sim = n_sim,
  #                props = c(0.1, 1),
  #                varthetas = c(0.5, 1, 2),
  #                rhos = c(0.7, 0.9, 0.99),
  #                locations = c(10, 50),
  #                precision_est_structs = c("correct", "banded"),
  #                est_bands = 0)
  # many_cpt_distr(out_file,
  #                init_data(n = 200, p = 100),
  #                method_params(b = 1),
  #                n_sim = n_sim,
  #                props = c(0.01, 1),
  #                varthetas = c(0.5, 1, 2),
  #                rhos = c(0.7, 0.9, 0.99),
  #                locations = c(10, 50, 100),
  #                precision_est_structs = c("correct", "banded"),
  #                est_bands = 0)
}

#' @export
plot_cpt_distr <- function(file_name,
                           data = init_data(n = 100, p = 10, rho = 0.9, band = 2,
                                            proportions = 0.1, vartheta = 1,
                                            locations = 50, durations = 50),
                           params = method_params(),
                           n_sim = 200, vars_in_title = NA) {
  all_params <- c(data, params, list("n_sim" = n_sim))
  res <- do.call("rbind",
                 Map(read_cpt_distr,
                     expand_list(all_params, list("cost" = c("iid", "cor"))),
                     MoreArgs = list(file_name = file_name)))
  title <- make_title(all_params, cpt_distr_title_parts(vars_in_title), "cpt")
  ggplot2::ggplot(res, ggplot2::aes(cpt, colour = cost)) +
    ggplot2::stat_bin(alpha = 0.8, ggplot2::aes(y=..density..), geom = "step",
                      position = "identity", binwidth = 1) +
    ggplot2::ggtitle(title) +
    ggplot2::scale_x_continuous(name = "Estimated cpt",
                                limits = c(0, res$n[1]),
                                breaks = res$location[1] - 1)
}

#' @export
boxplot_cpt_distr <- function(file_name,
                              data = init_data(n = 100, p = 10, rho = 0.9, band = 2,
                                               vartheta = 1, proportions = 0.1,
                                               locations = 50, durations = 50),
                              params = method_params(),
                              n_sim = 200, vars_in_title = NA) {
  data_list <- expand_list(data, list("vartheta" = c(0.5, 1, 2)))
  res <- do.call("rbind", lapply(data_list, function(data) {
    read_single_cpt_distr(file_name, data, params, n_sim)
  }))
  title <- make_title(c(data, params), cpt_distr_title_parts(vars_in_title), "cpt")
  ggplot2::ggplot(res, ggplot2::aes(as.factor(vartheta), cpt, fill = cost)) +
    ggplot2::geom_boxplot() +
    ggplot2::scale_x_discrete("Signal strength") +
    ggplot2::scale_y_continuous("Cpt estimate") +
    ggplot2::ggtitle(title)
}


#' @export
grid_plot_cpt_distr <- function(plot_func,
                                variables = list("rho" = c(0.7, 0.9, 0.99),
                                                 "proportions" = c(0.1, 1)),
                                data = init_data(n = 100, p = 10, band = 2,
                                                 vartheta = 1, locations = 50,
                                                 durations = 50),
                                params = method_params(),
                                n_sim = 200) {
  if (length(variables) > 2)
    stop("Maximum two variables at a time.")
  params_list <- split_lists(expand_list(c(data, params), variables),
                             list("data"   = names(data),
                                  "method" = names(params)))
  plots <- Map(plot_func, data = params_list$data,
               params = params_list$method,
               MoreArgs = list("file_name"     = "cpt_distr.csv",
                               "n_sim"         = n_sim,
                               "vars_in_title" = names(variables)))
  vars_in_title <- cpt_distr_title_parts(NA)
  vars_in_title <- vars_in_title[!vars_in_title %in% names(variables)]
  dims <- c(length(variables[[2]]), length(variables[[1]]))
  grid_plot(plots, dims, make_title(c(data, params), vars_in_title, "cpt"))
}
