###############################################################
## load libraries
###############################################################
library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggvis)
library(corrplot)

# glimpse(data)

###############################################
## create new variables
###############################################
#' @export
load_pump_data <- function() {
  ###############################################################
  ## read data
  ###############################################################
  data <- readRDS(file = "./data/SubseaPumpData/data.Rds")

  data <- data %>%
  dplyr::mutate(differential_pressure = discharge_pressure - suction_pressure,
                mpfm_gas_flow_rate = mpfm_gvf / 100 * mpfm_flow_rate,
                mpfm_oil_flow_rate =  0.08 * (1 - mpfm_gvf / 100) * mpfm_flow_rate,
                mpfm_water_flow_rate = 0.92 * (1 - mpfm_gvf / 100) * mpfm_flow_rate,
                running = as_factor(if_else(speed >= 1200, "Running", "Not Running")),
                year = lubridate::year(time),
                month = lubridate::month(time),
                year_month = as.numeric(paste0(lubridate::year(time),
                                               str_pad(lubridate::month(time), 2,
                                                       pad = "0"))),
                date = lubridate::date(time))

  ###############################################
  ## define variable for healthy vs. non-helathy data
  ###############################################
  ## create intervals
  a <- pump_anom_intervals()
  ## define variable
  data %>%
    dplyr::mutate(status = as_factor(if_else(date %within% a, "Not Healthy", "Healthy")))
}

pump_anoms <- function() {
  s1 <- ymd("1970-01-11", tz = "CET")
  e1 <- ymd("1970-04-12", tz = "CET")

  s2 <- ymd("1970-06-05", tz = "CET")
  e2 <- ymd("1970-06-14", tz = "CET")

  s3 <- ymd("1972-12-14", tz = "CET")
  e3 <- ymd("1973-07-27", tz = "CET")

  s4 <- ymd("1978-03-29", tz = "CET")
  e4 <- ymd("1979-01-14", tz = "CET")
  list(c(s1, e1), c(s2, e2), c(s3, e3), c(s4, e4))
}

pump_anom_intervals <- function() {
  lapply(pump_anoms(), function(x) lubridate::interval(x[1], x[2]))
}

failure_period <- function(nr) {
  a1 <- lubridate::interval(start = ymd("1970-01-01", tz = "CET"),
                            end = ymd("1970-04-13", tz = "CET"))
  a2 <- lubridate::interval(start = ymd("1970-05-12", tz = "CET"),
                            end = ymd("1970-06-15", tz = "CET"))
  a3 <- lubridate::interval(start = ymd("1970-10-15", tz = "CET"),
                            end = ymd("1973-07-28", tz = "CET"))
  a4 <- lubridate::interval(start = ymd("1973-11-09", tz = "CET"),
                            end = ymd("1979-01-15", tz = "CET"))
  a <- list(a1, a2, a3, a4)
  a[nr]
}

###############################################
## explorative statistics
###############################################
## data availability
## number of observations each week
# data %>%
#   dplyr::mutate(week = floor_date(x = time, unit = "weeks"),
#                 week = as_date(week)) %>%
#   dplyr::group_by(week) %>%
#   dplyr::summarise(n = n()) %>%
#   ggvis(x = ~week, y = ~n) %>%
#   layer_points() %>%
#   add_axis("x", title = "data") %>%
#   add_axis("y", title = "count") %>%
#   add_axis("x", orient = "top", ticks = 0, title = "Data availabily",
#            properties = axis_props(axis = list(stroke = "white")))
#
# ## display daily subsea barrier temp for periods running/not running
# data %>%
#   dplyr::select(date, subsea_barrier_temperature, running) %>%
#   dplyr::group_by(date, running) %>%
#   dplyr::summarise(mean_temp = mean(subsea_barrier_temperature)) %>%
#   ggvis(x = ~date, y = ~mean_temp, stroke = ~running) %>%
#   layer_lines() %>%
#   add_relative_scales() %>%
#   add_axis("x", title = "date") %>%
#   add_axis("x", orient = "top", ticks = 0, title = "Daily mean subsea barrier fluid temperature",
#            properties = axis_props(axis = list(stroke = "white"))) %>%
#   add_legend(scales = "stroke",
#              title = "",
#              properties = legend_props(
#                legend = list(
#                  x = scaled_value("x_rel", 0.45),
#                  y = scaled_value("y_rel", -0.15))))

###############################################
## descriptive statistics for running/status combinations
###############################################
# data %>%
#   dplyr::group_by(running, status) %>%
#   dplyr::summarise(mean(subsea_barrier_temperature))

# Variables: 22
# $ time
# $ mpfm_gvf
# $ discharge_pressure
# $ discharge_temperature
# $ mpfm_flow_rate
# $ mpfm_pressure
# $ mpfm_temperature
# $ output_current
# $ output_power
# $ solenoid_dump_valve_count
# $ solenoid_feed_valve_count
# $ speed
# $ subsea_barrier_pressure
# $ subsea_barrier_temperature
# $ suction_pressure
# $ tank_level
# $ differential_pressure
# $ running
# $ year
# $ month
# $ date
# $ status
## display daily subsea barrier temp for periods running/not running, healthy/not healty.
relevant_vars <- function(type = "all") {
  temperatures <- c("discharge_temperature",
                    "mpfm_temperature",
                    "subsea_barrier_temperature")
  pressures <- c("discharge_pressure",
                 "mpfm_pressure",
                 "subsea_barrier_pressure",
                 "suction_pressure",
                 "differential_pressure")
  other <- c("mpfm_gvf",
             "mpfm_gas_flow_rate",
             "mpfm_oil_flow_rate",
             "mpfm_flow_rate",
             "output_current",
             "output_power",
             "speed")
  if (type == "temperatures") return(temperatures)
  else if (type == "pressures") return(pressures)
  else return(c(temperatures, pressures, other))
}

daily_means <- function(data, vars = "all") {
  data <- data %>%
    dplyr::filter(running == "Running") %>%
    dplyr::select(c("date", "status", relevant_vars(vars))) %>%
    dplyr::group_by(date, status) %>%
    dplyr::summarise_all(mean) %>%
    dplyr::mutate(running = "Running")
  data$ind <- 1:nrow(data)
  data
}

anonymise <- function(p) {
    p + ggplot2::theme(axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       axis.ticks=element_blank(),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank()) +
    scale_colour_discrete(name = "", labels = c("A", "B"))
}

plot_daily_mean <- function(data, var, plot_type = "point",
                            colour_by = "status", anonymise = FALSE) {
  dat <- data %>%
    dplyr::filter(running == "Running") %>%
    dplyr::select(date, starts_with(var), status) %>%
    dplyr::group_by(date, status) %>%
    dplyr::summarise(mean_temp = mean(eval(parse(text = var))))
  dat$ind <- 1:nrow(dat)
  dat$nostatus <- "nostatus"

  if (plot_type == "point")
  p <- ggplot2::ggplot(data = dat, ggplot2::aes(x = ind, y = mean_temp,
                                                colour = get(colour_by))) +
      ggplot2::geom_point(size = 0.7)
  else if (plot_type == "hist")
  p <- ggplot2::ggplot(data = dat, ggplot2::aes(x = mean_temp, fill = get(colour_by))) +
      ggplot2::geom_histogram(ggplot2::aes(y = ..density..),
                              bins = 50, position = "dodge")
  else if (plot_type == "density")
  p <- ggplot2::ggplot(data = dat, ggplot2::aes(x = mean_temp, colour = get(colour_by))) +
      ggplot2::geom_density()
  if (anonymise) p <- anonymise(p)
  else p <- p + ggplot2::ggtitle(paste0("Daily mean ", var))
  p
}

#' @export
group_plot_daily_means <- function(data, var_type = "all", plot_type = "point",
                                   colour_by = "status", anonymise = FALSE) {
  vars <- relevant_vars(var_type)
  plots <- Map(plot_daily_mean, var = vars,
               MoreArgs = list(data = data, plot_type = plot_type,
                               colour_by = colour_by, anonymise = anonymise))
  ggpubr::ggarrange(plotlist = plots, nrow = length(vars), ncol = 1,
                    common.legend = TRUE, legend = "right")
}

plot_monthly_mean <- function(data, var) {
  data %>%
    dplyr::filter(running == "Running") %>%
    dplyr::select(year_month, starts_with(var), status) %>%
    dplyr::group_by(year_month, status) %>%
    dplyr::summarise(mean_temp = mean(eval(parse(text = var)))) %>%
    ggvis(x = ~year_month, y = ~mean_temp, fill = ~status) %>%
    layer_points(size = 0.1) %>%
    add_relative_scales() %>%
    add_axis("x", title = "month") %>%
    add_axis("x", orient = "top", ticks = 0, title = paste0("Daily mean ", var),
             properties = axis_props(axis = list(stroke = "white"))) %>%
    add_legend(scales = "fill",
               title = "",
               properties = legend_props(
                 legend = list(
                   x = scaled_value("x_rel", 0.45),
                   y = scaled_value("y_rel", -0.15))))
}

#' @export
mvcapa_pump <- function(x, est_band = 2, b = 1, b_point = 1, minsl = 20, nr = 3,
                        vars = "all", diff = FALSE, rank_transform = FALSE,
                        true_anoms = NULL, cost = "cor") {
  n <- nrow(x)
  p <- ncol(x)
  if (rank_transform) x <- apply(x, 2, function(x_i) qnorm(rank(x_i)/(n + 1)))
  x <- anomaly::robustscale(x)
  if (diff) {
    diff_x <- x[2:n, ] - x[1:(n - 1), ]
    Q_hat <- 2 * estimate_precision_mat(diff_x, adjacency_mat(banded_neighbours(est_band, ncol(x))))
  } else
    Q_hat <- estimate_precision_mat(x, adjacency_mat(banded_neighbours(est_band, p)))
  if (cost == "cor") {
    res <- mvcapa_cor(x, Q_hat, b = b, b_point = b_point, min_seg_len = minsl)
    plot_capa(list("x" = x, "anoms" = res), true_anoms = true_anoms)
  } else if (cost == "iid") {
    beta <- iid_penalty(n, p, b)
    beta_tilde <- iid_point_penalty(n, p, b_point)
    res <- anomaly::capa.mv(x,
                            beta        = beta,
                            beta_tilde  = beta_tilde,
                            min_seg_len = minsl,
                            type        = "mean")
    plot_capa(res, true_anoms = true_anoms, cost = "iid")
  }
}

plot_cor_mat <- function(x, est_band, diff = FALSE) {
  n <- nrow(x)
  p <- ncol(x)
  x <- anomaly::robustscale(x)
  Q_hat <- estimate_precision_mat(x, adjacency_mat(banded_neighbours(est_band, p)))
  # Q_hat_diff <- 2 * estimate_precision_mat(diff_x, adjacency_mat(banded_neighbours(est_band, ncol(x))))
  Sigma <- solve(standardise_precision_mat(Q_hat))
  # Sigma_diff <- solve(standardise_precision_mat(Q_hat_diff))
  corrplot(Sigma, method = "number", mar = c(0, 0, 0, 0))
}

save_cor_mat_plot <- function(x, est_band, diff) {
  png("./images/candidate_data_cor_mat.png", width = 6, height = 6, units = "in", res = 800)
  plot_cor_mat(x, est_band)
  dev.off()
}

save_mvcapa_pump_plot <- function(x, cost, b, b_point = 1, minsl = 20, band = 4, true_anoms = NULL) {
  pl <- mvcapa_pump(x, est_band = band, b = b, b_point = b_point, minsl = minsl,
                    true_anoms = true_anoms, cost = cost)
  file_name <- paste0("candidate_data_", cost,
                      "_penscale", b,
                      "_pointpenscale", b_point,
                      "_minsl", minsl,
                      ".png")
  ggsave(paste0("./images/", file_name), width = 8, height = 6, units = "in")
}

promising_plots <- function() {
  x <- make_all_residuals(pump_daily)
  band <- 4
  plot_cor_mat(x, band)

  true_anoms <- pump_anom_inds(pump_daily)
  mvcapa_pump(x, est_band = band, b = 7, true_anoms = true_anoms)
  mvcapa_pump(x, est_band = band, b = 10, true_anoms = true_anoms)
  mvcapa_pump(x, est_band = band, b = 25, true_anoms = true_anoms)
  mvcapa_pump(x, est_band = band, b = 20, true_anoms = true_anoms, b_point = 12, minsl = 20)

  mvcapa_pump(x, est_band = band, b = 10, true_anoms = true_anoms, cost = "iid")
  mvcapa_pump(x, est_band = band, b = 20, true_anoms = true_anoms, cost = "iid")
  mvcapa_pump(x, est_band = band, b = 30, true_anoms = true_anoms, cost = "iid")
  mvcapa_pump(x, est_band = band, b = 30, true_anoms = true_anoms, cost = "iid", b_point = 5, minsl = 20)

  save_mvcapa_pump_plot(x, "cor", 7, 1, 10, band, true_anoms)
  save_mvcapa_pump_plot(x, "cor", 10, 1, 10, band, true_anoms)
  save_mvcapa_pump_plot(x, "cor", 25, 1, 10, band, true_anoms)
  save_mvcapa_pump_plot(x, "cor", 20, 12, 20, band, true_anoms)
  save_mvcapa_pump_plot(x, "iid", 10, 1, 10, band, true_anoms)
  save_mvcapa_pump_plot(x, "iid", 20, 1, 10, band, true_anoms)
  save_mvcapa_pump_plot(x, "iid", 30, 1, 10, band, true_anoms)
}


# temperatures <- c("discharge_temperature",
#                   "mpfm_temperature",
#                   "subsea_barrier_temperature")
# pressures <- c("discharge_pressure",
#                "mpfm_pressure",
#                "subsea_barrier_pressure",
#                "suction_pressure",
#                "differential_pressure")
# other <- c("mpfm_gvf",
#            "mpfm_gas_flow_rate",
#            "mpfm_oil_flow_rate",
#            "mpfm_flow_rate",
#            "output_current",
#            "output_power",
#            "speed")

make_residuals <- function(pump_daily, var, BC = FALSE) {
  covariates <- " ~ mpfm_gas_flow_rate + mpfm_oil_flow_rate + output_current + output_power + speed"
  if (BC) {
    pump_daily[[var]][pump_daily[[var]] <= 0.01] <- 0.01
    pump_daily[[var]] <- predict(caret::BoxCoxTrans(pump_daily[[var]]), pump_daily[[var]])
  }
  model_formula <- formula(paste0(var, covariates))
  lm_obj <- lm(model_formula, data = pump_daily)
  pred_name <- paste0(var, "_give_flow_rate")
  pump_daily[[pred_name]] <- predict(lm_obj, newdata = pump_daily)
  # res_name <- paste0(var, "_residual")
  pump_daily[[var]] - pump_daily[[pred_name]]
  # pump_daily[[res_name]] <- pump_daily[[var]] - pump_daily[[pred_name]]
  # pump_daily[res_name]
}

make_all_residuals <- function(pump_daily, nr = c(3, 4), BC = FALSE) {
  vars <- c("discharge_temperature",
            "mpfm_temperature",
            # "subsea_barrier_temperature",
            "mpfm_pressure",
            "subsea_barrier_pressure",
            "differential_pressure")
  failure_nr <- failure_period(nr)
  pump_daily <- pump_daily %>% filter(date %within% failure_nr)
  do.call("cbind", lapply(vars, make_residuals, pump_daily = pump_daily))
}

pump_anom_inds <- function(pump_daily, nr = 3:4) {
  anoms <- do.call("c", pump_anoms()[nr])
  anoms <- lubridate::date(anoms)
  failure_nr <- failure_period(nr)
  pump_daily <- pump_daily %>% filter(date %within% failure_nr)
  which(pump_daily$date %in% anoms)
}



plot_all_residuals <- function(pump_daily) {
  res <- reshape2::melt(make_all_residuals(pump_daily))
  res$Var2 <- as.factor(res$Var2)
  ggplot2::ggplot(res, ggplot2::aes(x = Var1, y = value, colour = Var2)) +
    ggplot2::geom_line()
}

residual_plot <- function(pump_daily, var, colour_by = "status", nr = c(3, 4),
                          BC = FALSE) {
  failure_nr <- failure_period(nr)
  pump_daily <- pump_daily %>% filter(date %within% failure_nr)
  res_name <- paste0(var, "_residual")
  pump_daily[[res_name]] <- make_residuals(pump_daily, var, nr, BC)
  ggpubr::ggarrange(plot_daily_mean(pump_daily, res_name, "density", colour_by = colour_by),
                    plot_daily_mean(pump_daily, res_name, colour_by = colour_by),
                    plot_daily_mean(pump_daily, var, colour_by = colour_by),
                    nrow = 3, ncol = 1,
                    common.legend = TRUE, legend = "right")
  # pump_daily[[res_name]]
}
