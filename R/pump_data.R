###############################################################
## load libraries
###############################################################
library(tidyverse)
library(lubridate)
library(ggplot2)
library(ggvis)

# glimpse(data)

###############################################
## create new variables
###############################################
load_pump_data <- function() {
  ###############################################################
  ## read data
  ###############################################################
  data <- readRDS(file = "./data/SubseaPumpData/data.Rds")

  data <- data %>%
  dplyr::mutate(differential_pressure = discharge_pressure - suction_pressure,
                running = as_factor(if_else(speed >= 1200, "Running", "Not Running")),
                year = lubridate::year(time),
                month = lubridate::month(time),
                year_month = as.numeric(paste0(lubridate::year(time), str_pad(lubridate::month(time), 2, pad = "0"))),
                date = lubridate::date(time))

  ###############################################
  ## define variable for healthy vs. non-helathy data
  ###############################################
  ## create intervals
  a1 <- lubridate::interval(start = ymd("1970-01-01", tz = "CET"),
                            end = ymd("1970-01-11", tz = "CET"))
  a2 <- lubridate::interval(start = ymd("1970-05-12", tz = "CET"),
                            end = ymd("1970-06-05", tz = "CET"))
  a3 <- lubridate::interval(start = ymd("1970-10-15", tz = "CET"),
                            end = ymd("1972-12-14", tz = "CET"))
  a4 <- lubridate::interval(start = ymd("1973-11-09", tz = "CET"),
                            end = ymd("1978-03-29", tz = "CET"))
  a5 <- lubridate::interval(start = ymd("1979-03-01", tz = "CET"),
                            end = ymd("1980-06-04", tz = "CET"))
  a <- list(a1, a2, a3, a4, a5)
  ## define variable
  data %>%
    dplyr::mutate(status = as_factor(if_else(time %within% a, "Healthy", "Not Healthy")))
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
                 "subsea_barrier_temperature",
                 "suction_pressure",
                 "differential_pressure")
  other <- c("mpfm_gvf",
             "mpfm_flow_rate",
             "output_current",
             "output_power",
             "speed")
  if (type == "temperatures") return(temperatures)
  else if (type == "pressures") return(pressures)
  else return(c(temperatures, pressures, other))
}

daily_means <- function(data, vars = "all") {
  data %>%
    dplyr::filter(running == "Running") %>%
    dplyr::select(c("date", "status", relevant_vars(vars))) %>%
    dplyr::group_by(date, status) %>%
    dplyr::summarise_all(mean)
}

plot_daily_mean <- function(data, var) {
  data %>%
    dplyr::filter(running == "Running") %>%
    dplyr::select(date, starts_with(var), status) %>%
    dplyr::group_by(date, status) %>%
    dplyr::summarise(mean_temp = mean(eval(parse(text = var)))) %>%
    ggvis(x = ~date, y = ~mean_temp, fill = ~status) %>%
    layer_points(size = 0.1) %>%
    add_relative_scales() %>%
    add_axis("x", title = "date") %>%
    add_axis("x", orient = "top", ticks = 0, title = paste0("Daily mean ", var),
             properties = axis_props(axis = list(stroke = "white"))) %>%
    add_legend(scales = "fill",
               title = "",
               properties = legend_props(
                 legend = list(
                   x = scaled_value("x_rel", 0.45),
                   y = scaled_value("y_rel", -0.15))))
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

failure_period <- function(nr) {
  a1 <- lubridate::interval(start = ymd("1970-01-01", tz = "CET"),
                            end = ymd("1970-04-12", tz = "CET"))
  a2 <- lubridate::interval(start = ymd("1970-05-12", tz = "CET"),
                            end = ymd("1970-06-14", tz = "CET"))
  a3 <- lubridate::interval(start = ymd("1970-07-27", tz = "CET"),
                            end = ymd("1973-11-08", tz = "CET"))
  a4 <- lubridate::interval(start = ymd("1973-11-09", tz = "CET"),
                            end = ymd("1979-01-14", tz = "CET"))
  a <- list(a1, a2, a3, a4)
  a[nr]
}

mvcapa_pump <- function(est_band = 2, b = 1, b_point = 1, nr = 3, vars = "all", diff = FALSE, rank_transform = FALSE) {
  failure_nr <- failure_period(nr)
  daily_pump_data <- daily_means(data, vars) %>% filter(date %within% failure_nr)
  x <- as.matrix(daily_pump_data[, 3:ncol(daily_pump_data)])
  print(apply(x, 2, median))
  if (!diff) {
    if (rank_transform) x <- apply(x, 2, function(x_i) qnorm(rank(x_i)/(nrow(x) + 1)))
    x <- anomaly::robustscale(x)
    diff_x <- x[2:nrow(x), ] - x[1:(nrow(x) - 1), ]
    Q_hat <- estimate_precision_mat(x, adjacency_mat(banded_neighbours(est_band, ncol(x))))
    # Q_hat <- 2 * estimate_precision_mat(diff_x, adjacency_mat(banded_neighbours(est_band, ncol(x))))
    print(Q_hat)
    print(solve(Q_hat))
    res <- mvcapa_cor(x, Q_hat, b = b, b_point = b_point, min_seg_len = 5)
    collective_anoms <- collective_anomalies(list("anoms" = res))
    print(unique(collective_anoms$start))
    print(daily_pump_data$date[unique(collective_anoms$start)])
    plot_capa(list("x" = x, "anoms" = res))
  } else {
    diff_x <- x[2:nrow(x), ] - x[1:(nrow(x) - 1), ]
    diff_x <- anomaly::robustscale(diff_x)
    Q_hat <- estimate_precision_mat(diff_x, adjacency_mat(banded_neighbours(est_band, ncol(x))))
    print(Q_hat)
    res <- mvcapa_cor(diff_x, Q_hat, b = b, b_point = b_point, min_seg_len = 5)
    plot_capa(list("x" = diff_x, "anoms" = res))
  }
}
