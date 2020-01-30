load_as_dt <- function(file_name) {
  file_name <- paste0('./data/', file_name)
  fread(file_name)
}

# Water quality data.
load_water_quality_data <- function() {
  water_quality_data <- load_as_dt('beach_water_quality_automated_sensors_1.csv')
  water_quality_data[, 'measurement_timestamp' := anytime::anytime(measurement_timestamp)]
  water_quality_data
}

plot_wq_data <- function(y = 2018, comp = 'water_temperature') {
  wq_data <- load_water_quality_data()
  comp <- parse(text = comp)
  ggplot2::ggplot(wq_data[year(measurement_timestamp) == y],
                  ggplot2::aes(x = measurement_timestamp,
                               y = eval(comp),
                               color = beach_name)) +
    ggplot2::geom_line()
}

# Shuttle data.
load_shuttle_data <- function() {
  shuttle_data <- load_as_dt('shuttle-unsupervised-ad.csv')
  shuttle_data[, 'ind' := 1:.N]
  # water_quality_data[, 'measurement_timestamp' := anytime::anytime(measurement_timestamp)]
  shuttle_data <- as.data.table(tidyr::gather(shuttle_data, 'var', 'value', -c(ind, V10)))
  shuttle_data[, 's_value' := (value - mean(value)) / sd(value), var]
  shuttle_data
}

plot_data <- function(dataset = 'shuttle', comp = 'V1') {
  get_dataset <- function(dataset) {
    if (dataset == 'shuttle') return(load_shuttle_data())
  }

  dataset <- get_dataset(dataset)
  dataset <- dataset[var == comp]
  ggplot2::ggplot(dataset, ggplot2::aes(x = ind, y = s_value, color = var)) +
    ggplot2::geom_line()
}

# Temperature data.
load_global_temperature_data <- function() {
  temp_data <- load_as_dt('climate-change-earth-surface-temperature-data/GlobalLandTemperaturesByCity.csv')
  temp_data[, 'dt' := as.Date(dt)]
  temp_data[, 'CityCountry' := paste0(City, ', ', Country)]
  temp_data
}

sub_temp_data <- function(temp_data, country, city = 'all', month = 'all', from_year = 1700) {
  get_year <- function(date) {
    date <- as.character(date)
    unlist(extract_nested_element(1, strsplit(date, '-')))
  }

  get_month <- function(date) {
    date <- as.character(date)
    unlist(extract_nested_element(2, strsplit(date, '-')))
  }

  sub_data <- temp_data[vapply(Country, function(x) any(x == country), logical(1))]
  if (city != 'all') sub_data <- sub_data[City == city]
  sub_data <- sub_data[get_year(dt) >= from_year]
  if (month != 'all') sub_data <- sub_data[get_month(dt) == month]
  sub_data
}

plot_temp_data <- function(temp_data, city = 'Oslo', country = 'Norway', month = 'all', from_year = 1700) {
  sub_data <- sub_temp_data(temp_data, country, city, month, from_year)
  ggplot2::ggplot(sub_data, ggplot2::aes(x = dt, y = AverageTemperature)) +
    ggplot2::geom_line() +
    ggplot2::ggtitle(paste0(city, ', ', month))
}

country_video <- function(temp_data, country = 'Norway', month = 'all', from_year = 1700) {
  # Noteworthy:
  #   - Sweden (in combination with rest of Scandinavia?)
  #   - Italy in August.
  #
  # Compare months. Which months are the most anomalous?

  sub_data <- temp_data[Country == country]
  cities <- sub_data[, unique(City)]
  for (i in 1:length(cities)) {
    print(paste0(cities[i], ': ', i, ' of ', length(cities), '.'))
    show(plot_temp_data(sub_data, cities[i], country, month, from_year))
    Sys.sleep(2)
  }
}

noteworthy_temp_plots <- function() {
  plots <- list(plot_temp_data('Reykjavík', from_year = 1800, month = '01'),
                plot_temp_data('Peking', from_year = 1800, month = '01'),
                plot_temp_data('Peking', from_year = 1800, month = '07')
  )
  for(plot in plots) {
    show(plot)
    Sys.sleep(2)
  }
}

avg_temp_cols <- function(wide_temp_data) {
  split_col_names <- strsplit(names(wide_temp_data), '[.]')
  first_elems <- unlist(extract_nested_element(1, split_col_names))
  avg_temp_inds <- first_elems == 'AverageTemperature'
  avg_temp_data <- wide_temp_data[, avg_temp_inds, with = FALSE]
  names(avg_temp_data) <- unlist(extract_nested_element(2, split_col_names[avg_temp_inds]))
  avg_temp_data
}

mvcapa_cor_temp <- function(temp_data, country = 'Norway', city = 'all', month = '01', from_year = 1770,
                            lambda_min_ratio = 0.1, b = 1, min_seg_len = 3, max_seg_len = 50) {
  # Anomalies detected:
  #   - Switzerland august anomaly.

  sub_data <- sub_temp_data(temp_data, country, city, month, from_year)
  wide_sub_data <- reshape(sub_data[, .(dt, CityCountry, AverageTemperature)],
                           direction = 'wide', timevar = 'CityCountry',
                           idvar = 'dt', v.names = 'AverageTemperature')
  x <- as.matrix(avg_temp_cols(wide_sub_data))
  res <- mvcapa_cor(x, lambda_min_ratio = lambda_min_ratio, b = b,
                    min_seg_len = min_seg_len, max_seg_len = max_seg_len)
  plot_capa(res)
}
