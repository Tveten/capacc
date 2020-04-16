load_price_data <- function(type = "all") {
  if (type == "business_services")
    file_name <- "oslo_nyse_nasdaq____business_services__2020-03-03.parquet"
  else if (type == "consumer_cyclicals")
    file_name <- "oslo_nyse_nasdaq____consumer_cyclicals__2020-03-03.parquet"
  else if (type == "consumer_services")
    file_name <- "oslo_nyse_nasdaq____consumer_services__2020-03-03.parquet"
  else if (type == "energy")
    file_name <- "oslo_nyse_nasdaq____energy__2020-03-03.parquet"
  else if (type == "industrials")
    file_name <- "oslo_nyse_nasdaq____industrials__2020-03-03.parquet"
  else
    file_name <- "oslo_nyse_nasdaq____all__2020-03-04.parquet"

  print(file_name)
  arrow::read_parquet(paste0("./data/Prisdata/", file_name), as_tibble = TRUE)
}
