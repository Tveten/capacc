setup_parallel <- function() {
  n_cores <- parallel::detectCores()
  c <- parallel::makeCluster(n_cores, outfile = '', type = 'PSOCK')
  doParallel::registerDoParallel(c)
  c
}

stop_parallel <- function(c) {
  parallel::stopCluster(c)
}

extract_nested_element <- function(i, list_of_lists) {
  lapply(list_of_lists, function(list_in_list) list_in_list[[i]])
}
