
profile_capacc <- function(x, Q, n_sim = 100) {
  p <- ncol(Q)
  nbs <- lower_nbs(Q)
  extended_nbs <- extended_lower_nbs(nbs)
  sparse_penalty <- get_penalty('sparse', 1000, p)
  dense_penalty <- get_penalty('dense', 1000, p)
  start_profiler("./profile.out")
  for (i in 1:n_sim) optimise_mvnormal_saving(x, Q, nbs, extended_nbs,
                                              dense_penalty$alpha,
                                              sparse_penalty$beta,
                                              sparse_penalty$alpha)
  stop_profiler()
}
