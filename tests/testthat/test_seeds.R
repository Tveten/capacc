test_seeds <- function(variables) {
  data <- init_data()
  method <- method_params()

  params_list <- split_params(
    expand_list(c(data, method), variables),
    list("data"   = names(data),
         "method" = names(method))
  )
  get_sim_seeds(params_list, variables)
}

test_that('Seeds are generated correctly for multiple simulation runs', {
  variables1 <- list("cost"        = c("iid", "cor"),
                     "precision_est_struct" = c(NA, "banded"),
                     "est_band"    = c(1, 2))
  expect_true(length(unique(test_seeds(variables1))) == 1)
  variables2 <- list("cost"        = c("iid", "cor"),
                     "precision_est_struct" = c(NA, "banded"),
                     "est_band"    = c(1, 2),
                     "rho" = c(0.5, 0.9))
  expect_true(length(unique(test_seeds(variables2))) == 2)
  variables3 <- list("cost"        = c("iid", "cor"),
                     "est_band"    = c(1, 2),
                     "rho" = c(0.5, 0.9))
  expect_true(length(unique(test_seeds(variables3))) == 2)
  variables4 <- list("cost"        = c("iid", "cor"),
                     "est_band"    = c(1, 2),
                     "rho" = c(0.5, 0.9),
                     "proportions" = c(0.1, 0.2, 0.3, 0.4))
  expect_true(length(unique(test_seeds(variables4))) == 8)
  variables5 <- list("cost"        = c("iid", "cor"),
                     "rho"         = c(0.5, 0.9),
                     "est_band"    = c(1, 2))
  expect_error(test_seeds(variables5))
})
