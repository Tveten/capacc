setup_parallel <- function(cpus = 1) {
  if (is.null(cpus)) cpus <- parallel::detectCores() - 1
  c <- parallel::makeCluster(cpus, outfile = '', type = 'PSOCK')
  doParallel::registerDoParallel(c)
  c
}

stop_parallel <- function(c) {
  parallel::stopCluster(c)
}

extract_nested_element <- function(i, list_of_lists) {
  lapply(list_of_lists, function(list_in_list) list_in_list[[i]])
}

extract_nested_elements <- function(i, list_of_lists) {
  lapply(list_of_lists, function(list_in_list) list_in_list[i])
}

expand_list <- function(a, vars) {
  var_grid <- expand.grid(vars, stringsAsFactors = FALSE)
  var_names <- names(var_grid)

  # For loop to avoid "c" being used on each row of var_grid, which might change
  # the type of variables.
  out_list <- list()
  for (i in 1:nrow(var_grid)) {
    a_copy <- a
    for (j in 1:ncol(var_grid)) {
      a_copy[[var_names[j]]] <- unname(var_grid[i, j])
    }
    out_list[[length(out_list) + 1]] <- a_copy
  }
  out_list
}

split_lists <- function(a, grouping) {
  split_list <- lapply(grouping, extract_nested_elements, list_of_lists = a)
  names(split_list) <- names(grouping)
  split_list
}

combine_lists <- function(a) {
  Reduce(function(x, y) Map(c, x, y), a)
}

identity <- function(a) {return(a)}

off_diag_ind <- function(nr, p) {
  A <- matrix(FALSE, nrow = p, ncol = p)
  delta <- row(A) - col(A)
  for (n in nr) {
    A[delta == n] <- TRUE
  }
  A
}

off_diag <- function(nr, x) {
  x[off_diag_ind(nr, nrow(x))]
}

count_nz_columns <- function(x) {
  m <- nrow(x)
  sum(apply(x, 2, function(u) isTRUE(all.equal(u, rep(0, m)))))
}

is_zero <- function(u) {
  vapply(u, function(j) isTRUE(all.equal(j, 0)), logical(1))
}

indicator <- function(J, p) {
  x <- rep(0, p)
  x[J] <- 1
  x
}

inds_from_intervals <- function(starts, ends, n) {
  if (length(starts) != length(ends)) stop("starts and ends must be of the same length")
  n_intervals <- length(starts)
  inds <- rep(FALSE, n)
  if (is.na(starts) && length(starts) == 1) return(inds)
  else {
    interval_inds <- do.call("c", lapply(1:n_intervals, function(i) {
      if (!is.na(starts[i]) && !is.na(ends[i]))
        return(starts[i]:ends[i])
      else return(integer(0))
    }))
    inds[interval_inds] <- TRUE
    return(inds)
  }
}

dot_every <- function(n, f) {
  i <- 1
  function(...) {
    if (i %% n) cat(".")
    i <<- i + 1
    f(...)
  }
}

print_iterations <- function(f) {
  i <- 1
  function(...) {
    print(paste0("Iteration ", i))
    i <<- i + 1
    f(...)
  }
}

print_progress <- function(f, end) {
  i <- 1
  function(...) {
    print(paste0("Progress: ", round(i / end, 3) * 100, "%."))
    i <<- i + 1
    f(...)
  }
}

adjacent_dist <- function(x) {
  do.call("c", lapply(2:nrow(x), function(i) {
    as.numeric(dist(x[(i - 1):i, ]))
  }))
}

is_equal <- function(x_vec, y) {
  vapply(x_vec, function(x) isTRUE(all.equal(x, y)), logical(1))
}
