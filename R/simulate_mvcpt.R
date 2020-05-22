#'
#' rescale_variance <- function(x, by_col = TRUE) {
#'   if (by_col) {
#'     for (j in 1:ncol(x)){
#'       scale <- mad(diff(x[, j])) / sqrt(2)
#'       x[, j] <- x[, j] / scale
#'     }
#'   } else {
#'     for (j in 1:ncol(x)){
#'       scale <- mad(diff(x[j, ])) / sqrt(2)
#'       x[j, ] <- x[j, ] / scale
#'     }
#'   }
#'   x
#' }
#'
#' #' @export
#' simulate_mvcpt <- function(data = init_data(n = 200, p = 10, rho = 0.8, band = 2,
#'                                             locations = 80, durations = 120),
#'                            params = method_params(), seed = NULL) {
#'   get_adj_mat <- function(est_struct) {
#'     if (est_struct == "correct")
#'       return(data$Sigma_inv)
#'     else if (est_struct == "banded")
#'       return(adjacency_mat(banded_neighbours(params$est_band, data$p)))
#'   }
#'
#'   if (!is.null(seed)) set.seed(seed)
#'   x <- rescale_variance(simulate_cor_(data))
#'   diff_x <- x[2:nrow(x), ] - x[1:(nrow(x) - 1), ]
#'   Q_hat <- 2 * estimate_precision_mat(diff_x, get_adj_mat(params$precision_est_struct))
#'   if (grepl("inspect", params$cost))
#'     return(single_cor_inspect(t(x), Q_hat))
#'   else if (grepl("mvlrt", params$cost))
#'     return(single_mvnormal_changepoint(x, Q_hat,
#'                                        b = params$b,
#'                                        min_seg_len = params$minsl))
#' }
#'
