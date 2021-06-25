#' @export
simple_example <- function(p = 4, n = 200, vartheta = 10, method = "cor",
                           b = 2, b_point = 0.5, v = 1, seed = NULL, bw = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  if (v == 1) {
    true_anoms <- data.frame(start = c(20, n - 50), duration = c(round(n/10), 30))
    true_anoms$end <- true_anoms$start + true_anoms$duration
    point_locs <- c(round(n / 2)) + c(-5, 0, 5)
    x <- simulate_cor(n = n, p = p,
                      locations = true_anoms$start,
                      durations = true_anoms$duration,
                      proportions = c(2/p, 1/p),
                      vartheta = vartheta, change_type = "random",
                      point_locations = point_locs, point_proportions = 1/p,
                      point_mu = c(1.5 * vartheta, - vartheta, vartheta),
                      Sigma = car_precision_mat(banded_neighbours(3, p), rho = 0.9))$x
  } else if (v == 2) {
    true_anoms <- data.frame(start = c(50, round(n / 3), n - 100),
                             duration = c(round(n / 20), round(n / 40), round(n / 10)))
    true_anoms$end <- true_anoms$start + true_anoms$duration
    point_locs <- c(200, 210, 215, 260, 400, 406, 460, 630, 700, 705, 712, 800)
    x <- simulate_cor(n = 1000, p = p,
                      locations = true_anoms$start,
                      durations = true_anoms$duration,
                      proportions = c(2/p, 1, 1/p),
                      vartheta = c(0.8, 4, 0.8), change_type = "random",
                      point_locations = point_locs,
                      point_mu = 4, point_proportions = 2/p,
                      Sigma = constant_cor_mat(p, 0.5))$x
  }
  true_anoms <- rbind(true_anoms, data.frame(start = point_locs - 1, duration = 1, end = point_locs))
  if (method == "cor") {
    x <- centralise(x)
    Q_hat <- robust_sparse_precision(x, adjacency_mat(banded_neighbours(4, p), sparse = FALSE))
    res <- capacc(x, Q_hat, b = b, b_point = b_point)
    # Function from pump_data.R
    if (bw) plot2.capacc(list("x" = x, "anoms" = res), true_anoms = true_anoms, bw = TRUE)
    else plot_capa(list("x" = x, "anoms" = res), true_anoms = true_anoms)
  } else if (method == "iid") {
    beta <- iid_penalty(n, p, b)
    beta_tilde <- iid_point_penalty(n, p, b_point)
    res <- anomaly::capa.mv(x,
                            beta        = beta,
                            beta_tilde  = beta_tilde,
                            type        = "mean")
    if (bw) plot2.capacc(res, cost = "iid", true_anoms = true_anoms, bw = TRUE)
    else plot_capa(res, cost = "iid", true_anoms = true_anoms)
    # plot_capa(res, cost = "iid", true_anoms = true_anoms)
  }
}

save_simple_example <- function(p = 4, n = 200, vartheta = 10,
                                b = 2, b_point = 2, seed = NULL, bw = FALSE) {
  file_name <- paste0("simple_example",
                      "_p", p,
                      "_vartheta", vartheta,
                      "_b", b,
                      "_bpoint", b_point)
  if (bw) file_name <- paste0(file_name, "_bw")
  file_name <- paste0(file_name, ".png")
  show(simple_example(p, n, vartheta, b, b_point, seed, bw = bw))
  ggsave(paste0("./images/", file_name), width = 7, height = 5, units = "in")
}

simple_examples_presentation <- function() {
  save_simple_example(vartheta = 0, b = 10^6, b_point = 10^6, seed = 6)
  save_simple_example(b = 10^6, b_point = 10^6, seed = 6)
  save_simple_example(seed = 6)
}

simple_examples_paper <- function(p = 10, save = FALSE, bw = FALSE) {
  plots <- list(simple_example(p = p, n = 1000, method = "cor", b = 1, b_point = 0.5,
                               vartheta = 1, v = 2, seed = 103, bw = bw),
                simple_example(p = p, n = 1000, method = "iid", b = 1.68, b_point = 0.5,
                               vartheta = 1, v = 2, seed = 103, bw = bw))
  figure <- ggpubr::ggarrange(plotlist = plots, nrow = 1, ncol = 2,
                              common.legend = TRUE,
                              legend = "bottom")
  if (save) {
    show(figure)
    file_name <- paste0("./images/simple_example", p)
    if (bw) file_name <- paste0(file_name, "_bw")
    file_name <- paste0(file_name, ".png")
    ggsave(file_name, width = 8, height = 5, units = "in")
  } else return(figure)
}

