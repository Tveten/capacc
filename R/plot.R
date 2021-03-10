
capa_line_plot <- function(object, epoch = dim(object$x)[1],
                           subset = 1:ncol(object$x), variate_names = NULL,
                           true_anoms = NULL) {
    # creating null entries for ggplot global variables so as to pass CRAN checks
    x <- value <- ymin <- ymax <- x1 <- x2 <- y1 <- y2 <- x1 <- x2 <- y1 <- y2 <- NULL
    x <- object$x
    data_df <- as.data.frame(x)
    # names <- paste("y", 1:ncol(x), sep = "")
    if (is.null(variate_names)) {
      names <- as.character(1:ncol(x))
    } else names <- variate_names
    colnames(data_df) <- names
    data_df <- as.data.frame(data_df[, subset, drop = FALSE])
    n <- nrow(data_df)
    p <- ncol(data_df)
    data_df <- cbind(data.frame("x" = 1:n), data_df)
    data_df <- reshape2::melt(data_df, id = "x")
    out <- ggplot2::ggplot(data = data_df, ggplot2::aes(x = x, y = value)) +
      ggplot2::geom_point(size = 0.4)

    c_anoms <- collective_anomalies(object)
    p_anoms <- point_anomalies(object)
    c_anoms <- c_anoms[c_anoms$variate %in% subset, ]
    p_anoms <- p_anoms[p_anoms$variate %in% subset,]

    # Highlight the collective anomalies
    est_col <- "red"
    if(!any(is.na(c_anoms)) & nrow(c_anoms) > 0)
    {
      c_anoms_data_df <- c_anoms[, 1:3]
      names(c_anoms_data_df) <- c(names(c_anoms_data_df)[1:2], "variable")
      c_anoms_data_df$variable <- names[c_anoms_data_df$variable]
      c_anoms_data_df$ymin <- -Inf
      c_anoms_data_df$ymax <- Inf
      out <- out + ggplot2::geom_rect(data = c_anoms_data_df,
                                      inherit.aes = FALSE,
                                      mapping = ggplot2::aes(xmin = start,
                                                             xmax = end,
                                                             ymin = ymin,
                                                             ymax = ymax,
                                                             fill = "Estimated collective anomaly"),
                                      # colour = est_col,
                                      alpha=0.5)
    }

    # Highlight the point anomalies
    if(!any(is.na(p_anoms)) & nrow(p_anoms) > 0)
    {
      p_anoms_data_df <- Reduce(rbind, Map(function(a,b) {
        data_df[data_df$variable == names[a] & data_df$x == b, ]
      }, p_anoms$variate, p_anoms$location)
      )
      out <- out + ggplot2::geom_point(data = p_anoms_data_df,
                                       ggplot2::aes(colour = "Estimated point anomaly"),
                                       size = 0.9)
    }

    # Highlight true anomalous observations on x-axis if specified.
    if (!is.null(true_anoms)) {
      anoms <- unlist(lapply(1:nrow(true_anoms), function(i) (true_anoms$start[i] + 1):true_anoms$end[i]))
      true_anoms <- data.frame(x = 1:n, anom = "Known normal", variable = names[length(names)], stringsAsFactors = FALSE)
      true_anoms$anom[anoms] <- "Known anomaly"
      out <- out + ggplot2::geom_rug(data = true_anoms[true_anoms$anom == "Known anomaly", ], inherit.aes = FALSE,
                                     sides = "b", outside = TRUE, length = unit(0.5, "cm"),
                                     # size = 0.8,
                                     colour = "dodgerblue3",
                                     mapping = ggplot2::aes(x = x, linetype = anom)) +
        ggplot2::coord_cartesian(clip = "off")
    }

    # Highlight estimated anomalous observations on x-axis.
    starts <- unique(c_anoms$start)
    ends <- unique(c_anoms$end)
    point_locs <- unique(p_anoms$location)
    if (length(c(starts, ends, point_locs)) > 0) {
      all_anom_est_x <- c(point_locs, unlist(lapply(1:length(starts), function(i) {
        starts[i]:ends[i]
      })))
      est_anoms <- data.frame(x = 1:n, anom = "Estimated normal", variable = names[length(names)], stringsAsFactors = FALSE)
      est_anoms$anom[all_anom_est_x] <- "Estimated anomaly"
      out <- suppressMessages(out + ggplot2::geom_rug(data = est_anoms[est_anoms$anom == "Estimated anomaly", ],
                                                      inherit.aes = FALSE,
                                     sides = "b", outside = TRUE, length = unit(0.5 - 0.5 / 1.618, "cm"),
                                     # size = 0.3,
                                     colour = est_col,
                                     mapping = ggplot2::aes(x = x)) +
        ggplot2::coord_cartesian(clip = "off"))
    }

    out <- out +
      ggplot2::scale_colour_manual(name = "",
                                   values = c("Estimated point anomaly" = est_col)) +
      ggplot2::scale_fill_manual(name = "",
                                 values = c("Estimated collective anomaly" = est_col)) +
      ggplot2::scale_linetype_manual(name = "",
                                     values = c("Known normal"  = "blank",
                                                "Known anomaly" = "solid",
                                                "Estimated normal"  = "blank",
                                                "Estimated anomaly" = "solid")) +
      ggplot2::facet_grid(factor(variable, levels = names) ~ .,
                                 scales = "free_y")


    # grey out the data after epoch
    if(epoch != nrow(x))
    {
      d <- data.frame(variable = names[subset], x1 = epoch, x2 = n,
                      y1 = -Inf, y2 = Inf)
      out <- out + geom_rect(data = d, inherit.aes = FALSE,
                               mapping = aes(xmin = x1, xmax = x2,
                                             ymin = y1, ymax = y2),
                               fill = "yellow", alpha = 0.2)
    }
    # change background
    if (!is.null(true_anoms)) legend_position <- "bottom"
    else legend_position <- "none"
    out <- out + ggplot2::theme(
                                # Hide panel borders and remove grid lines
                                panel.border = ggplot2::element_blank(),
                                panel.grid.major = ggplot2::element_blank(),
                                panel.grid.minor = ggplot2::element_blank(),
                                # Change axis line
                                axis.line = ggplot2::element_line(colour = "black"),
                                axis.ticks.length.x=unit(0.65, "cm"),
                                legend.position = legend_position
                              ) +
      ggplot2::xlab("t") +
      ggplot2::ylab("Value")
    if (length(subset) > 8)
      out <- out + ggplot2::scale_y_continuous(breaks = NULL)

    return(out)
}

capa_tile_plot <- function(object, variate_names = NULL,
                           epoch = dim(object$x)[1], subset = 1:ncol(object$x),
                           true_anoms = NULL) {
    # nulling out variables used in ggplot to get the package past CRAN checks
    x1 <- y1 <- x2 <- y2 <- variable <- value <- NULL
    df <- as.data.frame(object$x)
    df <- as.data.frame(df[, subset, drop = FALSE])
    # normalise data
    for(i in 1:ncol(df))
    {
        df[, i] <- (df[, i] - min(df[, i])) / (max(df[, i]) - min(df[, i]))
    }
    n <- data.frame("n" = seq(1, nrow(df)))
    molten.data <- reshape2::melt(cbind(n, df), id = "n")
    out <- ggplot2::ggplot(molten.data, ggplot2::aes(n, variable))
    out <- out + ggplot2::geom_tile(ggplot2::aes(fill = value))
    # get any collective anomalies
    c_anoms <- collective_anomalies(object)
    c_anoms <- c_anoms[c_anoms$variate %in% subset, ]
    c_anoms <- unique(c_anoms[, 1:2])
    if(!any(is.na(c_anoms)) & nrow(c_anoms) > 0)
    {
        ymin <- 0
        ymax <- ncol(df)
        out <- out + ggplot2::annotate("rect", xmin = c_anoms$start,
                                       xmax = c_anoms$end, ymin = ymin,
                                       ymax = ymax + 1, alpha = 0.0,
                                       color = "red", fill = "yellow")
    }
    if(epoch != nrow(object$x))
    {
        d <- data.frame(x1 = epoch, x2 = nrow(object$x), y1 = -Inf, y2 = Inf)
        out <- out + ggplot2:::geom_rect(data = d, inherit.aes = FALSE,
                                         mapping = ggplot2::aes(xmin = x1,
                                                                xmax = x2,
                                                                ymin = y1,
                                                                ymax = y2),
                                         fill = "yellow", alpha=0.2)
    }
    if(is.null(variate_names))
    {
        out <- out + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                    axis.title = ggplot2::element_blank())
    }
    out <- out + ggplot2::theme(
                                # Hide panel borders and remove grid lines
                                panel.border = ggplot2::element_blank(),
                                panel.grid.major = ggplot2::element_blank(),
                                panel.grid.minor = ggplot2::element_blank(),
                                # Change axis line
                                axis.line = ggplot2::element_line(colour = "black"),
                                )
    return(out)
}

#' Plotting function for capacc objects
#'
#' @param object A capacc object to be plotted.
#' @param epoch Plots data from 1:epoch.
#' @param subset Plots variables in subset.
#' @param variate_names An optional vector with the names of the variables.
#' @param true_anoms An optional data frame with columns "start" and "end" to illustrate known locations of anomalies.
#'
#' @return A ggplot of the detected collective and point anomalies.
#'
#' @export
plot.capacc <- function(object, epoch = n,
                        subset = 1:p, variate_names = NULL,
                        true_anoms = NULL) {
# plot.capacc <- function(object, epoch = n,
#                         subset = 1:p, variate_names = NULL,
#                         true_anoms = NULL, cost = "cor") {
  # if (cost == "cor") {
  n <- nrow(object$x)
  p <- ncol(object$x)
  # } else if (cost == "iid") {
  #   n <- nrow(object@data)
  #   p <- ncol(object@data)
  # }
  if (length(subset) <= 30)
    return(capa_line_plot(object, epoch, subset, variate_names, true_anoms))
  else
    return(capa_tile_plot(object, variate_names, epoch, subset, true_anoms))
}

