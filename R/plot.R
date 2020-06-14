
capa_line_plot <- function(object, epoch = dim(object$x)[1],
                           subset = 1:ncol(object$x), variate_names = NULL,
                           true_anoms = NULL, cost = "cor") {
    # creating null entries for ggplot global variables so as to pass CRAN checks
    x <- value <- ymin <- ymax <- x1 <- x2 <- y1 <- y2 <- x1 <- x2 <- y1 <- y2 <- NULL
    if (cost == "cor") x <- object$x
    else if (cost == "iid") x <- object@data
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
    out <- ggplot2::ggplot(data = data_df)
    out <- out + ggplot2::aes(x = x, y = value)
    out <- out + ggplot2::geom_point(size = 0.4)

    # highlight the collective anomalies
    if (cost == "cor") c_anoms <- collective_anomalies(object)
    else if (cost == "iid") c_anoms <- anomaly::collective_anomalies(object)
    c_anoms <- c_anoms[c_anoms$variate %in% subset, ]
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
                                      alpha=0.5) +
        ggplot2::scale_fill_manual(name = "",
                                   values = c("Estimated collective anomaly" = "red"))
      # values = c("Known anomaly" = "green",
      #            "Estimated anomaly" = "blue"))
    }

    # highlight the point anomalies
    if (cost == "cor") p_anoms <- point_anomalies(object)
    else if (cost == "iid") p_anoms <- anomaly::point_anomalies(object)
    p_anoms <- p_anoms[p_anoms$variate %in% subset,]
    if(!any(is.na(p_anoms)) & nrow(p_anoms) > 0)
    {
      p_anoms_data_df <- Reduce(rbind, Map(function(a,b) {
        data_df[data_df$variable == names[a] & data_df$x == b, ]
      }, p_anoms$variate, p_anoms$location)
      )
      out <- out + ggplot2::geom_point(data = p_anoms_data_df,
                                       ggplot2::aes(colour = "Estimated point anomaly"),
                                       size = 0.8) +
        ggplot2::scale_colour_manual(name = "",
                                     values = c("Estimated point anomaly" = "red"))

    }

    out <- out + ggplot2::facet_grid(factor(variable, levels = names) ~ .,
                                     scales = "free_y")

    # Highlight true anomalies if specified.
    if (!is.null(true_anoms)) {
      n_true_anoms <- length(true_anoms) / 2
      even <- 2 * 1:n_true_anoms
      odd <- even - 1
      true_anoms_df <- data.frame("variable" = names[rep(1:p, n_true_anoms)],
                                  "start"    = rep(true_anoms[odd], each = p),
                                  "end"      = rep(true_anoms[even], each = p))
      true_anoms_df$ymin <- -Inf
      true_anoms_df$ymax <- Inf
      out <- out + ggplot2::geom_rect(data = true_anoms_df,
                                      inherit.aes = FALSE,
                                      mapping = ggplot2::aes(xmin = start,
                                                             xmax = end,
                                                             ymin = ymin,
                                                             ymax = ymax,
                                                             linetype = "Known collective anomaly"),
                                      alpha = 0,
                                      size = 1,
                                      colour = "green3") +
        ggplot2::scale_linetype_manual(name = "",
                                       values = c("Known collective anomaly" = 1))
    }

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
                                legend.position = legend_position
                              ) +
      ggplot2::xlab("t") +
      ggplot2::ylab("Value")

    return(out)
}

capa_tile_plot <- function(object, variate_names = FALSE,
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
    if(variate_names == FALSE)
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

plot_capa <- function(object, epoch = n,
                      subset = 1:p, variate_names = NULL,
                      true_anoms = NULL, cost = "cor") {
  if (cost == "cor") {
    n <- nrow(object$x)
    p <- ncol(object$x)
  } else if (cost == "iid") {
    n <- nrow(object@data)
    p <- ncol(object@data)
  }
  if (length(subset) <= 30)
    return(capa_line_plot(object, epoch, subset, variate_names, true_anoms, cost))
  else
    return(capa_tile_plot(object, variate_names, epoch, subset, true_anoms))
}

