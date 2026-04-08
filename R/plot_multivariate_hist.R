#' Prepare data for visualizing a multivariate histogram in two chosen dimensions
#'
#' @description
#' Prepares the data needed to visualise a multivariate adaptive histogram
#' (as returned by \code{\link{build_adaptive_histogram}} or \code{\link{BuildHist}}) in two chosen
#' dimensions.
#'
#' @param X A numeric matrix of dimension \eqn{n \times d}, the original data.
#' @param hist A numeric matrix of dimension \eqn{m \times (2d + 5)} as
#'   returned by \code{\link{build_adaptive_histogram}} or \code{\link{BuildHist}}.
#' @param lower_constraint A numeric vector of length \eqn{d} specifying the
#'   lower bound for each coordinate. If provided,
#'   computes the conditional distribution of \eqn{X} within the lower and upper constraints.
#'   Otherwise, computes the marginal distribution in \code{plot_coord}.
#' @param upper_constraint A numeric vector of length \eqn{d} specifying the
#'   upper bound for each coordinate. Analogously overridden for
#'   \code{plot_coord}.
#' @param plot_coord An integer vector of length 2 giving the column indices
#'   of \code{X} to use as the x- and y-axes of the plot.
#' @param nint A positive integer. The number of intervals along each plotting
#'   axis. The grid contains \code{nint^2} cells in total.
#' @param show_data Logical. If \code{TRUE} (default), a subset of at most
#'   \code{ndat} data points inside the constraint region is returned.
#' @param ndat A positive integer. Maximum number of data points to return
#'   when \code{show_data = TRUE}. If more than \code{ndat} observations lie
#'   within the constraints, a random sample of size \code{ndat} is drawn
#'   without replacement. Defaults to \code{2000}.
#' @param ... Currently unused.
#'
#' @details
#' Aggregates the estimated probability content in a regular grid over the two coordinates
#' specified for plotting. For example, in three-dimensional setting, let \code{lower_constraint = upper_constraint = NULL},
#' \code{plot_coord = c(1,2)} and \code{nint=10}, then we divide the first two coordinates into 10 equally sized intervals, creating
#' a grid with 100 rectangles in \eqn{X_1} and \eqn{X_2}. Then, calculate the probability content in each rectangle by aggregating all the histogram regions that intersect
#' them. This calculates the estimated marginal distribution in \eqn{(X_1, X_2)}.
#'
#' If \code{lower_constraint} and \code{upper_constraint} are specified, then compute the
#' aggregated probability content within the 2-d rectangle in first two dimensions and further within the lower and upper constraints of the other coordinates.
#' This effectively calculates the conditional distribution of \eqn{X_1} and \eqn{X_2} within the constraints.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{\code{grid_lower}}{A numeric matrix of dimension
#'     \eqn{\texttt{nint}^2 \times d} giving the lower bounds of each grid
#'     cell in all \eqn{d} coordinates.}
#'   \item{\code{grid_upper}}{A numeric matrix of the same dimension giving
#'     the upper bounds.}
#'   \item{\code{plot_coord}}{The input coordinates for plotting.}
#'   \item{\code{est_prob}}{A numeric matrix of dimension
#'     \eqn{\texttt{nint}^2 \times 3}. The first column contains the estimated
#'     probability mass in each rectangle; the second and third columns contain the lower and upper
#'     confidence bounds of the probability mass.}
#'   \item{\code{data}}{A numeric matrix of the (at most
#'     \code{ndat}) data points within the lower and upper constraints in the
#'     \code{plot_coord}. Only present when \code{show_data = TRUE}.}
#' }
#'
#' @examples
#' \dontrun{
#' n <- 10000
#' p <- 4
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' adaptive_beta_tree <- build_adaptive_histogram(X, alpha = 0.1)
#'
#' plot_dat <- plot_histogram_data(
#'   X          = X,
#'   hist       = adaptive_beta_tree$hist,
#'   plot_coord = c(1, 2),
#'   nint       = 20,
#'   show_data  = TRUE,
#'   ndat       = 500
#' )
#' }
#'
#' @seealso
#' \code{\link{plot_histogram}} for rendering the output of this function.
#'
#' @export
plot_histogram_data <- function(X, hist,
                                lower_constraint = NULL,
                                upper_constraint = NULL,
                                plot_coord,
                                nint,
                                show_data = TRUE,
                                ndat      = 2000,
                                ...) {

  if (!is.matrix(X) || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
  if (!is.matrix(hist) || !is.numeric(hist))
    stop("'hist' must be a numeric matrix.")
  if (length(plot_coord) != 2L || any(plot_coord < 1) || any(plot_coord > ncol(X)))
    stop("'plot_coord' must be an integer vector of length 2 with valid column indices.")
  if (!is.numeric(nint) || nint < 1)
    stop("'nint' must be a positive integer.")

  d <- ncol(X)

  if (is.null(lower_constraint)) lower_constraint <- rep(-Inf, d)
  if (is.null(upper_constraint)) upper_constraint <- rep( Inf, d)

  if (length(lower_constraint) != d || length(upper_constraint) != d)
    stop("'lower_constraint' and 'upper_constraint' must each have length ncol(X).")

  # plot all data in the specified coordinates for plotting
  lower_constraint[plot_coord] <- apply(X, 2, min)[plot_coord] - 0.1
  upper_constraint[plot_coord] <- apply(X, 2, max)[plot_coord] + 0.1

  # build a grid in the plot coordinates
  xval <- seq(lower_constraint[plot_coord[1]],
              upper_constraint[plot_coord[1]],
              length.out = nint + 1)
  yval <- seq(lower_constraint[plot_coord[2]],
              upper_constraint[plot_coord[2]],
              length.out = nint + 1)

  grid_lower <- matrix(rep(lower_constraint, nint^2), ncol = d, byrow = TRUE)
  grid_lower[, plot_coord] <- as.matrix(expand.grid(xval[1:nint],
                                                    yval[1:nint]))
  grid_upper <- matrix(rep(upper_constraint, nint^2), ncol = d, byrow = TRUE)
  grid_upper[, plot_coord] <- as.matrix(expand.grid(xval[2:(nint + 1)],
                                                    yval[2:(nint + 1)]))

  # compute log-volumes of all histogram regions
  log_volume <- rowSums(log(hist[, (d + 1):(2 * d)] - hist[, 1:d]))

  # estimated probability mass in each grid cell
  est_prob <- matrix(NA_real_, nrow(grid_lower), ncol = 3)

  for (i in seq_len(nrow(grid_lower))) {

    le <- grid_lower[i, ]
    ue <- grid_upper[i, ]

    # compute log-volume of the intersection of each histogram region with grid cell i
    log_volume_intersect <- rep(NA_real_, nrow(hist))

    for (j in seq_len(nrow(hist))) {
      # check if histogram region j overlaps grid cell i in all coordinates
      overlaps <- all(hist[j, 1:d] < ue) && all(hist[j, (d + 1):(2*d)] > le)

      if (overlaps) {
        # compute volume of intersection
        intersect_lengths <- pmin(hist[j, (d + 1):(2 * d)], ue) - pmax(hist[j, 1:d], le)

        if (all(intersect_lengths > 0))
          log_volume_intersect[j] <- sum(log(intersect_lengths))
      }
    }

    # est_prob[i, k] = sum_j density_k[j] * volume_intersect[j, i]
    est_prob[i, 1] <- sum(exp(log(hist[, 2*d + 1]) + log_volume_intersect),
                          na.rm = TRUE)
    est_prob[i, 2] <- sum(exp(log(pmax(hist[, 2*d + 2], .Machine$double.eps)) +
                                log_volume_intersect), na.rm = TRUE)
    est_prob[i, 3] <- sum(exp(log(pmax(hist[, 2*d + 3], .Machine$double.eps)) +
                                log_volume_intersect), na.rm = TRUE)
  }

  if (!show_data) {
    return(list(grid_lower = grid_lower,
                grid_upper = grid_upper,
                plot_coord = plot_coord,
                est_prob   = est_prob))
  }

  # identify/sub-sample obs within the constraint region
  inside <- apply(
    (t(X[,-plot_coord]) < upper_constraint[-plot_coord]) * (t(X[,-plot_coord]) > lower_constraint[-plot_coord]), 2, prod) == 1
  index_inside <- which(inside)

  index <- if (length(index_inside) > ndat){
    sample(index_inside, ndat, replace = FALSE)
  }else{
    index_inside
  }

  list(grid_lower = grid_lower,
       grid_upper = grid_upper,
       plot_coord = plot_coord,
       est_prob   = est_prob,
       data       = X[index, plot_coord, drop = FALSE])
}


#' Plot a multivariate histogram in two dimensions
#'
#' @description
#' Visualize an adaptive beta-tree histogram in two-dimensions marginally or conditional on inside a hyperrectangle.
#'
#' @param plot_data A named list as returned by
#'   \code{\link{plot_histogram_data}}.
#' @param show_data Logical. If \code{TRUE} (default), shows the data
#'   points in \code{plot_data$data}.
#'
#' @details \code{\link{plot_histogram_data}} computes the probability content in each region,
#' \code{plot_histogram} further computes the conditional distribution after dividing by the total probability content.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @examples
#'
#' \dontrun{
#' n <- 10000
#' p <- 4
#' X <- matrix(rnorm(n*p), nrow = n, ncol = p)
#' adaptive_beta_tree <- build_adaptive_histogram(X, alpha = 0.1)
#'
#' plot_dat <- plot_histogram_data(
#'   X          = X,
#'   hist       = adaptive_beta_tree$hist,
#'   plot_coord = c(1, 2),
#'   nint       = 20,
#'   show_data  = TRUE,
#'   ndat       = 500
#' )
#'
#' plot_histogram(plot_dat)
#' }
#'
#' @seealso
#' \code{\link{plot_histogram_data}} for computing the input \code{plot_data},
#' \code{\link{build_adaptive_histogram}} for producing the histogram.
#'
#' @importFrom ggplot2 ggplot geom_rect aes scale_fill_gradient2
#'   scale_x_continuous scale_y_continuous theme_classic labs theme
#'   expansion element_blank element_text geom_point
#'
#' @export
plot_histogram <- function(plot_data, show_data = TRUE) {

  required_names <- c("grid_lower", "grid_upper", "plot_coord", "est_prob")
  missing_names  <- setdiff(required_names, names(plot_data))
  if (length(missing_names) > 0)
    stop("'plot_data' is missing required elements: ",
         paste(missing_names, collapse = ", "), ".")
  if (show_data && is.null(plot_data$data))
    stop("'plot_data$data' is NULL but show_data = TRUE. ",
         "Re-run plot_histogram_data() with show_data = TRUE.")

  # grid locations for the two plotting coordinates
  c1 <- plot_data$plot_coord[1]
  c2 <- plot_data$plot_coord[2]

  grid_lower_x <- plot_data$grid_lower[, c1]
  grid_upper_x <- plot_data$grid_upper[, c1]
  grid_lower_y <- plot_data$grid_lower[, c2]
  grid_upper_y <- plot_data$grid_upper[, c2]

  # compute conditional distribution
  grid_volume  <- (grid_upper_x - grid_lower_x) * (grid_upper_y - grid_lower_y)
  total_prob   <- sum(plot_data$est_prob[, 1])
  if (total_prob <= 0)
    warning("Total estimated probability mass is <= 0; density cannot be computed.")
  grid_density <- plot_data$est_prob[, 1] / max(total_prob, .Machine$double.eps) /
    pmax(grid_volume, .Machine$double.eps)

  mid_density <- quantile(grid_density, 0.999)

  # variable names for axis labels
  var_names <- if (!is.null(colnames(plot_data$grid_lower))){
    colnames(plot_data$grid_lower)[plot_data$plot_coord]
  }else{
    paste0("X", plot_data$plot_coord)
  }


  g <- ggplot() +
    geom_rect(aes(xmin = grid_lower_x,
                  xmax = grid_upper_x,
                  ymin = grid_lower_y,
                  ymax = grid_upper_y,
                  fill = grid_density)) +
    scale_fill_gradient2(low      = "white",
                         mid      = "#6baed6",
                         high     = "#08306b",
                         midpoint = mid_density,
                         name     = "Density") +
    scale_x_continuous(expand = expansion(0)) +
    scale_y_continuous(expand = expansion(0)) +
    theme_classic() +
    labs(x = var_names[1], y = var_names[2]) +
    theme(panel.background = element_blank(),
          panel.border      = element_blank(),
          text              = element_text(size = 18))

  if (show_data) {
    pt_alpha <- max(nrow(plot_data$data), 1L)^(-0.35)
    g <- g + geom_point(aes(x = plot_data$data[, 1],
                            y = plot_data$data[, 2]),
                        size  = 0.4,
                        alpha = pt_alpha)
  }

  return(g)
}
