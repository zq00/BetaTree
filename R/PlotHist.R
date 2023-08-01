#' Plot a 2D Beta tree histogram
#'
#' Draw a Beta tree histogram.
#'
#' \code{PlotHist} plots the histogram superimposed on the original dataset.
#' If the number of data points exceeds \code{nsample}, randomly select \code{nsample} obs. to show in the graph. Default value is \eqn{n_sample=10^4}
#'
#' @param X data matrix of size n by 2.
#' @param B matrix of rectangles returned by the function \link[BetaTrees]{SelectNodes} or \link[BetaTrees]{BuildHist}.
#' @param nsample maximum number of data points to display in the plot.
#' @param ... additional parameter for the graph: \code{low}, \code{mid}, \code{high} specifies the colors in \link[ggplot2]{scale_fill_gradient2} function.
#' @return A \code{ggplot} object.
#' @importFrom stats median
#' @importFrom ggplot2 ggplot
#' @export
#' @examples
#' X <-  matrix(rnorm(2000, 0, 1), ncol = 2)
#' rect <- BuildHist(X, plot = FALSE)
#' PlotHist(X, rect)
PlotHist <- function(X, B, nsample = 10000, ...){
  n <- nrow(X)
  boundary <- data.frame(
    xmin = B[,1], xmax = B[,3], ymin = B[,2],ymax = B[,4],
    nobs = B[,8], density = B[,5]
  )
  # take a subsample if n is large
  if(n > nsample){ x_small <- X[sample(1:n,  nsample, replace = FALSE), ]
  }else{ x_small <- X}

  # User can input the two colors in the color gradients of the histogram
  params <-  list(...)
  if(is.null(params$low)){low <- "#f7fbff"}
  if(is.null(params$mid)){ mid <- "#6baed6"}
  if(is.null(params$high)){ high <- "#08306b"}
  g <- ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = boundary$xmin, xmax = boundary$xmax,ymin = boundary$ymin, ymax = boundary$ymax, fill = boundary$density), color = "black")  +
    ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high, midpoint = stats::median(boundary$density), aesthetics = "fill") +
    ggplot2::geom_point(ggplot2::aes(x = x_small[,1], y = x_small[,2]), size = 0.4, alpha =n^(-0.25)) +
    ggplot2::xlab("x") +  ggplot2::ylab("y") +
    ggplot2::theme(text=ggplot2::element_text(size =18))

  return(g)
}

