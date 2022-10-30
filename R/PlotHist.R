#' Plot a 2D Beta Tree histogram
#'
#' Draw a Beta tree histogram from the rectangles returned by \code{SelectNodes}.
#'
#' \code{PlotHist} plots the histogram superimposed on the original dataset.
#' If the number of data points exceeds \code{nsample}, randomly select \code{nsample} obvervations to show in the graph. Default value is \eqn{n_sample=10^4}
#'
#' @param X data matrix of size n by 2.
#' @param B matrix of rectangles returned by the function \link[BetaTrees]{SelectNodes} or \link[BetaTrees]{BuildHist}.
#' @param nsample maximum number of data points to display in the plot.
#'
#' @return A \code{ggplot} object.
#' @export
#' @examples
#' X <-  matrix(rnorm(2000, 0, 1), ncol = 2)
#' rect <- BuildHist(X, plot = F)
#' PlotHist(X, rect)
PlotHist <- function(X, B, nsample = 10000){
  n <- nrow(X)
  boundary <- data.frame(
    xmin = B[,1], xmax = B[,3], ymin = B[,2],ymax = B[,4],n = B[,8],
    density = B[,5]
  )
  # take a subsample if N is large
  if(n > nsample){ x_small <- X[sample(1:n,  nsample, replace = FALSE),]
  }else{ x_small <- X}
  g <- ggplot2::ggplot() +
    ggplot2::geom_rect(ggplot2::aes(xmin = xmin, xmax = xmax,ymin = ymin, ymax = ymax,fill = density), color = "black",data = boundary)  +
    ggplot2::scale_fill_gradient2(low = "#e0f3db", mid = "#a8ddb5", high = "#43a2ca",aesthetics = "fill") +
    ggplot2::geom_point(ggplot2::aes(x = x_small[,1], y = x_small[,2]), size = 0.4, alpha =n^(-0.25)) +
    ggplot2::xlab("X") +  ggplot2::ylab("Y") +
    ggplot2::theme_bw() +  ggplot2::theme(text=ggplot2::element_text(size =18))

  return(g)
}

