#' Build a Beta Tree
#'
#' Compute a Beta Tree histogram for multivariate data at given confidence level.
#'
#' Construct the Beta Tree histogram by iteratively partitioning sample space along the sample median.
#' Then, construct *simultaneous* confidence interval (at level \eqn{\alpha}) for all of the regions using Bonferroni
#' or weighted Bonferroni adjustment. Finally, select the largest bounded regions (w.r.t. inclusion) that that pass a goodness of fit (GOF) test, i.e.,
#' the empirical density lies within the confidence interval of both a node and both of its children.
#' The function can plot the histogram if the data is two-dimensional.
#'
#' @section Multiple testing correction:
#'
#' If Bonferroni correction is selected, then the confidence level at each region is defined as
#'
#' @param X a data matrix of size n*p.
#' @param alpha confidence level, default level is \code{alpha = 0.1}.
#' @param method use \code{method = "bonferroni"} or \code{method = "weighted_bonferroni"}(default) to adjust for multiple hypothesis testing.
#' @param plot if \code{TRUE} (default), plot the histogram if the data is two-dimensional.
#' @param ... Additional parameter if \code{bounded = T}.
#' @inherit BuildKDTree
#' @returns A matrix of rectangles in the histogram of obs. in the input node.
#' A region in the histogram is represented by a row, which contains the following columns: the lower bounds, upper bounds, empirical density,
#' lower confidence bounds, upper confidence bounds of this node, number of obs. and the depth of the node.
#' @export
#' @examples
#' X <-  matrix(rnorm(2000, 0, 1), ncol = 2)
#' rect <- BuildHist(X)
#' # Plot the histogram
#' rect <- BuildHist(X, plot = TRUE)
#' # set confidence level and multiple hypothesis testing correction
#' rect <- BuildHist(X, alpha = 0.05, method = "weighted_bonferroni")
BuildHist <- function(X, alpha = 0.1, method = "weighted_bonferroni", bounded = F, plot = F, ...){
  n <- nrow(X); d <- ncol(X)
  kdtree <- BuildKDTree(X, bounded, ...)
  nd <- kdtree$nd
  ahat <- ConfLevel(nd, alpha, method, n, d)
  kdtree <- SetBounds(kdtree$kdtree, ahat, n)
  B <- matrix(nrow = 0, ncol = (2*d + 5))         # matrix that will hold selected regions
  B <- SelectNodes(kdtree, B, ahat, n)          # traverse tree and select maximal nodes that are bounded and pass GOF

  if(plot){
    if(ncol(X) == 2){
      PlotHist(X, B)
    }else{
      cat("The data dimension is not 2!")
    }
  }

  B
}
