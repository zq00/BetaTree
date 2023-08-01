#' Build a Beta tree histogram
#'
#' Compute a Beta tree histogram for multivariate data at a given confidence level.
#'
#' Construct the Beta tree histogram by (1) iteratively partitioning sample space along the sample median (similar to a k-d tree).
#' (2) construct *simultaneous* confidence interval (at level \eqn{\alpha}) for all of the regions (adjust for multiple testing with Bonferroni
#' or weighted Bonferroni adjustment). (3) select the largest bounded regions (w.r.t. inclusion) that that pass a goodness of fit (GOF) test, i.e.,
#' the empirical density lies within the confidence interval of all of the further partitions of the region (i.e., descendants of the node in the k-d tree).
#'
#' @param X A data matrix of size n by d.
#' @param alpha Significance level, default is \code{alpha = 0.1}. The CI covers the average density in every region of the Beta tree histogram simultaneously with probability \eqn{1-\alpha}.
#' @param method Use \code{method = "bonferroni"} or \code{method = "weighted_bonferroni"}(default) to adjust for multiple hypothesis testing.
#' @param plot if \code{TRUE} (default), plot the histogram if the data is two-dimensional.
#' @param ... Additional parameter if \code{bounded = T}.
#' @inherit BuildKDTree
#' @returns A matrix. Each row represents one regions in the histogram.
#' The first \eqn{d} columns are the lower bounds of the region; the next \eqn{d} columns are the upper bounds (\eqn{d} is data dimension);
#' the \eqn{2d+1} column stores the empirical density in the region; the next two columns  are
#' the lower and upper confidence bounds of the average density in the region; the last two columns are the number of obs. inside the region and the depth of the node in the k-d tree.
#'
#' @export
#' @examples
#' X <-  matrix(rnorm(2000, 0, 1), ncol = 2)
#' B <- BuildHist(X)
#' # Plot the histogram
#' B <- BuildHist(X, plot = TRUE)
#' # Set significance level and multiple hypothesis testing correction
#' B <- BuildHist(X, alpha = 0.05, method = "weighted_bonferroni")
#' # Initialize a bounding box
#' B <- BuildHist(X, alpha = 0.05, method = "weighted_bonferroni",
#'      bounded = TRUE, option = "ndat", qt = c(1,1))
BuildHist <- function(X, alpha = 0.1, method = "weighted_bonferroni", bounded = F, plot = F, ...){
  n <- nrow(X); d <- ncol(X)
  kdtree <- BuildKDTree(X, bounded, ...)
  nd <- kdtree$nd
  ahat <- ConfLevel(nd, alpha, method)
  kdtree <- SetBounds(kdtree$kdtree, ahat, n)
  B <- matrix(nrow = 0, ncol = (2*d + 5))         # matrix that will hold selected regions
  B <- SelectNodes(kdtree, B, ahat, n)          # traverse tree and select maximal nodes that are bounded and pass GOF

  if(plot == T){
    if(ncol(X) == 2){
      g <- PlotHist(X, B, ...)
      print(g)
      return(B)
    }else{
      cat("The data dimension is not 2!")
      return(B)
    }
  }else{
    return(B)
  }
}
