#'Select nodes in the histogram
#'
#'Extract maximal (w.r.t. inclusion, i.e. tree order) nodes that are
#'bounded and whose empirical density is within confidence bounds.
#'
#'A bounded node is selected if the empirical density is within the confidence bounds.
#'The empirical density is
#'\deqn{
#'h = \frac{n+1}{N |R|},
#'}
#'where \eqn{N} is the number of obs., and \eqn{n} is the number of obs. in the region,
#' and \eqn{|R|} is the volume of the region.
#'Otherwise, if the node is not a leaf, move on to the left and right children of the node.
#'
#'@param node A k-d tree with confidence thresholds, output from \code{SetBounds}.
#'@param B A matrix.
#'@param ahat A vector of confidence levels at each depth, starting at depth = 1. Output of the funtion \code{ConfLevel}
#'@param n Total number of observations.
#'@return A matrix of rectangles in the histogram of obs. in the input node.
#'    A region in the histogram is represented by a row, which contains the lower and upper bounds, empirical density,
#'    lower and upper confidence bounds of this node, number of obs. and the depth of the node.
SelectNodes <- function(node, B, ahat, n){
  recurse <- TRUE
  if(node$bounded){
    h <- (node$ndat +1)/n/prod(node$up - node$low)
    # cat(node$depth, ",", h, ",", node$lower, ",", node$upper, "\n")
    if ((node$lower <= h) & (h <= node$upper)) { # take this node
      recurse <- FALSE
      lower <- stats::qbeta(ahat[node$depth] / 2, node$ndat + 1, n - node$ndat) / prod(node$up - node$low)
      upper <- stats::qbeta(1 - ahat[node$depth] / 2, node$ndat + 1, n - node$ndat) / prod(node$up - node$low)
      B <- rbind(B,c(node$low, node$up,h, lower,upper, node$ndat ,node$depth))
    }
  }
  if (recurse & !node$leaf) {B <- SelectNodes(node$leftchild, B, ahat, n); B <- SelectNodes(node$rightchild, B,ahat, n)}
  return(B)
}
