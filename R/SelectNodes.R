#'Select nodes in the histogram
#'
#'Extract maximal (w.r.t. inclusion) nodes that are
#'bounded and whose empirical density is within confidence bounds.
#'
#'A bounded node is selected if the empirical density is inside the confidence bounds of all of its children.
#'The empirical density is defined as
#'\deqn{
#'h = \frac{n+1}{N |R|},
#'}
#'where \eqn{N} is the total number of obs., \eqn{n} is the number of obs. in the region,
#' and \eqn{|R|} is the volume of the region.
#'If the node is not selected and is not a leaf, then move on to the left and right children of the node.
#'
#'@param node A k-d tree output from \code{SetBounds}.
#'@param B A matrix containing nodes that have been selected.
#'@param ahat A vector of significance levels at each depth, starting at depth = 0. Output of the funtion \code{ConfLevel}
#'@param n Total number of observations.
#'@return A matrix of selected regions within and including the input node.
#'    Each row represents one selected region and has \eqn{2d + 5} columns, containing the lower and upper bounds of the region,
#'    the empirical density, lower and upper confidence bounds of the average density in this region,
#'    number of obs. inside and the depth of the node.
#' @importFrom stats qbeta
#' @export
SelectNodes <- function(node, B, ahat, n){
  recurse <- TRUE
  if(node$bounded){
    h <- (node$ndat +1)/n/prod(node$up - node$low)
    # cat(node$depth, ",", h, ",", node$lower, ",", node$upper, "\n")
    if ((node$lower <= h) & (h <= node$upper)) { # take this node
      recurse <- FALSE
      lower <- stats::qbeta(ahat[node$depth + 1] / 2, node$ndat + 1, n - node$ndat) / prod(node$up - node$low)
      upper <- stats::qbeta(1 - ahat[node$depth + 1] / 2, node$ndat + 1, n - node$ndat) / prod(node$up - node$low)
      B <- rbind(B,c(node$low, node$up,h, lower,upper, node$ndat ,node$depth))
    }
  }
  if (recurse & !node$leaf) {B <- SelectNodes(node$leftchild, B, ahat, n); B <- SelectNodes(node$rightchild, B,ahat, n)}
  return(B)
}
