#' Add confidence bounds
#'
#' Recursively set confidence bounds at bounded nodes in the k-d tree.
#'
#' The upper and lower confidence bounds of a bounded region at depth \eqn{D} is given by
#' \deqn{
#' \text{lower} = \frac{qBeta(\hat{\alpha}_D/2, N+1, N-n)}{|R|}\\
#' \text{upper} = \frac{qBeta(1-\hat{\alpha}_D/2, N+1, N-n)}{|R|}
#' }
#' where \eqn{\hat{\alpha}_D} is the significance level at depth \eqn{D}, \eqn{N} is the total number of obs.,  \eqn{n} is the number of obs. in this region,
#' and |R| is the volume of the region. If there are no children nodes, ignore the lower and upper bounds of the children.
#' We set the lower and upper bound of each region as
#' \deqn{
#' \widetilde{\text{lower}} = \max\left(\text{lower}, \text{lower}(\text{2nd child})\right) \\
#' \widetilde{\text{upper}} = \max\left(\text{upper}, \text{upper}(\text{1st child}), \text{upper}(\text{2nd child})\right),
#' }
#' so that when we select nodes, we pick the largest nodes such that the empirical density is between \eqn{\widetilde{\text{lower}}} and  \eqn{\widetilde{\text{upper}}}.
#' The confidence bounds are computed in the function \code{ConfLevel}.
#' @param node A list element representing the region.
#' @param ahat A vector of confidence levels at each depth, starting at depth = 1. Output of the funtion \code{ConfLevel}
#' @param n Total number of observations.
#' @return An updated k-d tree starting at input node.
#' @export
SetBounds <- function(node, ahat, n){
  if (node$bounded){
    lower <- stats::qbeta(ahat[node$depth + 1] / 2, node$ndat + 1, n - node$ndat) / prod(node$up - node$low)
    upper <- stats::qbeta(1 - ahat[node$depth + 1] / 2, node$ndat + 1, n - node$ndat) / prod(node$up - node$low)
  }
  if(node$leaf){if (node$bounded) {node$lower <- lower; node$upper <- upper}
  }else{
    node$leftchild <- SetBounds(node$leftchild, ahat, n)
    node$rightchild <- SetBounds(node$rightchild, ahat, n)
  }
  if(!node$leaf & node$bounded){
    node$lower <- max(lower, node$leftchild$lower,node$rightchild$lower)
    node$upper <- min(upper, node$leftchild$upper, node$rightchild$upper)
  }
  return(node)
}
