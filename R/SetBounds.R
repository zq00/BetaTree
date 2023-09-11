#' Add confidence bounds to bounded regions in the k-d tree
#'
#' Recursively set confidence bounds for bounded nodes in the k-d tree.
#'
#' Because
#' \deqn{
#' F(R_k)\sim \mathrm{Binomial}(n_k+1,n-n_k),
#' }
#' where \eqn{R_k} is a region, \eqn{n} is the total number of obs. and \eqn{n_k} is the number of obs. in \eqn{R_k},
#' the upper and lower confidence bounds of the average density in a bounded region \eqn{R_k} at depth \eqn{D} is given by
#' \deqn{
#' \text{lower}(R_k) = \frac{qBeta(\hat{\alpha}_D/2, n_k+1, n-n_k)}{|R_k|}\\
#' \text{upper}(R_k) = \frac{qBeta(1-\hat{\alpha}_D/2, n_k+1, n-n_k)}{|R_k|}
#' }
#' where \eqn{\hat{\alpha}_D} is the significance levels at depth \eqn{D},
#' and \eqn{|R_k|} is the volume of the region.
#'
#' We select nodes in which data are approximately uniformly distributed, i.e., for which the empirical density falls within the confidence intervals of all of its children.
#' We set the lower and upper bound of each region as
#' \deqn{
#' \widetilde{\text{lower}} = \max\left(\text{lower}, \text{lower}(\text{1st child}), \text{lower}(\text{2nd child})\right) \\
#' \widetilde{\text{upper}} = \min\left(\text{upper}, \text{upper}(\text{1st child}), \text{upper}(\text{2nd child})\right),
#' }
#' so that when we select nodes, we pick the largest nodes such that the empirical density is between \eqn{\widetilde{\text{lower}}} and  \eqn{\widetilde{\text{upper}}}.
#' The significance levels are computed with the function \link[BetaTrees]{ConfLevel}.
#'
#' @param node A list storing information about a region.
#' @param ahat A vector of significance level at each depth, starting at depth = 0. Output of the funtion \code{ConfLevel}.
#' @param n Total number of observations.
#' @return An updated k-d tree starting at input node with \code{lower} and \code{upper} filled in.
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
