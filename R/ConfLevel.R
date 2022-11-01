#' Set confidence level to each depth
#'
#' Computes the confidence level at every level of the tree.
#'
#' If \code{method = Bonferroni}, the confidence level at depth d is
#' \deqn{
#' \hat{\alpha}_d = \frac{\alpha}{(d_\max - d_\min + 1) * n_d}, D\geq D_\min,
#' }
#' where \eqn{d_\max} is the maximum depth of the tree, \eqn{d_\min} is the smallest depth where there are bounded regions, and \eqn{n_d} is number of bounded regions at depth \eqn{d}.
#'
#' If \code{method = weighted_bonferroni}, the confidence level at depth d is
#' \deqn{
#' \hat{\alpha}_d = \frac{\alpha}{(d_\max - d + 2) n_d \sum_{B=2}^{d_\max - d_\min + 2}1/B},D\geq D_\min,
#' }
#' that is, the leaves get weights alpha/2 etc.
#'
#' @param nd A vector. Number of bounded regions at every level.
#' @param alpha Numeric. confidence level.
#' @param method Use \code{bonferroni} or \code{weighted_bonferroni} method to adjust for multiple hypothesis testing.
#' @return A vector of confidence level at every depth (starting at depth = 1)
#' @export
ConfLevel <- function(nd, alpha, method){
  Dmax <- length(nd)  # depth of tree
  Dmin <- min(which(!is.na(nd)))
  # set sig. level to each depth
  ahat <- numeric(Dmax)
  for(l in Dmin:Dmax){
    switch (method,
            "bonferroni" = {ahat[l] <- alpha / (Dmax - Dmin + 1) / nd[l]},
            "weighted_bonferroni" = {
              wsum <- sum(1/(2:(Dmax-Dmin + 2)))
              ahat[l] <- alpha / (Dmax - l + 2) / nd[l] / wsum
            }
    )
  }
  if(Dmin>1) {  # dummy value: first few will not be used as there are no bounded regions
    ahat[1:(Dmin-1)] <- 0.5
  }

  return(ahat)
}

