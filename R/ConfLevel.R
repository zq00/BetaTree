#' Set confidence level for regions at each tree depth
#'
#' Computes the confidence level at each tree depth.
#'
#' If \code{method = "bonferroni"}, the significance level at depth d is
#' \deqn{
#' \hat{\alpha}_d =  \frac{\alpha}{(d_\max - d_\min + 1) * n_d}, D\geq D_\min,
#' }
#' where \eqn{d_\max} is the maximum tree depth, \eqn{d_\min} is the minimum depth where there are bounded regions, and \eqn{n_d} is number of bounded regions at depth \eqn{d}.
#'
#' If \code{method = weighted_bonferroni}, the confidence level at depth d is
#' \deqn{
#' \hat{\alpha}_d = \frac{\alpha}{(d_\max - d + 2) n_d \sum_{B=2}^{d_\max - d_\min + 2}1/B},D\geq D_\min.
#' }
#'
#'
#' @param nd A vector. Number of bounded regions at every depth (starting at depth 0).
#' @param alpha Numeric. The significance level (of the average density) of all the regions adds to \eqn{\alpha}.
#' @param method Use \code{method = "bonferroni"} or \code{method = "weighted_bonferroni"} method to adjust for multiple hypothesis testing.
#' @return A vector of significance level \eqn{\alpha_D} at every tree depth (starting at depth = 0)
#' @export
ConfLevel <- function(nd, alpha, method){
  Dmax <- length(nd)  # depth of tree
  Dmin <- min(which(nd > 0))
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
    ahat[1:(Dmin-1)] <- 0
  }

  return(ahat)
}

