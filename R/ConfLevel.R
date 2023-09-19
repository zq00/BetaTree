#' Set significance level for regions at each tree depth
#'
#' Computes the significance level at each tree depth.
#'
#' The root node automatically has significance level 0.
#' If \code{method = "bonferroni"}, the significance level at depth d is
#' \deqn{
#' \hat{\alpha}_d =  \frac{\alpha}{(d_\max - d_\min + 1) * n_d}, d\geq d_\min,
#' }
#' where \eqn{d_\max} is the maximum tree depth, \eqn{d_\min} is the minimum depth (greater than 0) where there are bounded regions, and \eqn{n_d} is number of bounded regions at depth \eqn{d}.
#'
#' If \code{method = weighted_bonferroni}, the confidence level at depth d is
#' \deqn{
#' \hat{\alpha}_d = \frac{\alpha}{(d_\max - d + 2) n_d \sum_{B=2}^{d_\max - d_\min + 2}1/B},d\geq d_\min.
#' }
#'
#'
#' @param nd A vector. Number of bounded regions at every depth (starting at depth 0).
#' @param alpha Numeric. The significance level (of the average density) of all the regions sums to \eqn{\alpha}.
#' @param method Use \code{method = "bonferroni"} or \code{method = "weighted_bonferroni"} method to adjust for multiple hypothesis testing.
#' @return A vector of significance level \eqn{\alpha_D} at every tree depth (starting at depth = 0)
#' @export
ConfLevel <- function(nd, alpha, method){
  DmaxP1 <- length(nd)   # tree depth + 1
  DminP1 <- max(min(which(nd > 0)), 2) # (minimum tree depth where there are bounded regions) + 1; if the root node is bounded, then set to be 2
  # set sig. level to each depth (ahat[1] corresponds to tree depth 0, and has level = 0)
  ahat <- numeric(DmaxP1)
  for(l in DminP1:DmaxP1){
    switch (method,
            "bonferroni" = {ahat[l] <- alpha / (DmaxP1 - DminP1 + 1) / nd[l]},
            "weighted_bonferroni" = {
              wsum <- sum(1/(2:(DmaxP1-DminP1 + 2)))
              ahat[l] <- alpha / (DmaxP1 - l + 2) / nd[l] / wsum
            }
    )
  }

  return(ahat)
}

