#' Adjacency matrix of regions in the Beta Tree
#'
#' Computes the adjacency matrix of all of the regions in a Beta Tree histogram
#'
#' Two regions are neighbors if and only if the interval (defined by the boundary) in each dimension intersects.
#' For example, let the two regions \eqn{X} defined by \eqn{[X_{i,\mathrm{low}}], X_{i,\mathrm{up}}]}
#' and \eqn{Y} defined by \eqn{[Y_{i,\mathrm{low}}], Y_{i,\mathrm{up}}]} (for \eqn{i=1,\ldots, d}), where \eqn{[X_{i,\mathrm{low}}]} and \eqn{X_{i,\mathrm{up}}]}
#' are the lower and upper bounds in dimension \eqn{i} for the region \eqn{X}. Then, \eqn{X} and \eqn{Y} are neighbors if and only if
#' \eqn{[X_{i,\mathrm{low}}], X_{i,\mathrm{up}}] \cap [Y_{i,\mathrm{low}}], Y_{i,\mathrm{up}}] \neq \emptyset} for every \eqn{i=1,\ldots, d}.
#'
#' @inherit  FindModes
#' @returns A matrix of size N*N, where N is the number of regions in the histogram, i.e. number of rows in the histogram.
#' @export
#' @examples
#' X <-  matrix(rnorm(2000, 0, 1), ncol = 2)
#' rect <- BuildHist(X)
#' adj <- compute_adjacency_mat(hist = rect, d = 2)

compute_adjacency_mat <- function(hist, d){
  N <- nrow(hist)
  adj <- matrix(0, N, N)
  # Find neighbors of each rectangle in the histogram (by removing rectangles NOT in the nbd)
  nbd <- list()
  for(i in 1:N){
    nbd[[i]] <- (1:N)[-i]# initiate neighbors by every rectangle that is NOT itself
    for(dim in 1:d){ # iterate through each dimension
      if(sum(nbd[[i]]) == 0) break; # region i does not have any neighors
      x <- hist[i, c(dim, d + dim)]; # lower and upper bounds of the region
      y <- rect[nbd[[i]], c(dim, d + dim), drop = F] # lower and upper bounds of candidate neighbors
      indices <- which(y[,2] < x[1] | y[,1] > x[2]) # these are NOT in the nbd
      nbd[[i]] <- nbd[[i]][-indices]
    }
    if(sum(nbd[[i]]) > 0) {
      adj[i, nbd[[i]]] <- 1 # fill in adjacency matrix
      }
  }

  return(adj)
}
