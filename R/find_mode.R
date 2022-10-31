#' Compute modes in the Beta Tree histogram
#'
#' Identifies distinct modes in the histogram
#'
#' Two regions \code{i} and \code{j} are distinct modes if along at least one path connecting them,
#'   there exists one region whose upper confidence limit is below the lower confidence limit
#'   of both \code{i} and \code{j}.
#' This function iterates through all of the regions in the Beta Tree histogram according to the
#'   descending order of empirical density and marks a region as a distinct mode if it is a distinct mode
#'   compared to every current modes.
#'
#' @param hist A matrix. Each row represents one region in the histogram.
#'   If the histogram is constructed for \code{d} dimensional data, then the first \code{d} columns contain the lower bound of each dimension,
#'   and the next \code{d} columns contain the upper bound. The \code{2d + 1} column contains the empirical density.
#'    \code{hist} can be the output of \link[BetaTrees]{BuildHist} function.
#' @param d A numeric value of the data dimension.
#' @param cutoff Maximum path length. By default, we only evaluate paths with length at most \code{cutoff = 6}.
#' @returns
#' \describe{
#' \item{mode}{A vector of candidate modes.}
#' \item{hist}{A matrix of the input Beta Tree histogram.}
#' \item{g}{A graph object constructed using the adjacency matrix of the regions in the histogram.}
#' }
#' @export
#' @examples
#' # A mixture of two Gaussians
#' X <-  matrix(rnorm(4000, 0, 1), ncol = 2)
#' X[1:1000, 1] <- X[1:1000, 1] + 4
#' rect <- BuildHist(X) # construct a Beta Tree histogram
#' modes <- FindModes(hist = rect, d = 2)$mode

FindModes <- function(hist, d, cutoff = 6){
  density <- hist[,(2*d+1)] # empirical density of each region
  ci <- hist[, (2*d+2):(2*d+3)] # lower and upper confidence bounds of each region
  node_order <- order(density, decreasing = T) # order nodes by decreasing order of empirical density
  candidate_mode <- node_order[1] #

  # compute adjacency matrix
  adj <- compute_adjacency_mat(hist, d)
  # compute a graph from the adj
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  for(i in 2:nrow(hist)){ # iterate through every region
    flag  <-  F # If T, then the region i is not a mode
    for(j in candidate_mode){ # iterate through every current mode
      if(adj[node_order[i], j] == 1) {flag = T; break;} # connected to current mode
      new_val <- is_connected(node_order[i], j, g, ci, cutoff)
      if(new_val == "connected") {flag = T; break; }
    }
    if(!flag) candidate_mode <- c(candidate_mode, node_order[i])
  }
  return(list(mode = candidate_mode,
              hist = hist,
              g = g))

}