#' Are two regions distinct modes?
#'
#' Two regions \code{i} and \code{j} are two *distinct* modes if
#' there *exists* one path connecting \code{i} and \code{j}, there exists at least one \code{R} along the path,
#' such that the upper endpoint of confidence interval of \code{R} is *below* the lower endpoints of that of
#' both \code{i} and \code{j}.
#'
#' For computational reasons, we only check paths whose length is at most \code{cutoff}.
#'
#' @param i,j Test if the \eqn{i}-th and the \eqn{j}-th regions are distinct modes or not.
#' @param g A graph returned by \link[igraph]{graph_from_adjacency_matrix}
#' @param ci A matrix of size N*2 of the confidence intervals of empirical densities.
#'   Each row represents one region in the Beta Tree and the two columns correspond to the lower and upper confidence limits.
#' @inherit FindMode
#' @returns \code{unconnected} if the two regions are distinct modes and \code{connected} otherwise.
#' @export

is_connected <- function(i, j, g, ci, cutoff = 6){
  all_path <- igraph::all_simple_paths(g, i, j, mode = "all", cutoff = cutoff)
  if(length(all_path) == 0) return("unconnected") # no path connecting the two regions

  for(k in 1:length(all_path)){
    path <- all_path[[k]]
    if(length(path) == 2){ # if only i and j are in the path, then they are neighbors
      return("connected")
    }
    val <- sapply(path, function(t) ci[t,2] < ci[i,1] & ci[t,2] < ci[j,1]) # T if the CI of a rectangle along the path dips below the CI of both i AND j
    if(all(!val)){
      return("connected")
    }
  }
  return("unconnected")
}
