#' Identify Modes in a Beta-tree Histogram via Monte Carlo approximation
#'
#' @description
#' Implements an approximate mode-hunting algorithm by using randomly sampled paths
#' to estimate whether a region is a distinct mode.
#'
#' @param hist A numeric matrix of dimension \eqn{m\times (2d+5)}.
#'    Each row represents one region in the histogram.
#'    If the histogram is constructed for \code{d} dimensional data,
#'    then the \eqn{1:d} columns contain the lower bounds for each coordinate,
#'    columns \eqn{(d+1):(2d)} contain the upper bounds.
#'    The \eqn{2d + 1} column contains the empirical density.
#'    The \eqn{(2d+2):(2d+3)} columns contain the lower and upper confidence bounds of the empirical density.
#'    \code{hist} can be the output of \link[BetaTree]{BuildHist} or \link[BetaTree]{build_adaptive_histogram} function.
#' @param d A positive integer corresponding to the data dimension.
#' @param L A positive integer. The length of random path.
#' @param B A positive integer. The number of independent random paths
#'   used in \code{\link{connected_to_modes}} to assess whether a candidate region
#'   is connected to any of the current modes.
#' @param verbose Prints progress  if \code{TRUE}.
#'
#' @details
#' Evaluate every region in order of decreasing density whether it is a distinct mode compared to current modes.
#' In the exact algorithm \link{FindModes}, a region \eqn{R} is considered a distinct mode from a current model \eqn{M} if
#' if along every path connecting \eqn{R} and \eqn{M}, there exists at least one region whose upper confidence bound of density
#' is below the lower confidence bounds of both \eqn{R} and \eqn{M}. The approximate algorithm here evaluates \code{B} random paths of length \code{L}
#' instead of evaluating every path.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{\code{mode}}{An integer vector of row indices of \code{hist}
#'     identifying the detected modes.}
#'   \item{\code{hist}}{The input \code{hist} matrix.}
#'   \item{\code{cluster}}{A list of length \code{length(mode)}. The
#'     \eqn{k}-th element is an integer vector of row indices of histogram regions
#'     that belonging to the neighborhood of mode \eqn{k}. That is, there is a path from this region to the mode
#'     that the upper confidence bound of every region along the path (including the region itself) is higher than the lower confidence bound of the mode. }
#' }
#'
#' @examples
#' \dontrun{
#' # mixture of two Gaussian distribution
#' X <- rbind(matrix(rnorm(500, mean = -2), ncol = 2),
#'            matrix(rnorm(500, mean =  2), ncol = 2))
#' beta_tree_hist <- BuildHist(X)
#'
#' # find modes using approximate algorithm
#' modes <- FindModesApproximate(hist  = beta_tree_hist,
#'                               d     = 2,
#'                               L     = 50,
#'                               B     = 1000,
#'                               verbose = F)
#'
#' # number of detected modes
#' length(modes$mode)
#' }
#'
#' @export
FindModesApproximate <- function(hist, d, L, B, verbose = FALSE) {

  if (!is.matrix(hist) || !is.numeric(hist))
    stop("'hist' must be a numeric matrix.")
  if (ncol(hist) != 2 * d + 5)
    stop("'hist' must have 2*d + 5 columns.")
  if (!is.numeric(L) || L < 1)
    stop("'L' must be a positive integer (random walk length).")
  if (!is.numeric(B) || B < 1)
    stop("'B' must be a positive integer (number of random walks).")


  if (nrow(hist) == 1L) {
    cat("The histogram contains a single region!\n")
    return(list(mode    = 1L,
                hist    = hist,
                cluster = list(1L)))
  }

  cluster <- list()  # stores regions in the neighbourhood of mode k

  # extract empirical density and confidence bounds of each region
  density    <- hist[, 2 * d + 1]           # empirical density estimate
  ci         <- hist[, (2*d + 2):(2*d + 3)] # lower and upper CI bounds

  node_order <- order(density, decreasing = TRUE)

  # compute adjacency graph from adjacency matrix
  adj <- compute_adjacency_mat(hist, d)
  g   <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  rm(adj)

  # first mode is the region with the highest density
  mode       <- node_order[1]
  cluster[[1]] <- find_neighbors(mode[1], g, ci)
  candidates <- setdiff(node_order[-1], cluster[[1]])

  if (verbose) {
    cat("# neighbors of mode 1 = ", length(cluster[[1]]), "\n")
    cat("# candidate modes      = ", length(candidates),  "\n")
  }

  nmodes <- 1L

  while (length(candidates) > 0) { # evaluate every candidate mode

    connected <- connected_to_modes(candidates[1], mode, g, cluster, ci, L, B)

    if (connected$connected == "unconnected") {

      nmodes       <- nmodes + 1L
      mode         <- c(mode, candidates[1])
      cluster[[nmodes]] <- find_neighbors(mode[nmodes], g, ci)
      candidates   <- setdiff(candidates[-1], cluster[[nmodes]])

      if (verbose) {
        cat("# neighbors of mode ", nmodes, "=", length(cluster[[nmodes]]), "\n")
        cat("# candidate modes  =", length(candidates), "\n")
      }

    } else {
      candidates <- candidates[-1]
    }
  }

  return(list(mode    = mode,
              hist    = hist,
              cluster = cluster))
}


#' Find Regions in the Neighborhood of a Mode
#'
#' @description
#' Computes the set of histogram regions that form the neighborhood of a
#' given mode.
#'
#' @param mode A single positive integer giving the row index of the mode
#'   region in the histogram matrix.
#' @param g An \code{igraph} undirected graph object representing the
#'   adjacency structure of the histogram regions, as produced by
#'   \code{\link[igraph]{graph_from_adjacency_matrix}} from the output of
#'   \code{\link{compute_adjacency_mat}}.
#' @param ci A numeric matrix with two columns. Row \eqn{i} gives the
#'   lower (column 1) and upper (column 2) confidence bounds for the
#'   average density of histogram region \eqn{i}.
#'
#' @details
#' A region \eqn{R} is considered to be in the neighborhood of a model \eqn{M} if
#' there exists a path connecting \eqn{R} and \eqn{M} such that the upper confidence
#' bound of every region along the path is above the lower confidence bound of \eqn{M}.
#'
#' @return An integer vector of row indices (including \code{mode}) of
#'   the regions belonging to the neighbourhood of \code{mode}.
#'
#' @seealso
#' \code{\link{FindModesApproximate}} for the function that calls
#' \code{find_neighbors}.
#'
#' @export
find_neighbors <- function(mode, g, ci) {

  if (!igraph::is.igraph(g))
    stop("'g' must be an igraph graph object.")
  if (!is.matrix(ci) || ncol(ci) != 2)
    stop("'ci' must be a two-column numeric matrix.")
  if (mode < 1 || mode > nrow(ci))
    stop("'mode' must be a valid row index of 'ci'.")

  # initialize with direct neighbors of the mode region whose upper CI > the mode's lower CI
  neighbor_nodes <- as.integer(igraph::neighbors(g, mode))
  x <- neighbor_nodes
  x <- x[ci[x, 2] > ci[mode, 1]]

  # add x's neighbors whose upper CI > the mode's lower CI
  while (length(x) > 0) {

    new_nbd <- unique(unlist(lapply(x, function(t) as.integer(igraph::neighbors(g, t)))))
    x <- new_nbd[ci[new_nbd, 2] > ci[mode, 1]]

    x <- setdiff(x, neighbor_nodes)
    neighbor_nodes <- unique(c(neighbor_nodes, new_nbd))
  }

  # append the mode region itself
  return(c(neighbor_nodes, mode))
}


#' Is a Region a Distinct Mode? (Approximate Algorithm)
#'
#' @description
#' Tests whether a histogram region is connected to any of the
#' current modes using Monte Carlo algorithm.
#'
#' @param candidate_mode A positive integer giving the row index of
#'   the candidate region.
#' @param curr_mode An integer vector of row indices identifying the
#'   currently identified mode regions in the histogram.
#' @param g An \code{igraph} undirected graph object representing the
#'   adjacency graph of the histogram regions.
#' @param cluster A list of length \code{length(curr_mode)}, where the
#'   \eqn{k}-th element is a vector indicating regions in the
#'   neighborhood of mode \eqn{k}, as returned by
#'   \code{\link{find_neighbors}}.
#' @param ci A numeric matrix with two columns giving the lower
#'   and upper confidence bounds for the average density of each histogram region.
#' @param L A positive integer. The length of random path.
#' @param B A positive integer. The number of independent random paths
#'   used to assess whether a candidate region is connected to any of the current modes.
#'
#' @details Generate \code{B} random paths of length \code{L} starting from
#' \code{candidate_mode} and evaluate (1) whether it connects the \code{candidate_mode}
#' to any \code{curr_mode}; if so, (2) whether there exists a region along the path whose upper
#' confidence bound is below the lower confidence bounds of both the \code{candidate_mode} and
#' the \code{curr_mode}. If such a region does not exist, then the \code{candidate_mode} *cannot* be a distinct mode.
#'
#' @return A named list with one element:
#' \describe{
#'   \item{\code{connected}}{A character string: \code{"connected"} if
#'     the candidate is *not* a distinct mode and \code{"unconnected"} otherwise.}
#' }
#'
#' @seealso
#' \code{\link{FindModesApproximate}} calls this function to test whether a candidate region is a distince mode.
#'
#' @export
connected_to_modes <- function(candidate_mode, curr_mode, g, cluster,
                               ci, L, B) {

  if (!igraph::is.igraph(g))
    stop("'g' must be an igraph graph object.")
  if (!is.matrix(ci) || ncol(ci) != 2)
    stop("'ci' must be a two-column numeric matrix.")
  if (L < 1)
    stop("'L' must be a positive integer.")
  if (B < 1)
    stop("'B' must be a positive integer.")

  connected <- FALSE
  nmode     <- length(curr_mode)

  for (b in seq_len(B)) {

    # generate a random path starting from the candidate_mode
    path <- as.integer(igraph::random_walk(g, candidate_mode, steps = L))

    for (j in seq_len(nmode)) {
      intersections <- intersect(path, cluster[[j]]) # does this path intersect neighbors of j-th current mode?

      if (length(intersections) == 0L) { # if not, check next mode
        next;
      } else {
        ind <- min(which(path %in% cluster[[j]])) # first region of intersection

        val <- vapply(path[seq_len(ind)],
                      function(t)
                        ci[t, 2] < ci[curr_mode[j], 1] &&
                        ci[t, 2] < ci[candidate_mode, 1],
                      logical(1))

        if (all(!val)) {
          # if no region with upper CI < lower CI of both current mode and candidate mode
          # then candidate mode cannot be a distinct mode
          connected <- TRUE
          break
        }
      }
    }

    if (connected) break
  }

  return(list(
    connected = if (connected) "connected" else "unconnected"
  ))
}
