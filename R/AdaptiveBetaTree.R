#' Choose partition dimension
#'
#' @description
#' Selects partition dimension by testing marginal uniformity and pairwise independence.
#'
#' @param x A numeric matrix dimension \eqn{n \times d} containing observations falling inside the current
#'   region. Each row corresponds to one observation.
#' @param thresh_marginal A numeric value between \eqn{(0, 1)}. The significance
#'   threshold for the marginal Anderson-Darling uniformity test.
#' @param thresh_interaction A positive numeric value. The threshold on the
#'   combined interaction score \eqn{\sum_{j\neq i} -\log p_{ij}} for coordinate
#'   \eqn{i}, where \eqn{p_{ij}} is the p-value of testing independence between \eqn{X_i} and \eqn{X_j}
#' @param lower A numeric vector of length \eqn{d} specifying the lower
#'   bounds of the current region.
#' @param upper A numeric vector of length \eqn{d} specifying the upper
#'   bounds of the current region.
#'
#' @details
#' Adaptively choosing partition dimension by first testing marginal uniformity
#' using Anderson-Darling test for each coordinate \eqn{i \in \{1, \ldots, d\}}, i.e.
#' we test if \eqn{x_{1i}, \ldots, x_{ni}} are uniformly distributed within
#' \eqn{(\texttt{lower}_i, \texttt{upper}_i)}. If the minimum p-value is less than \code{thresh_marginal},
#' randomly sample one coordinate with p-value less than \code{thresh_marginal} as the partition dimension.
#'
#' Otherwise, test pairwise independence between all \eqn{\binom{d}{2}} coordinate pairs \eqn{(i, j)} using a Fisher's exact test.
#' To test if \eqn{X_i} is independent of \eqn{X_j} in the region \eqn{R}, we divide \eqn{R} into four
#' quadrants along the middle of \eqn{i} and \eqn{j} coordinates and apply Fisher's exact test to counts of observations
#' in each quadrant. Compute a combined p-value for \eqn{i}th coordinate
#' as \eqn{p_i = \sum_{j \neq i} -\log p_{ij}} and randomly sample one coordinate with \eqn{p_i > }\code{thresh_interaction} as partition dimension.
#'
#' If no coordinate is significant in either tests, randomly pick one coordinate to partition.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{below_marginal}}{A vector of integers corresponding to coordinates
#'     whose marginal Anderson-Darling p-value is below
#'     \code{thresh_marginal}, \code{NA} if no rejection.}
#'   \item{\code{below_interaction}}{A vector of integers corresponding to coordinates
#'      whose combined interaction score \eqn{p_i} exceeds
#'     \code{thresh_interaction}, \code{NA} if no rejection.}
#'   \item{\code{p}}{Selected partition coordinate.}
#' }
#'
#' @examples
#' set.seed(42)
#'
#' # Uniform observations
#' x <- matrix(runif(300), nrow = 100, ncol = 3)
#' get_partition_dim(x,
#'                  thresh_marginal    = 0.05,
#'                  thresh_interaction = qgamma(p=1-0.05, shape=3-1, rate=1),
#'                  lower = rep(0, 3),
#'                  upper = rep(1, 3))
#'
#' # Non-uniform observations
#' y <- matrix(runif(300), nrow = 100, ncol = 3)
#' y[, 1] <- rbeta(100, 3, 3)
#' get_partition_dim(y,
#'                  thresh_marginal    = 0.05,
#'                  thresh_interaction = qgamma(p=1-0.05, shape=3-1, rate=1),
#'                  lower = rep(0, 3),
#'                  upper = rep(1, 3))
#'
#' @seealso
#' \code{\link{build_adaptive_histogram}} for the function that calls
#' \code{get_partition_dim} to adaptively choose partition locations.
#'
#' @export
get_partition_dim <- function(x, thresh_marginal, thresh_interaction,
                              lower, upper) {
  # check inputs
  if (!is.matrix(x) || !is.numeric(x))
    stop("'x' must be a numeric matrix.")
  if (nrow(x) < 2L)
    stop("'x' must have at least 2 rows.")
  if (!is.numeric(thresh_marginal) || length(thresh_marginal) != 1 ||
      thresh_marginal <= 0 || thresh_marginal >= 1)
    stop("'thresh_marginal' must be a single numeric value in (0, 1).")
  if (!is.numeric(thresh_interaction) || length(thresh_interaction) != 1 ||
      thresh_interaction <= 0)
    stop("'thresh_interaction' must be a single positive numeric value.")
  if (length(lower) != ncol(x) || length(upper) != ncol(x))
    stop("'lower' and 'upper' must each have length ncol(x).")

  d <- ncol(x)

  if(d == 1){
    return(list(below_marginal    = NA,
                below_interaction = NA,
                p                 = 1))
  }

  # test marginal uniformity
  p_marginal <- numeric(d)
  for (i in seq_len(d)) {
    p_marginal[i] <- ad.test(x[, i], null = "punif",
                             min = lower[i], max = upper[i])$p
  }

  # if reject uniform marginals for some coordinate, return one of the significant coordinates as partition dimension
  if (min(p_marginal) < thresh_marginal) {
    rejected <- which(p_marginal < thresh_marginal)
    p        <- sample(rejected, size = 1)
    return(list(below_marginal    = rejected,
                below_interaction = NA,
                p                 = p))
  }

  # otherwise, test pairwise independence
  mlog_p_interaction <- matrix(NA_real_, d, d)

  if(d > 2){
    for (i in seq_len(d - 1)) {
      for (j in (i + 1):d) {

        # split each coordinate in the middle of the region
        mid <- (lower[c(i, j)] + upper[c(i, j)]) / 2

        # count observations in each quadrant
        n_upper_left  <- sum(x[, i] < mid[1] & x[, j] > mid[2])
        n_upper_right <- sum(x[, i] > mid[1] & x[, j] > mid[2])
        n_lower_right <- sum(x[, i] > mid[1] & x[, j] < mid[2])
        n_lower_left  <- sum(x[, i] < mid[1] & x[, j] < mid[2])

        # compute 2x2 contingency table
        dat <- matrix(c(n_upper_left, n_upper_right,
                        n_lower_left, n_lower_right), 2, 2, byrow = TRUE)

        # compute Fisher's exact test p-value
        mlog_p_interaction[i, j] <-
          mlog_p_interaction[j, i] <- -log(fisher.test(dat)$p.value)
      }
    }
  }

  # combined interaction score for each coordinate
  mlog_p_interaction_combined <- rowSums(mlog_p_interaction, na.rm = TRUE)

  # if the interaction score exceeds the threshold for some coordinate, return one of the significant coordinates as partition dimension
  if (max(mlog_p_interaction_combined) > thresh_interaction) {
    rejected <- which(mlog_p_interaction_combined > thresh_interaction)
    p        <- sample(rejected, size = 1)
    return(list(below_marginal    = NA,
                below_interaction = rejected,
                p                 = p))
  }

  # if no rejection in either tests, randomly sample one coordinate as partition dimension
  p <- sample(seq_len(d), size = 1)
  return(list(below_marginal    = NA,
              below_interaction = NA,
              p                 = p))
}


#' Build an Adaptive Beta-Tree Histogram
#'
#' @description
#' Constructs an adaptive multivariate histogram using a k-d tree structure,
#' where the splitting dimension at each node is chosen adaptively by
#' \code{\link{get_partition_dim}} based on statistical tests for marginal
#' non-uniformity and pairwise interaction. The histogram produces a set of
#' rectangular regions, each with a confidence interval on the average
#' density, calibrated to control the simultaneous coverage probability
#' at level \eqn{1 - \alpha}. The extreme observations in each coordinate
#' are used as the bounding box of the histogram and removed from the
#' data before tree construction.
#'
#' @param X A numeric matrix of dimension \eqn{n \times d}, the data.
#'   Rows are observations and columns are coordinates.
#' @param thresh_marginal A numeric scalar in \eqn{(0, 1)}. The significance
#'   threshold for the marginal Anderson-Darling uniformity test inside
#'   \code{\link{get_partition_dim}}. Passed directly at every node.
#' @param thresh_interaction A positive numeric scalar. The threshold on
#'   the combined interaction score inside \code{\link{get_partition_dim}}.
#'   Passed directly at every node.
#' @param alpha A numeric scalar in \eqn{(0, 1)}. The overall significance
#'   level for the simultaneous confidence intervals on the average density
#'   in each histogram bin. Passed to \code{\link{ConfLevel}}. Defaults to
#'   \code{0.1}.
#' @param method A character string specifying the multiple testing
#'   correction method used in \code{\link{ConfLevel}}. One of
#'   \code{"bonferroni"} or \code{"weighted_bonferroni"} (default).
#'
#' @details
#' \strong{Bounding box construction:}
#'
#' Before tree construction, the minimum and maximum observation in each
#' coordinate \eqn{i} are extracted from the current data matrix and used
#' as the lower and upper bounds \eqn{(\texttt{lower}_i, \texttt{upper}_i)}
#' of the histogram's bounding box. The two corresponding observations are
#' then removed from the data. This is repeated sequentially for each
#' coordinate \eqn{i = 1, \ldots, d}, so that \eqn{2d} observations are
#' removed in total. Note that after removing the extreme observations for
#' coordinate \eqn{i}, the remaining data is used to find the extremes for
#' coordinate \eqn{i+1}; thus the bounds for later coordinates are the
#' extremes of the already-pruned dataset.
#'
#' \strong{K-d tree construction:}
#'
#' The tree is built recursively via \code{add_node}. At each node:
#' \enumerate{
#'   \item If the number of observations in the node is less than
#'     \eqn{4 \log n} (where \eqn{n} is the total sample size), the node
#'     is marked as a leaf and no further splitting is performed.
#'   \item Otherwise, \code{\link{get_partition_dim}} selects the splitting
#'     coordinate \eqn{p}, and the data are split at the median of
#'     coordinate \eqn{p} (i.e., at quantile \eqn{qt = 0.5}). The median
#'     observation itself is excluded from both children.
#'   \item The left child receives observations strictly below the median
#'     and the right child receives observations strictly above the median.
#' }
#' The depth-wise region counts \code{nd} are accumulated during traversal
#' and passed to \code{\link{ConfLevel}} to determine the per-node
#' confidence level \code{ahat} adjusted for multiple testing.
#'
#' \strong{Confidence bounds and region selection:}
#'
#' After tree construction, \code{\link{SetBounds}} attaches lower and
#' upper confidence bounds on the average density to each node.
#' \code{\link{SelectNodes}} then extracts the selected regions into a
#' matrix \code{B} of dimension \eqn{m \times (2d + 5)}, where each row
#' encodes the bounds, density estimate, confidence bounds, observation
#' count, and depth of one selected region.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{\code{split_dir}}{An integer vector recording the splitting
#'     coordinate chosen at each internal node of the k-d tree, in the
#'     order the nodes were processed. The first element is always
#'     \code{NA} (initialisation artefact).}
#'   \item{\code{hist}}{A numeric matrix of dimension
#'     \eqn{m \times (2d + 5)}, where each row corresponds to a selected
#'     histogram region. The first \eqn{d} columns give the lower bounds,
#'     the next \eqn{d} columns give the upper bounds, and the final 5
#'     columns give the estimated average density, lower confidence bound,
#'     upper confidence bound, observation count, and tree depth.}
#'   \item{\code{lower}}{A numeric vector of length \eqn{d} giving the
#'     lower bounds of the bounding box (the minimum observation in each
#'     coordinate before sequential pruning).}
#'   \item{\code{upper}}{A numeric vector of length \eqn{d} giving the
#'     upper bounds of the bounding box.}
#' }
#'
#' @examples
#' set.seed(42)
#'
#' # Bivariate standard normal data
#' X <- matrix(rnorm(1000), nrow = 500, ncol = 2)
#'
#' hist_obj <- build_adaptive_histogram(
#'   X                 = X,
#'   thresh_marginal   = 0.05,
#'   thresh_interaction = 10,
#'   alpha             = 0.1,
#'   method            = "weighted_bonferroni"
#' )
#'
#' # Bounding box
#' hist_obj$lower
#' hist_obj$upper
#'
#' # Number of selected histogram regions
#' nrow(hist_obj$hist)
#'
#' # Splitting directions used
#' table(hist_obj$split_dir)
#'
#' @seealso
#' \code{\link{get_partition_dim}} for the adaptive splitting rule,
#' \code{\link{ConfLevel}} for the multiple-testing-adjusted confidence
#' level,
#' \code{\link{SetBounds}} for attaching density confidence bounds to
#' nodes,
#' \code{\link{SelectNodes}} for extracting selected regions from the
#' tree.
#'
#' @export
build_adaptive_histogram <- function(X, thresh_marginal, thresh_interaction,
                                     alpha  = 0.1,
                                     method = "weighted_bonferroni") {

  # ---- Input validation ----
  if (!is.matrix(X) || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
  if (nrow(X) < 4)
    stop("'X' must have at least 4 rows.")
  if (!is.numeric(thresh_marginal) || thresh_marginal <= 0 || thresh_marginal >= 1)
    stop("'thresh_marginal' must be a numeric value in (0, 1).")
  if (!is.numeric(thresh_interaction) || thresh_interaction <= 0)
    stop("'thresh_interaction' must be a single positive numeric value.")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a numeric value in (0, 1).")
  if (!method %in% c("bonferroni", "weighted_bonferroni"))
    stop("'method' must be one of \"bonferroni\" or \"weighted_bonferroni\".")

  nd        <- NA_integer_   # region count at each tree depth level
  split_dir <- NA_integer_   # splitting coordinate at each internal node
  d         <- ncol(X)
  n         <- nrow(X)

  # ---- Step 1: Extract bounding box and remove extreme observations ----
  # For each coordinate i, the min and max observations define the bounding
  # box and are removed from the dataset. Removal is sequential: bounds for
  # coordinate i+1 are found after removing the extremes for coordinate i.
  lower <- numeric(d)
  upper <- numeric(d)
  for (i in seq_len(d)) {
    Xi      <- X[, i]
    ind_min <- which.min(Xi)
    ind_max <- which.max(Xi)
    lower[i] <- Xi[ind_min]
    upper[i] <- Xi[ind_max]
    X        <- X[-c(ind_min, ind_max), , drop = FALSE]
  }

  # ---- Step 2: Initialise root node ----
  rootnode <- list(
    leftchild  = NULL,
    rightchild = NULL,
    ndat       = nrow(X),   # number of observations in this node
    depth      = 0L,
    low        = lower,     # lower bounds of this node's region
    up         = upper,     # upper bounds of this node's region
    lower      = NA_real_,  # lower confidence bound for average density
    upper      = NA_real_,  # upper confidence bound
    leaf       = FALSE,
    bounded    = TRUE,
    dir        = NA_integer_
  )

  # ---- Step 3: Recursively build the k-d tree ----
  add_node <- function(x, node, d, thresh_marginal, thresh_interaction) {

    # Accumulate region count at this depth level
    # Fixed bug: original code only set nd[depth+1] = 1 the first time a
    # depth level was encountered, never incrementing for subsequent nodes
    # at the same depth. Fix: increment if already initialised.
    if (is.na(nd[node$depth + 1])) {
      nd[node$depth + 1] <<- 1L
    } else {
      nd[node$depth + 1] <<- nd[node$depth + 1] + 1L
    }

    if (node$ndat < 4 * log(n)) {

      # Too few observations: mark as leaf and return
      node$leaf <- TRUE

    } else {

      # Enough observations: split this node
      leftnode  <- node
      rightnode <- node

      # Choose splitting coordinate adaptively
      dir        <- get_partition_dim(x,
                                      lower              = node$low,
                                      upper              = node$up,
                                      thresh_marginal    = thresh_marginal,
                                      thresh_interaction = thresh_interaction)
      node$dir   <- dir$p
      p          <- dir$p
      qt         <- 0.5   # always split at the median

      split_dir <<- c(split_dir, dir$p)

      # Update depth for children
      depth             <- node$depth + 1L
      leftnode$depth    <- depth
      rightnode$depth   <- depth

      # Sort data along the splitting coordinate and split at the median
      x <- x[order(x[, p]), , drop = FALSE]
      m <- node$ndat
      median_pos <- ceiling(qt * m)   # index of the median observation

      leftnode$ndat  <- median_pos - 1L
      rightnode$ndat <- m - median_pos

      # Median observation is excluded from both children
      xleft  <- x[seq_len(median_pos - 1),        , drop = FALSE]
      xright <- x[(median_pos + 1):m,              , drop = FALSE]

      # Update region bounds: the median value becomes the new boundary
      leftnode$up[p]   <- x[median_pos, p]
      rightnode$low[p] <- x[median_pos, p]

      node$leftchild  <- add_node(xleft,  leftnode,  d,
                                  thresh_marginal, thresh_interaction)
      node$rightchild <- add_node(xright, rightnode, d,
                                  thresh_marginal, thresh_interaction)
    }

    return(node)
  }

  kdtree <- add_node(X, rootnode, d,
                     thresh_marginal    = thresh_marginal,
                     thresh_interaction = thresh_interaction)

  # ---- Step 4: Set confidence bounds and select histogram regions ----
  ahat   <- ConfLevel(nd, alpha, method)
  kdtree <- SetBounds(kdtree, ahat, n)

  # Each row of B encodes one selected region: [lower(d), upper(d), density, CI_lower, CI_upper, n, depth]
  B <- matrix(nrow = 0, ncol = (2 * d + 5))
  B <- SelectNodes(kdtree, B, ahat, n)

  return(list(split_dir = split_dir,
              hist      = B,
              lower     = lower,
              upper     = upper))
}
