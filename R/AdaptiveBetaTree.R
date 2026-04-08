#' Choose Partition Dimension
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
#' we test if \eqn{x_{1i}, \ldots, x_{ni}} are uniformly distributed within the lower and upper bounds of coordinate \eqn{i}.
#' If the minimum p-value is less than \code{thresh_marginal},
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
#' Constructs an Beta-tree histogram there partition coordinates are chosen adaptively
#' based on tests of marginal uniformity and pairwise independence.
#'
#' @param X A numeric matrix dimension \eqn{n \times d} containing observations falling inside the current
#'   region. Each row corresponds to one observation.
#' @param thresh_marginal A numeric value between \eqn{(0, 1)}. The significance
#'   threshold for the marginal Anderson-Darling uniformity test. Default value is \code{alpha / d}.
#' @param thresh_interaction A positive numeric value. The threshold on the
#'   combined interaction score \eqn{\sum_{j\neq i} -\log p_{ij}} for coordinate
#'   \eqn{i}, where \eqn{p_{ij}} is the p-value of testing independence between \eqn{X_i} and \eqn{X_j}.
#'   Default value is \eqn{(1-\alpha)} quantile of a Gamma variable with shape \eqn{d-1} and rate 1.
#' @param alpha Significance level, default is \code{alpha = 0.1}.
#'   The Beta-tree histogram provides simultaneous confidence intervals for average densities
#'   in each histogram region with coverage probability \eqn{1-\alpha}.
#' @param method A character string specifying whether to use \code{method = "bonferroni"} or \code{method = "weighted_bonferroni"} method to
#'  adjust for multiple hypothesis testing when computing significance levels.
#'
#' @details The adaptive Beta-tree histogram uses the same algorithm as the Beta-tree histogram (see \code{\link{BuildHist}} function
#' and package Vignettes \code{beta_trees}) with two differences:
#'
#' First, the adaptive Beta-tree automatically computes a bounded histogram by placing the first splits at the
#' minimum and maximum of each coordinates.
#'
#' Second, the adaptive Beta-tree uses \code{\link{get_partition_dim}} to select partition dimension
#' at each step. The same set of thresholds \code{thresh_marginal} and \code{thresh_interaction} are used at every step.
#'
#' @return A named list with three components:
#' \describe{
#'   \item{\code{split_dir}}{An integer vector recording the splitting
#'     coordinate chosen at each step. The first element is always \code{NA}. }
#'   \item{\code{hist}}{ A numeric matrix of dimension \eqn{m \times (2d + 5)}.
#'   Each row represents one region in the Beta-tree histogram.
#'    The first \eqn{d} columns are the lower bounds for the region; the next \eqn{d} columns are the upper bounds.
#'    The final 5 columns give the empirical density, lower and upper confidence bounds of average density, ,
#'     upper confidence bound, number of observations inside, and tree depth of each region. See also \code{\link{BuildHist}}.}
#'   \item{\code{lower}}{A numeric vector of length \eqn{d} corresponding to the
#'     the minimum observation in each coordinate.}
#'   \item{\code{upper}}{A numeric vector of length \eqn{d} corresponding to the
#'     the maximum observation in each coordinate.}
#' }
#'
#' @examples
#' set.seed(42)
#'
#' # Bivariate standard normal data
#' X <- matrix(rnorm(1000), nrow = 500, ncol = 2)
#'
#' adaptive_beta_tree <- build_adaptive_histogram(
#'   X                 = X,
#'   alpha             = 0.1,
#'   method            = "weighted_bonferroni"
#' )
#'
#' # Number of histogram regions
#' nrow(adaptive_beta_tree$hist)
#'
#' # Splitting directions
#' table(adaptive_beta_tree$split_dir)
#'
#' @seealso
#' \code{\link{get_partition_dim}} adaptively selectx partition coordinates.
#' \code{\link{BuildHist}} builds a Beta-tree histogram.
#'
#' @export
build_adaptive_histogram <- function(X,
                                     alpha  = 0.1,
                                     thresh_marginal = alpha / ncol(X),
                                     thresh_interaction = qgamma(shape = ncol(X) - 1,rate = 1, p = 1 - alpha),
                                     method = "weighted_bonferroni") {

  # ---- Input validation ----
  if (!is.matrix(X) || !is.numeric(X))
    stop("'X' must be a numeric matrix.")
  if (!is.numeric(thresh_marginal) || thresh_marginal <= 0 || thresh_marginal >= 1)
    stop("'thresh_marginal' must be a numeric value in (0, 1).")
  if (!is.numeric(thresh_interaction) || thresh_interaction <= 0)
    stop("'thresh_interaction' must be a single positive numeric value.")
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a numeric value in (0, 1).")
  if (!method %in% c("bonferroni", "weighted_bonferroni"))
    stop("'method' must be one of \"bonferroni\" or \"weighted_bonferroni\".")

  nd        <- NA_integer_
  split_dir <- NA_integer_
  d         <- ncol(X)
  n         <- nrow(X)

  # sequentially set the bounding box in each coordinate
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

  # initialize root node
  rootnode <- list(
    leftchild  = NULL,
    rightchild = NULL,
    ndat       = nrow(X),
    depth      = 0L,
    low        = lower,
    up         = upper,
    lower      = NA_real_,  # lower confidence bound of average density
    upper      = NA_real_,  # upper confidence bound
    leaf       = FALSE,
    bounded    = TRUE,
    dir        = NA_integer_
  )

  # recursive partitioning based on adaptively chosen partition coordinate
  add_node <- function(x, node, d, thresh_marginal, thresh_interaction) {

    if (is.na(nd[node$depth + 1])) {
      nd[node$depth + 1] <<- 1L
    } else {
      nd[node$depth + 1] <<- nd[node$depth + 1] + 1L
    }

    if (node$ndat < 4 * log(n)) {
      node$leaf <- TRUE
    } else {
      leftnode  <- node
      rightnode <- node

      # adaptively choose partition coordinate
      dir        <- get_partition_dim(x,
                                      lower              = node$low,
                                      upper              = node$up,
                                      thresh_marginal    = thresh_marginal,
                                      thresh_interaction = thresh_interaction)
      node$dir   <- dir$p # store partition coordinate
      p          <- dir$p # partition coordinate
      qt         <- 0.5   # partition at the median

      split_dir <<- c(split_dir, dir$p)

      # update depth for children
      depth             <- node$depth + 1L
      leftnode$depth    <- depth
      rightnode$depth   <- depth

      # sort data along the splitting coordinate and split at the median
      x <- x[order(x[, p]), , drop = FALSE]
      m <- node$ndat
      median_pos <- ceiling(qt * m)

      leftnode$ndat  <- median_pos - 1L
      rightnode$ndat <- m - median_pos

      # keep only observations inside each region
      xleft  <- x[seq_len(median_pos - 1),        , drop = FALSE]
      xright <- x[(median_pos + 1):m,              , drop = FALSE]

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

  # set confidence bounds for each region
  ahat   <- ConfLevel(nd, alpha, method)
  kdtree <- SetBounds(kdtree, ahat, n)

  # select histogram regions
  B <- matrix(nrow = 0, ncol = (2 * d + 5))
  B <- SelectNodes(kdtree, B, ahat, n)

  return(list(split_dir = split_dir,
              hist      = B,
              lower     = lower,
              upper     = upper))
}
