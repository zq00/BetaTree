#' Build a k-d tree
#'
#' Grow a k-d tree by iteratively partitioning along the sample median of each coordinate.
#'
#' \code{BuildKDTree} constructs a kd-tree from data by iteratively partitioning along the sample median of the \eqn{j}-th coordinate, where \eqn{j} cycles through \eqn{1,2,\ldots, d}.
#'      Stop partitioning a node when the number of obs. inside is less than \eqn{4\log(n)}.
#'
#'@param X  A data matrix of size n by d.
#'@param bounded If \code{bounded = TRUE}, then create an initial bounding box according to user input parameters:
#'\describe{
#'\item{option}{Define the bounding box by the order statistics (\code{option = ndat}) or quantiles (\code{option = qt}) of observations to be excluded.}
#'\item{qt}{A vector of length \code{d} equal to the data dimension. Order statistics or quantiles of obs. to be excluded in each dimension. E.g. for 2-dim data, \code{option = "ndat"} and \code{qt = c(50, 40)}, then the boundary in the x-axis are the 50-th and (n-50 +1)-th order statistics in the x-axis.
#'    The boundary in the y-axis is defined by the 40-th and \eqn{n_1 - 40 + 1} order statistics in the y-axis for the obs. inside the boundary defined by the x-axis in the previous step.
#'    Another way to put it is to exclude 50 obs. at the lower and upper end of x-coordinate, then for the obs. within, exclude 40 obs. at the lower and upper end of y-coordinate. If \code{option = "qt"} and \code{qt = c(0.05, 0.1)},
#'    then the boundaries are defined by the \eqn{\lceil 0.05 n \rceil} and \eqn{n - \lceil 0.05 n \rceil + 1} quantiles in the x-axis, and then the
#'     \eqn{\lceil 0.1 n_1 \rceil} and \eqn{n_1 - \lceil 0.1 n_1 \rceil + 1} quantiles in the y-axis among observations in the first boundary defined by x-axis.
#'     Equivalently, exclude  \eqn{\lceil 0.05 n \rceil} obs. at the lower and upper end of x-coordinate, then for the obs. within, exclude  \eqn{\lceil 0.1 n_1 \rceil} obs. at the lower and upper end of the y-coordinate.
#'    }
#'}
#'@param ... Additional parameters if \code{bounded = TRUE}.
#'
#' @return \code{build_tree} returns two values:
#'\describe{
#'\item{kdtree}{the k-d tree built from the root node. Each node contains the following information:
#' \describe{
#' \item{depth}{Tree depth of the region. The depth of root node is 0.}
#' \item{ndat}{number of obs. in this region.}
#' \item{low, up}{lower and upper bound of the region, if it exists.}
#' \item{lower, upper}{lower and upper confidence bound for the average density.}
#' \item{bounded}{whether the region is bounded. A region is bounded if all of (\code{low}, \code{up}) exist.}
#' \item{leaf}{whether the node is leaf or not. A node is a leaf if the number of obs. is less than \eqn{4 \log n}, where \eqn{n} is the total number of obs. Stop partitioning if a node is leaf. }
#' \item{leftchild, rightchild}{pointers to the left and right child of a node. \code{leftchild} points to the node representing the region less than the median in the partitioning dimension, and \code{rightchild} points to the node representing the region greater than the median.}
#' }
#' }
#'\item{nd}{The number of bounded regions at every level, starting at level 0.}
#'}
#'
#'@export
#'@examples
#'# build a k-d tree for standard normal samples
#' X <- matrix(rnorm(1000), ncol = 2)
#' kd <- BuildKDTree(X)
#' # Initialize a bounding box
#' kd <-  BuildKDTree(X, bounded = TRUE, option = "ndat", qt = c(1,1))
BuildKDTree <- function(X, bounded = F, ...){

  n <- nrow(X) # No. obs
  d <- ncol(X) # No. dim

  # Initialize a root node
  if(bounded){
    params <-  list(...)
    if(is.null(params$option)){stop("Please define bounding box (option = 'ndat' or 'qt')")}
    if(is.null(params$q)){stop("Please input number/quantile of points excluded in the bounding box")}

    option <- params$option
    q <- params$q
    Xinside <- X; nn <- nrow(Xinside) # No. obs inside the bounding box
    bounds <- matrix(NA, nrow = d, ncol = 2) # lower bound and upper bound of the bounding box (Columns)
    for(i in 1:d){ # iterate through each dimension
      Xsorted <- Xinside[order(Xinside[,i]),] # sort X according to i-th dimension
      bounds[i, 1] <-  switch(option,
                              "ndat" = Xsorted[q[i], i],
                              "qt" = Xsorted[ceiling(nn*q[i]), i])
      bounds[i, 2] <-  switch(option,
                              "ndat" = Xsorted[nn-q[i] + 1, i],
                              "qt" = Xsorted[nn - ceiling(nn*q[i]) + 1, i])
      Xinside <- switch(option,
                        "ndat" = Xsorted[(q[i] + 1):(nn-q[i]), , drop = F],
                        "qt" = Xsorted[(ceiling(nn*q[i]) + 1):(nn - ceiling(nn*q[i])), ,drop = F])
      nn <- nrow(Xinside)
    }

    rootnode <- list(
      leftchild = NULL,
      rightchild = NULL,
      ndat = nrow(Xinside),   # number of observ. in this node
      depth = 0,
      low = bounds[,1],   # vector of lower bound on region
      up  =  bounds[,2],
      lower = NA,  # lower confidence bound for average density
      upper = NA,  # upper confidence bound
      leaf  = FALSE,
      bounded = TRUE)      # node is leaf or not
  }else{
    rootnode <- list(
      leftchild = NULL,
      rightchild = NULL,
      ndat = n,
      depth = 0,
      low = rep(NA, d),
      up  =  rep(NA, d),
      lower = NULL,
      upper = NULL,
      bounded = FALSE,
      leaf  = FALSE)

    Xinside <- X
  }

  # Inner function to add children of the input node
  # input: X - observ. in current node
  #        node - a list element
  #        n - total number of observ.
  #        d - dimension
  # output: add two elements to kdtree
  add_node <- function(x, node, n, d){
    node$bounded = prod(sapply(c(node$low, node$up), function(t) !is.na(t)))

    if(is.na(nd[node$depth + 1])){ # the first region at depth d
      nd[node$depth + 1] <<- 0 # initialize with no bounded region
    }
    if(node$bounded){ # if the node is bounded, then count the region
      nd[node$depth + 1] <<- nd[node$depth + 1] + 1
    }

    if (node$ndat < 4*log(n)) { node$leaf = TRUE   # node is leaf, return it marked as such
    } else {  # split this node:
      leftnode = node; rightnode = node   # initialize by taking info from parent
      p = node$depth %% d + 1 # partition dimension

      depth = node$depth + 1
      leftnode$depth = depth; rightnode$depth = depth

      x = x[order(x[,p]),, drop = F] # sort according to the pth coordinate
      m = node$ndat  # number of data in x
      leftnode$ndat = ceiling(m/2)-1; rightnode$ndat = m-ceiling(m/2)
      xleft = x[1:(ceiling(m/2)-1),, drop = F]; xright = x[(ceiling(m/2)+1):m,, drop = F]

      leftnode$up[p] <- x[ceiling(m/2),p, drop = F]
      rightnode$low[p] <- x[ceiling(m/2),p, drop = F]

      node$leftchild = add_node(xleft,leftnode, n, d )
      node$rightchild = add_node(xright,rightnode, n, d)
    }
    return(node)
  }

  nd <- NA # No. regions in each level
  kdtree <- add_node(Xinside, rootnode, n, d)

  return(list(kdtree = kdtree, nd = nd))
}
