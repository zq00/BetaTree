#' Build a k-d tree
#'
#' Grow a kd-tree by iteratively partitioning along the sample median of each coordinate.
#' #' NOTE: Here I assume we start at depth 1 at all times!!
#'
#' \code{BuildKDTree} constructs a kd-tree from data by iteratively partitioning along the sample median each coordinate.
#' Each node in the kd-tree represents a region, and it consists of the following information:
#' \describe{
#' \item{depth}{level of the region. The depth of root node is 1.}
#' \item{ndat}{number of samples in this region.}
#' \item{low, up}{lower and upper bound on region, if it exists.}
#' \item{lower, upper}{lower and upper confidence bound for average density.}
#' \item{bounded}{whether the region is bounded. A region is bounded if all of (\code{low}, \code{up}) exist.}
#' \item{leaf}{whether the node is leaf or not. A node is a leaf if the number of observations is less than \eqn{4 \log n}, where \eqn{n} is the total number of samples. A node is partitioned if it is not a leaf. }
#' \item{leftchild, rightchild}{pointers to the left and right child of a node. \code{leftchild} points to the node representing the region of observations \emph{below} median in partition dimension, and \code{rightchild} points to the node representing the region of observations \emph{above} the median.}
#' }
#'
#'@param X data matrix with n rows and d columns
#'@param bounded if bounded = TRUE, then creates an initial bounding box according to user input parameters:
#'\describe{
#'\item{option}{Define the bounding box by the number of obs. (option = ndat) or quantiles (option = qt) of observations included.}
#'\item{q}{number/quantile of points to be excluded in each end of a dimension. For example, if data is two-dimensional, option = "ndat" and qt = 50, then exclude 50 obs. in the lower and upper end of x-coordinate, then for the observations within, exclude 50 obs. in the lower and upper end of y-coordinate.  }
#'}
#'
#'@return \code{build_tree} returns two values.
#'\describe{
#'\item{kdtree} the k-d tree built from root node.
#'\item{nd} The number of bounded regions at every level, starting at level 1.
#'}
#'
#'@examples
#'# build a k-d tree for standard normal samples
#'X = matrix(rnorm(20000, 0, 1), ncol = 2)
#'kd = BuildKDTree(X)
#'
#'
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
      depth = 1,
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
      depth = 1,
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

    if(is.na(nd[node$depth])){nd[node$depth] <<- 1
    }else{nd[node$depth] <<- nd[node$depth] + 1}

    if (node$ndat < 4*log(n)) { node$leaf = TRUE   # node is leaf, return it marked as such
    } else {  # split this node:
      leftnode = node; rightnode = node   # initialize by taking info from parent
      depth = node$depth + 1
      leftnode$depth = depth; rightnode$depth = depth

      p = depth %% d +1 # partition dimension
      x = x[order(x[,p]),, drop = F] # sort according to the pth coordinate
      m = node$ndat  # number of data in x
      leftnode$ndat = ceiling(m/2)-1; rightnode$ndat = m-ceiling(m/2)
      xleft = x[1:(ceiling(m/2)-1),, drop = F]; xright = x[(ceiling(m/2)+1):m,, drop = F]

      leftnode$up[p] <- x[ceiling(m/2),p, drop = F]; rightnode$low[p] <- x[ceiling(m/2),p, drop = F]

      node$leftchild = add_node(xleft,leftnode, n, d )
      node$rightchild = add_node(xright,rightnode, n, d)
    }
    return(node)
  }

  nd <- NA # No. regions in each level
  kdtree <- add_node(Xinside, rootnode, n, d)

  return(list(kdtree = kdtree, nd = nd))
}
