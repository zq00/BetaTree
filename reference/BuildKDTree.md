# Build a k-d tree

Grow a k-d tree by iteratively partitioning along the sample median of
each coordinate.

## Usage

``` r
BuildKDTree(X, bounded = F, ...)
```

## Arguments

- X:

  A data matrix of size n by d.

- bounded:

  If `bounded = TRUE`, then create an initial bounding box according to
  user input parameters:

  option

  :   Define the bounding box by the order statistics (`option = ndat`)
      or quantiles (`option = qt`) of observations to be excluded.

  qt

  :   A vector of length `d` equal to the data dimension. Order
      statistics or quantiles of obs. to be excluded in each dimension.
      E.g. for 2-dim data, `option = "ndat"` and `qt = c(50, 40)`, then
      the boundary in the x-axis are the 50-th and (n-50 +1)-th order
      statistics in the x-axis. The boundary in the y-axis is defined by
      the 40-th and \\n_1 - 40 + 1\\ order statistics in the y-axis for
      the obs. inside the boundary defined by the x-axis in the previous
      step. Another way to put it is to exclude 50 obs. at the lower and
      upper end of x-coordinate, then for the obs. within, exclude 40
      obs. at the lower and upper end of y-coordinate. If
      `option = "qt"` and `qt = c(0.05, 0.1)`, then the boundaries are
      defined by the \\\lceil 0.05 n \rceil\\ and \\n - \lceil 0.05 n
      \rceil + 1\\ quantiles in the x-axis, and then the \\\lceil 0.1
      n_1 \rceil\\ and \\n_1 - \lceil 0.1 n_1 \rceil + 1\\ quantiles in
      the y-axis among observations in the first boundary defined by
      x-axis. Equivalently, exclude \\\lceil 0.05 n \rceil\\ obs. at the
      lower and upper end of x-coordinate, then for the obs. within,
      exclude \\\lceil 0.1 n_1 \rceil\\ obs. at the lower and upper end
      of the y-coordinate.

- ...:

  Additional parameters if `bounded = TRUE`.

## Value

`build_tree` returns two values:

- kdtree:

  the k-d tree built from the root node. Each node contains the
  following information:

  depth

  :   Tree depth of the region. The depth of root node is 0.

  ndat

  :   number of obs. in this region.

  low, up

  :   lower and upper bound of the region, if it exists.

  lower, upper

  :   lower and upper confidence bound for the average density.

  bounded

  :   whether the region is bounded. A region is bounded if all of
      (`low`, `up`) exist.

  leaf

  :   whether the node is leaf or not. A node is a leaf if the number of
      obs. is less than \\4 \log n\\, where \\n\\ is the total number of
      obs. Stop partitioning if a node is leaf.

  leftchild, rightchild

  :   pointers to the left and right child of a node. `leftchild` points
      to the node representing the region less than the median in the
      partitioning dimension, and `rightchild` points to the node
      representing the region greater than the median.

- nd:

  The number of bounded regions at every level, starting at level 0.

## Details

`BuildKDTree` constructs a kd-tree from data by iteratively partitioning
along the sample median of the \\j\\-th coordinate, where \\j\\ cycles
through \\1,2,\ldots, d\\. Stop partitioning a node when the number of
obs. inside is less than \\4\log(n)\\.

## Examples

``` r
# build a k-d tree for standard normal samples
X <- matrix(rnorm(1000), ncol = 2)
kd <- BuildKDTree(X)
# Initialize a bounding box
kd <-  BuildKDTree(X, bounded = TRUE, option = "ndat", qt = c(1,1))
```
