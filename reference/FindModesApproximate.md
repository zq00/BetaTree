# Identify Modes in a Beta-tree Histogram via Monte Carlo approximation

Implements an approximate mode-hunting algorithm by using randomly
sampled paths to estimate whether a region is a distinct mode.

## Usage

``` r
FindModesApproximate(hist, d, L, B, verbose = FALSE)
```

## Arguments

- hist:

  A numeric matrix of dimension \\m\times (2d+5)\\. Each row represents
  one region in the histogram. If the histogram is constructed for `d`
  dimensional data, then the \\1:d\\ columns contain the lower bounds
  for each coordinate, columns \\(d+1):(2d)\\ contain the upper bounds.
  The \\2d + 1\\ column contains the empirical density. The
  \\(2d+2):(2d+3)\\ columns contain the lower and upper confidence
  bounds of the empirical density. `hist` can be the output of
  [BuildHist](https://zq00.github.io/BetaTree/reference/BuildHist.md) or
  [build_adaptive_histogram](https://zq00.github.io/BetaTree/reference/build_adaptive_histogram.md)
  function.

- d:

  A positive integer corresponding to the data dimension.

- L:

  A positive integer. The length of random path.

- B:

  A positive integer. The number of independent random paths used in
  [`connected_to_modes`](https://zq00.github.io/BetaTree/reference/connected_to_modes.md)
  to assess whether a candidate region is connected to any of the
  current modes.

- verbose:

  Prints progress if `TRUE`.

## Value

A named list with three components:

- `mode`:

  An integer vector of row indices of `hist` identifying the detected
  modes.

- `hist`:

  The input `hist` matrix.

- `cluster`:

  A list of length `length(mode)`. The \\k\\-th element is an integer
  vector of row indices of histogram regions that belonging to the
  neighborhood of mode \\k\\. That is, there is a path from this region
  to the mode that the upper confidence bound of every region along the
  path (including the region itself) is higher than the lower confidence
  bound of the mode.

## Details

Evaluate every region in order of decreasing density whether it is a
distinct mode compared to current modes. In the exact algorithm
[FindModes](https://zq00.github.io/BetaTree/reference/FindModes.md), a
region \\R\\ is considered a distinct mode from a current model \\M\\ if
if along every path connecting \\R\\ and \\M\\, there exists at least
one region whose upper confidence bound of density is below the lower
confidence bounds of both \\R\\ and \\M\\. The approximate algorithm
here evaluates `B` random paths of length `L` instead of evaluating
every path.

## Examples

``` r
if (FALSE) { # \dontrun{
# mixture of two Gaussian distribution
X <- rbind(matrix(rnorm(500, mean = -2), ncol = 2),
           matrix(rnorm(500, mean =  2), ncol = 2))
beta_tree_hist <- BuildHist(X)

# find modes using approximate algorithm
modes <- FindModesApproximate(hist  = beta_tree_hist,
                              d     = 2,
                              L     = 50,
                              B     = 1000,
                              verbose = F)

# number of detected modes
length(modes$mode)
} # }
```
