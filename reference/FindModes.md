# Identifies distinct modes in Beta-trees histogram

Compute modes in the histogram.

## Usage

``` r
FindModes(hist, d, cutoff = 6)
```

## Arguments

- hist:

  A matrix. Each row represents one region in the histogram. If the
  histogram is constructed for `d` dimensional data, then the first `d`
  columns contain the lower bounds for each coordinate, and the next `d`
  columns contain the upper bounds. The `2d + 1` column contains the
  empirical density. `hist` can be the output of
  [BuildHist](https://zq00.github.io/BetaTree/reference/BuildHist.md) or
  [SelectNodes](https://zq00.github.io/BetaTree/reference/SelectNodes.md)
  function.

- d:

  A numeric value of the data dimension.

- cutoff:

  Maximum path length. By default, we only evaluate paths with length at
  most `cutoff = 6`.

## Value

- mode:

  A vector of candidate modes.

- hist:

  The input Beta-tree histogram.

- g:

  A graph object constructed using the adjacency matrix of the regions
  in the histogram.

## Details

Two regions `i` and `j` are distinct modes if along every path
connecting them, there exists one region whose upper confidence limit is
below the lower confidence limit of both `i` and `j`. This function
iterates through all the regions in the Beta-tree histogram in
descending order of empirical density and marks a region as a mode if it
is a distinct mode from all currently discovered distinct modes.

## Examples

``` r
# A mixture of two Gaussians
X <-  matrix(rnorm(6000, 0, 1), ncol = 2)
X[1:1500, 1] <- X[1:1500, 1] + 4
rect <- BuildHist(X) # construct a Beta-tree histogram
modes <- FindModes(hist = rect, d = 2)$mode
rect[modes, ] # show information of the modes
#>            [,1]        [,2]        [,3]      [,4]       [,5]       [,6]
#> [1,]  4.0002715 -0.02163661  4.56182126 0.6199489 0.08696886 0.06169925
#> [2,] -0.6307679  0.02554109 -0.02458833 0.7152718 0.07494208 0.05316696
#>           [,7] [,8] [,9]
#> [1,] 0.1173330   93    5
#> [2,] 0.1011072   93    5
```
