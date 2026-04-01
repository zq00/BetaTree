# Compute the adjacency matrix of regions in the histogram

Computes the adjacency matrix of regions in a Beta tree histogram

## Usage

``` r
compute_adjacency_mat(hist, d)
```

## Arguments

- hist:

  A matrix. Each row represents one region in the histogram. If the
  histogram is constructed for `d` dimensional data, then the first `d`
  columns contain the lower bounds for each coordinate, and the next `d`
  columns contain the upper bounds. The `2d + 1` column contains the
  empirical density. `hist` can be the output of
  [BuildHist](BuildHist.md) or [SelectNodes](SelectNodes.md) function.

- d:

  A numeric value of the data dimension.

## Value

A matrix A of size N\*N, where N is the number of regions in the
histogram, i.e. number of rows in `hist`. \\A\_{i,j} = 1\\ if the i-th
and j-th regions are adjacent and 0 otherwise.

## Details

Two regions are neighbors if and only if the closures of the intervals
in each direction intersect. For example, let the two regions \\X\\
defined by \\(X\_{i,\mathrm{low}}, X\_{i,\mathrm{up}})\\ and \\Y\\
defined by \\(Y\_{i,\mathrm{low}}, Y\_{i,\mathrm{up}})\\ (for
\\i=1,\ldots, d\\), where \\X\_{i,\mathrm{low}}\\ and
\\X\_{i,\mathrm{up}}\\ are the lower and upper bounds in dimension \\i\\
for the region \\X\\. Then, \\X\\ and \\Y\\ are neighbors if and only if
\\\[X\_{i,\mathrm{low}}, X\_{i,\mathrm{up}}\] \cap
\[Y\_{i,\mathrm{low}}, Y\_{i,\mathrm{up}}\] \neq \emptyset\\ for every
\\i=1,\ldots, d\\.

## Examples

``` r
X <-  matrix(rnorm(2000, 0, 1), ncol = 2)
B <- BuildHist(X)
adj <- compute_adjacency_mat(hist = B, d = 2)
```
