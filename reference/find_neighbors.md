# Find Regions in the Neighborhood of a Mode

Computes the set of histogram regions that form the neighborhood of a
given mode.

## Usage

``` r
find_neighbors(mode, g, ci)
```

## Arguments

- mode:

  A single positive integer giving the row index of the mode region in
  the histogram matrix.

- g:

  An `igraph` undirected graph object representing the adjacency
  structure of the histogram regions, as produced by
  [`graph_from_adjacency_matrix`](https://r.igraph.org/reference/graph_from_adjacency_matrix.html)
  from the output of
  [`compute_adjacency_mat`](https://zq00.github.io/BetaTree/reference/compute_adjacency_mat.md).

- ci:

  A numeric matrix with two columns. Row \\i\\ gives the lower
  (column 1) and upper (column 2) confidence bounds for the average
  density of histogram region \\i\\.

## Value

An integer vector of row indices (including `mode`) of the regions
belonging to the neighbourhood of `mode`.

## Details

A region \\R\\ is considered to be in the neighborhood of a model \\M\\
if there exists a path connecting \\R\\ and \\M\\ such that the upper
confidence bound of every region along the path is above the lower
confidence bound of \\M\\.

## See also

[`FindModesApproximate`](https://zq00.github.io/BetaTree/reference/FindModesApproximate.md)
for the function that calls `find_neighbors`.
