# Are two regions distinct modes?

Two regions `i` and `j` are two *distinct* modes if for *every* path
connecting `i` and `j`, there exists at least one region `R` along the
path, such that the upper endpoint of the confidence interval of the
average density of `R` is *below* the lower endpoints for both `i` and
`j`.

## Usage

``` r
is_connected(i, j, g, ci, cutoff = 6)
```

## Arguments

- i, j:

  Test if the \\i\\-th and the \\j\\-th regions are distinct modes or
  not.

- g:

  A graph from the adjacency matrix returned by
  [graph_from_adjacency_matrix](https://r.igraph.org/reference/graph_from_adjacency_matrix.html).

- ci:

  A matrix of size N\*2 of the confidence intervals for the average
  densities. Each row stores the lower and upper confidence limits of
  one region in the Beta tree.

- cutoff:

  Maximum path length. By default, we only evaluate paths with length at
  most `cutoff = 6`.

## Value

`unconnected` if the two regions are distinct modes and `connected`
otherwise.

## Details

For computational reasons, we only check paths whose length is at most
`cutoff`.
