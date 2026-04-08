# Is a Region a Distinct Mode? (Approximate Algorithm)

Tests whether a histogram region is connected to any of the current
modes using Monte Carlo algorithm.

## Usage

``` r
connected_to_modes(candidate_mode, curr_mode, g, cluster, ci, L, B)
```

## Arguments

- candidate_mode:

  A positive integer giving the row index of the candidate region.

- curr_mode:

  An integer vector of row indices identifying the currently identified
  mode regions in the histogram.

- g:

  An `igraph` undirected graph object representing the adjacency graph
  of the histogram regions.

- cluster:

  A list of length `length(curr_mode)`, where the \\k\\-th element is a
  vector indicating regions in the neighborhood of mode \\k\\, as
  returned by
  [`find_neighbors`](https://zq00.github.io/BetaTree/reference/find_neighbors.md).

- ci:

  A numeric matrix with two columns giving the lower and upper
  confidence bounds for the average density of each histogram region.

- L:

  A positive integer. The length of random path.

- B:

  A positive integer. The number of independent random paths used to
  assess whether a candidate region is connected to any of the current
  modes.

## Value

A named list with one element:

- `connected`:

  A character string: `"connected"` if the candidate is *not* a distinct
  mode and `"unconnected"` otherwise.

## Details

Generate `B` random paths of length `L` starting from `candidate_mode`
and evaluate (1) whether it connects the `candidate_mode` to any
`curr_mode`; if so, (2) whether there exists a region along the path
whose upper confidence bound is below the lower confidence bounds of
both the `candidate_mode` and the `curr_mode`. If such a region does not
exist, then the `candidate_mode` *cannot* be a distinct mode.

## See also

[`FindModesApproximate`](https://zq00.github.io/BetaTree/reference/FindModesApproximate.md)
calls this function to test whether a candidate region is a distince
mode.
