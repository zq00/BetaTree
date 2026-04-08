# Select nodes in the histogram

Extract maximal (w.r.t. inclusion) nodes that are bounded and whose
empirical density is within confidence bounds.

## Usage

``` r
SelectNodes(node, B, ahat, n)
```

## Arguments

- node:

  A k-d tree output from `SetBounds`.

- B:

  A matrix containing nodes that have been selected.

- ahat:

  A vector of significance levels at each depth, starting at depth = 0.
  Output of the funtion `ConfLevel`

- n:

  Total number of observations.

## Value

A matrix of selected regions within and including the input node. Each
row represents one selected region and has \\2d + 5\\ columns,
containing the lower and upper bounds of the region, the empirical
density, lower and upper confidence bounds of the average density in
this region, number of obs. inside and the depth of the node.

## Details

A bounded node is selected if the empirical density is inside the
confidence bounds of all of its children. The empirical density is
defined as \$\$ h = \frac{n+1}{N \|R\|}, \$\$ where \\N\\ is the total
number of obs., \\n\\ is the number of obs. in the region, and \\\|R\|\\
is the volume of the region. If the node is not selected and is not a
leaf, then move on to the left and right children of the node.
