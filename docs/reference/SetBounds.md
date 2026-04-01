# Add confidence bounds to bounded regions in the k-d tree

Recursively set confidence bounds for bounded nodes in the k-d tree.

## Usage

``` r
SetBounds(node, ahat, n)
```

## Arguments

- node:

  A list storing information about a region.

- ahat:

  A vector of significance level at each depth, starting at depth = 0.
  Output of the funtion `ConfLevel`.

- n:

  Total number of observations.

## Value

An updated k-d tree starting at input node with `lower` and `upper`
filled in.

## Details

Because \$\$ F(R_k)\sim \mathrm{Binomial}(n_k+1,n-n_k), \$\$ where
\\R_k\\ is a region, \\n\\ is the total number of obs. and \\n_k\\ is
the number of obs. in \\R_k\\, the upper and lower confidence bounds of
the average density in a bounded region \\R_k\\ at depth \\D\\ is given
by \$\$ \text{lower}(R_k) = \frac{qBeta(\hat{\alpha}\_D/2, n_k+1,
n-n_k)}{\|R_k\|}\\ \text{upper}(R_k) = \frac{qBeta(1-\hat{\alpha}\_D/2,
n_k+1, n-n_k)}{\|R_k\|} \$\$ where \\\hat{\alpha}\_D\\ is the
significance levels at depth \\D\\, and \\\|R_k\|\\ is the volume of the
region.

We select nodes in which data are approximately uniformly distributed,
i.e., for which the empirical density falls within the confidence
intervals of all of its children. We set the lower and upper bound of
each region as \$\$ \widetilde{\text{lower}} = \max\left(\text{lower},
\text{lower}(\text{1st child}), \text{lower}(\text{2nd child})\right) \\
\widetilde{\text{upper}} = \min\left(\text{upper},
\text{upper}(\text{1st child}), \text{upper}(\text{2nd child})\right),
\$\$ so that when we select nodes, we pick the largest nodes such that
the empirical density is between \\\widetilde{\text{lower}}\\ and
\\\widetilde{\text{upper}}\\. The significance levels are computed with
the function ConfLevel.
