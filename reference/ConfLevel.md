# Set significance level for regions at each tree depth

Computes the significance level at each tree depth.

## Usage

``` r
ConfLevel(nd, alpha, method)
```

## Arguments

- nd:

  A vector. Number of bounded regions at every depth (starting at depth
  0).

- alpha:

  Numeric. The significance level (of the average density) of all the
  regions sums to \\\alpha\\.

- method:

  Use `method = "bonferroni"` or `method = "weighted_bonferroni"` method
  to adjust for multiple hypothesis testing.

## Value

A vector of significance level \\\alpha_D\\ at every tree depth
(starting at depth = 0)

## Details

The root node automatically has significance level 0. If
`method = "bonferroni"`, the significance level at depth d is \$\$
\hat{\alpha}\_d = \frac{\alpha}{(d\_\max - d\_\min + 1) \* n_d}, d\geq
d\_\min, \$\$ where \\d\_\max\\ is the maximum tree depth, \\d\_\min\\
is the minimum depth (greater than 0) where there are bounded regions,
and \\n_d\\ is number of bounded regions at depth \\d\\.

If `method = weighted_bonferroni`, the confidence level at depth d is
\$\$ \hat{\alpha}\_d = \frac{\alpha}{(d\_\max - d + 2) n_d
\sum\_{B=2}^{d\_\max - d\_\min + 2}1/B},d\geq d\_\min. \$\$
