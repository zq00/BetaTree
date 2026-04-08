# Build an Adaptive Beta-Tree Histogram

Constructs an Beta-tree histogram there partition coordinates are chosen
adaptively based on tests of marginal uniformity and pairwise
independence.

## Usage

``` r
build_adaptive_histogram(
  X,
  alpha = 0.1,
  thresh_marginal = alpha/ncol(X),
  thresh_interaction = stats::qgamma(shape = ncol(X) - 1, rate = 1, p = 1 - alpha),
  method = "weighted_bonferroni"
)
```

## Arguments

- X:

  A numeric matrix dimension \\n \times d\\ containing observations
  falling inside the current region. Each row corresponds to one
  observation.

- alpha:

  Significance level, default is `alpha = 0.1`. The Beta-tree histogram
  provides simultaneous confidence intervals for average densities in
  each histogram region with coverage probability \\1-\alpha\\.

- thresh_marginal:

  A numeric value between \\(0, 1)\\. The significance threshold for the
  marginal Anderson-Darling uniformity test. Default value is
  `alpha / d`.

- thresh_interaction:

  A positive numeric value. The threshold on the combined interaction
  score \\\sum\_{j\neq i} -\log p\_{ij}\\ for coordinate \\i\\, where
  \\p\_{ij}\\ is the p-value of testing independence between \\X_i\\ and
  \\X_j\\. Default value is \\(1-\alpha)\\ quantile of a Gamma variable
  with shape \\d-1\\ and rate 1.

- method:

  A character string specifying whether to use `method = "bonferroni"`
  or `method = "weighted_bonferroni"` method to adjust for multiple
  hypothesis testing when computing significance levels.

## Value

A named list with three components:

- `split_dir`:

  An integer vector recording the splitting coordinate chosen at each
  step. The first element is always `NA`.

- `hist`:

  A numeric matrix of dimension \\m \times (2d + 5)\\. Each row
  represents one region in the Beta-tree histogram. The first \\d\\
  columns are the lower bounds for the region; the next \\d\\ columns
  are the upper bounds. The final 5 columns give the empirical density,
  lower and upper confidence bounds of average density, , upper
  confidence bound, number of observations inside, and tree depth of
  each region. See also
  [`BuildHist`](https://zq00.github.io/BetaTree/reference/BuildHist.md).

- `lower`:

  A numeric vector of length \\d\\ corresponding to the the minimum
  observation in each coordinate.

- `upper`:

  A numeric vector of length \\d\\ corresponding to the the maximum
  observation in each coordinate.

## Details

The adaptive Beta-tree histogram uses the same algorithm as the
Beta-tree histogram (see
[`BuildHist`](https://zq00.github.io/BetaTree/reference/BuildHist.md)
function and package Vignettes `beta_trees`) with two differences:

First, the adaptive Beta-tree automatically computes a bounded histogram
by placing the first splits at the minimum and maximum of each
coordinates.

Second, the adaptive Beta-tree uses
[`get_partition_dim`](https://zq00.github.io/BetaTree/reference/get_partition_dim.md)
to select partition dimension at each step. The same set of thresholds
`thresh_marginal` and `thresh_interaction` are used at every step.

## See also

[`get_partition_dim`](https://zq00.github.io/BetaTree/reference/get_partition_dim.md)
adaptively selectx partition coordinates.
[`BuildHist`](https://zq00.github.io/BetaTree/reference/BuildHist.md)
builds a Beta-tree histogram.

## Examples

``` r
set.seed(42)

# Bivariate standard normal data
X <- matrix(rnorm(1000), nrow = 500, ncol = 2)

adaptive_beta_tree <- build_adaptive_histogram(
  X                 = X,
  alpha             = 0.1,
  method            = "weighted_bonferroni"
)

# Number of histogram regions
nrow(adaptive_beta_tree$hist)
#> [1] 19

# Splitting directions
table(adaptive_beta_tree$split_dir)
#> 
#>  1  2 
#> 19 12 
```
