# Choose Partition Dimension

Selects partition dimension by testing marginal uniformity and pairwise
independence.

## Usage

``` r
get_partition_dim(x, thresh_marginal, thresh_interaction, lower, upper)
```

## Arguments

- x:

  A numeric matrix dimension \\n \times d\\ containing observations
  falling inside the current region. Each row corresponds to one
  observation.

- thresh_marginal:

  A numeric value between \\(0, 1)\\. The significance threshold for the
  marginal Anderson-Darling uniformity test.

- thresh_interaction:

  A positive numeric value. The threshold on the combined interaction
  score \\\sum\_{j\neq i} -\log p\_{ij}\\ for coordinate \\i\\, where
  \\p\_{ij}\\ is the p-value of testing independence between \\X_i\\ and
  \\X_j\\

- lower:

  A numeric vector of length \\d\\ specifying the lower bounds of the
  current region.

- upper:

  A numeric vector of length \\d\\ specifying the upper bounds of the
  current region.

## Value

A named list with three elements:

- `below_marginal`:

  A vector of integers corresponding to coordinates whose marginal
  Anderson-Darling p-value is below `thresh_marginal`, `NA` if no
  rejection.

- `below_interaction`:

  A vector of integers corresponding to coordinates whose combined
  interaction score \\p_i\\ exceeds `thresh_interaction`, `NA` if no
  rejection.

- `p`:

  Selected partition coordinate.

## Details

Adaptively choosing partition dimension by first testing marginal
uniformity using Anderson-Darling test for each coordinate \\i \in \\1,
\ldots, d\\\\, i.e. we test if \\x\_{1i}, \ldots, x\_{ni}\\ are
uniformly distributed within the lower and upper bounds of coordinate
\\i\\. If the minimum p-value is less than `thresh_marginal`, randomly
sample one coordinate with p-value less than `thresh_marginal` as the
partition dimension.

Otherwise, test pairwise independence between all \\\binom{d}{2}\\
coordinate pairs \\(i, j)\\ using a Fisher's exact test. To test if
\\X_i\\ is independent of \\X_j\\ in the region \\R\\, we divide \\R\\
into four quadrants along the middle of \\i\\ and \\j\\ coordinates and
apply Fisher's exact test to counts of observations in each quadrant.
Compute a combined p-value for \\i\\th coordinate as \\p_i = \sum\_{j
\neq i} -\log p\_{ij}\\ and randomly sample one coordinate with \\p_i \>
\\`thresh_interaction` as partition dimension.

If no coordinate is significant in either tests, randomly pick one
coordinate to partition.

## See also

[`build_adaptive_histogram`](https://zq00.github.io/BetaTree/reference/build_adaptive_histogram.md)
for the function that calls `get_partition_dim` to adaptively choose
partition locations.

## Examples

``` r
set.seed(42)

# Uniform observations
x <- matrix(runif(300), nrow = 100, ncol = 3)
get_partition_dim(x,
                 thresh_marginal    = 0.05,
                 thresh_interaction = qgamma(p=1-0.05, shape=3-1, rate=1),
                 lower = rep(0, 3),
                 upper = rep(1, 3))
#> $below_marginal
#> [1] 3
#> 
#> $below_interaction
#> [1] NA
#> 
#> $p
#> [1] 1
#> 

# Non-uniform observations
y <- matrix(runif(300), nrow = 100, ncol = 3)
y[, 1] <- rbeta(100, 3, 3)
get_partition_dim(y,
                 thresh_marginal    = 0.05,
                 thresh_interaction = qgamma(p=1-0.05, shape=3-1, rate=1),
                 lower = rep(0, 3),
                 upper = rep(1, 3))
#> $below_marginal
#> [1] 1
#> 
#> $below_interaction
#> [1] NA
#> 
#> $p
#> [1] 1
#> 
```
