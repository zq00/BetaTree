# Plot a multivariate histogram in two dimensions

Visualize an adaptive beta-tree histogram in two-dimensions marginally
or conditional on inside a hyperrectangle.

## Usage

``` r
plot_histogram(plot_data, show_data = TRUE)
```

## Arguments

- plot_data:

  A named list as returned by
  [`plot_histogram_data`](https://zq00.github.io/BetaTree/reference/plot_histogram_data.md).

- show_data:

  Logical. If `TRUE` (default), shows the data points in
  `plot_data$data`.

## Value

A [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.

## Details

[`plot_histogram_data`](https://zq00.github.io/BetaTree/reference/plot_histogram_data.md)
computes the probability content in each region, `plot_histogram`
further computes the conditional distribution after dividing by the
total probability content.

## See also

[`plot_histogram_data`](https://zq00.github.io/BetaTree/reference/plot_histogram_data.md)
for computing the input `plot_data`,
[`build_adaptive_histogram`](https://zq00.github.io/BetaTree/reference/build_adaptive_histogram.md)
for producing the histogram.

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 10000
p <- 4
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
adaptive_beta_tree <- build_adaptive_histogram(X, alpha = 0.1)

plot_dat <- plot_histogram_data(
  X          = X,
  hist       = adaptive_beta_tree$hist,
  plot_coord = c(1, 2),
  nint       = 20,
  show_data  = TRUE,
  ndat       = 500
)

plot_histogram(plot_dat)
} # }
```
