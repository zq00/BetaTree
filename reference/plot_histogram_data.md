# Prepare data for visualizing a multivariate histogram in two chosen dimensions

Prepares the data needed to visualise a multivariate adaptive histogram
(as returned by
[`build_adaptive_histogram`](https://zq00.github.io/BetaTree/reference/build_adaptive_histogram.md)
or
[`BuildHist`](https://zq00.github.io/BetaTree/reference/BuildHist.md))
in two chosen dimensions.

## Usage

``` r
plot_histogram_data(
  X,
  hist,
  lower_constraint = NULL,
  upper_constraint = NULL,
  plot_coord,
  nint,
  show_data = TRUE,
  ndat = 2000,
  ...
)
```

## Arguments

- X:

  A numeric matrix of dimension \\n \times d\\, the original data.

- hist:

  A numeric matrix of dimension \\m \times (2d + 5)\\ as returned by
  [`build_adaptive_histogram`](https://zq00.github.io/BetaTree/reference/build_adaptive_histogram.md)
  or
  [`BuildHist`](https://zq00.github.io/BetaTree/reference/BuildHist.md).

- lower_constraint:

  A numeric vector of length \\d\\ specifying the lower bound for each
  coordinate. If provided, computes the conditional distribution of
  \\X\\ within the lower and upper constraints. Otherwise, computes the
  marginal distribution in `plot_coord`.

- upper_constraint:

  A numeric vector of length \\d\\ specifying the upper bound for each
  coordinate. Analogously overridden for `plot_coord`.

- plot_coord:

  An integer vector of length 2 giving the column indices of `X` to use
  as the x- and y-axes of the plot.

- nint:

  A positive integer. The number of intervals along each plotting axis.
  The grid contains `nint^2` cells in total.

- show_data:

  Logical. If `TRUE` (default), a subset of at most `ndat` data points
  inside the constraint region is returned.

- ndat:

  A positive integer. Maximum number of data points to return when
  `show_data = TRUE`. If more than `ndat` observations lie within the
  constraints, a random sample of size `ndat` is drawn without
  replacement. Defaults to `2000`.

- ...:

  Currently unused.

## Value

A named list with the following components:

- `grid_lower`:

  A numeric matrix of dimension \\\texttt{nint}^2 \times d\\ giving the
  lower bounds of each grid cell in all \\d\\ coordinates.

- `grid_upper`:

  A numeric matrix of the same dimension giving the upper bounds.

- `plot_coord`:

  The input coordinates for plotting.

- `est_prob`:

  A numeric matrix of dimension \\\texttt{nint}^2 \times 3\\. The first
  column contains the estimated probability mass in each rectangle; the
  second and third columns contain the lower and upper confidence bounds
  of the probability mass.

- `data`:

  A numeric matrix of the (at most `ndat`) data points within the lower
  and upper constraints in the `plot_coord`. Only present when
  `show_data = TRUE`.

## Details

Aggregates the estimated probability content in a regular grid over the
two coordinates specified for plotting. For example, in
three-dimensional setting, let
`lower_constraint = upper_constraint = NULL`, `plot_coord = c(1,2)` and
`nint=10`, then we divide the first two coordinates into 10 equally
sized intervals, creating a grid with 100 rectangles in \\X_1\\ and
\\X_2\\. Then, calculate the probability content in each rectangle by
aggregating all the histogram regions that intersect them. This
calculates the estimated marginal distribution in \\(X_1, X_2)\\.

If `lower_constraint` and `upper_constraint` are specified, then compute
the aggregated probability content within the 2-d rectangle in first two
dimensions and further within the lower and upper constraints of the
other coordinates. This effectively calculates the conditional
distribution of \\X_1\\ and \\X_2\\ within the constraints.

## See also

[`plot_histogram`](https://zq00.github.io/BetaTree/reference/plot_histogram.md)
for rendering the output of this function.

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
} # }
```
