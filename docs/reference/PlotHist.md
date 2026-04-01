# Plot a 2D Beta tree histogram

Draw a Beta tree histogram.

## Usage

``` r
PlotHist(X, B, nsample = 10000, ...)
```

## Arguments

- X:

  data matrix of size n by 2.

- B:

  matrix of rectangles returned by the function
  [SelectNodes](SelectNodes.md) or [BuildHist](BuildHist.md).

- nsample:

  maximum number of data points to display in the plot.

- ...:

  additional parameter for the graph: `low`, `mid`, `high` specifies the
  colors in
  [scale_fill_gradient2](https://ggplot2.tidyverse.org/reference/scale_gradient.html)
  function.

## Value

A `ggplot` object.

## Details

`PlotHist` plots the histogram superimposed on the original dataset. If
the number of data points exceeds `nsample`, randomly select `nsample`
obs. to show in the graph. Default value is \\nsample=10^4\\

## Examples

``` r
X <-  matrix(rnorm(2000, 0, 1), ncol = 2)
rect <- BuildHist(X, plot = FALSE)
PlotHist(X, rect)
```
