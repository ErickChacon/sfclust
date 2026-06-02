# Plot log marginal likelihood convergence trace

Plots the log marginal likelihood across MCMC samples for an `sfclust`
object, highlighting the selected sample.

## Usage

``` r
plot_log_mlik(x, sample = x$clust$id, ...)
```

## Arguments

- x:

  An `sfclust` object.

- sample:

  Integer specifying the clustering sample to highlight. Defaults to the
  last sample.

- ...:

  Additional arguments passed to `geom_line()`.

## Value

A `ggplot2` object.
