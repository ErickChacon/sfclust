# Plot observed time series by cluster

Plots individual region time series faceted by cluster, overlaid with
the cluster mean, for a given variable in an `sfclust` object.

## Usage

``` r
plot_clusters_series(x, var, clusters = NULL, sort = FALSE, fnames = NULL, ...)
```

## Arguments

- x:

  An `sfclust` object.

- var:

  An unquoted variable name to plot on the y-axis.

- clusters:

  Optional vector of cluster IDs to include. If `NULL`, all clusters are
  shown.

- sort:

  Logical; if `TRUE`, clusters are relabeled by decreasing size. Default
  is `FALSE`.

- fnames:

  Character. Name of the column to use as the x-axis (functional
  dimension). If `NULL`, taken from the result's stored args (set
  automatically by [`sfclust()`](sfclust.md)).

- ...:

  Additional arguments passed to `geom_line()` for individual region
  series.

## Value

A `ggplot2` object with one facet per cluster.
