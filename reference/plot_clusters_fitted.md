# Plot functional shapes of cluster linear predictors

Plots the estimated mean functional shape (linear predictor or
inverse-link scale) for each cluster in a given MCMC sample of an
`sfclust` object.

## Usage

``` r
plot_clusters_fitted(
  x,
  sample = x$clust$id,
  clusters = NULL,
  sort = FALSE,
  legend = FALSE,
  inv_link = TRUE,
  fnames = NULL,
  ...
)
```

## Arguments

- x:

  An `sfclust` object.

- sample:

  Integer specifying the clustering sample to display. Defaults to the
  last sample.

- clusters:

  Optional vector of cluster IDs to include. If `NULL`, all clusters are
  shown.

- sort:

  Logical; if `TRUE`, clusters are relabeled by decreasing size. Default
  is `FALSE`.

- legend:

  Logical; if `TRUE`, a color legend is included. Default is `FALSE`.

- inv_link:

  Logical; if `TRUE` (default), values are shown on the inverse-link
  (mean) scale.

- fnames:

  Character. Name of the column to use as the x-axis (functional
  dimension). If `NULL`, taken from the result's stored args (set
  automatically by [`sfclust()`](sfclust.md)).

- ...:

  Additional arguments passed to `geom_line()`.

## Value

A `ggplot2` object.
