# Plot a spatial map of cluster assignments

Produces a `ggplot2` map of spatial regions colored by their cluster
assignment for a given MCMC sample of an `sfclust` object.

## Usage

``` r
plot_clusters_map(
  x,
  sample = x$clust$id,
  clusters = NULL,
  sort = FALSE,
  legend = FALSE,
  geom_before = NULL,
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

  Logical; if `TRUE`, a fill legend is included. Default is `FALSE`.

- geom_before:

  An optional `ggplot2` geom layer to add before the cluster fill layer.
  Default is `NULL`.

- ...:

  Additional arguments passed to `geom_sf()`.

## Value

A `ggplot2` object.
