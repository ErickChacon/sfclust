# Plot function for `sfclust` objects

Visualizes fitted cluster functions (plot 1) and a log marginal
likelihood traceplot (plot 2). For `sfclust_stars` objects, also
includes a spatial map of cluster assignments (plot 1 in that method,
shifting the others to 2 and 3).

## Usage

``` r
# S3 method for class 'sfclust'
plot(
  x,
  sample = x$clust$id,
  which = 1:2,
  clusters = NULL,
  sort = FALSE,
  legend = FALSE,
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

- which:

  Integer vector indicating which plots to show. For `sfclust`: 1 =
  cluster functions, 2 = log marginal likelihood. For `sfclust_stars`: 1
  = map, 2 = cluster functions, 3 = log marginal likelihood.

- clusters:

  Optional vector of cluster IDs to include. If `NULL`, all clusters are
  shown.

- sort:

  Logical; if `TRUE`, clusters are relabeled by decreasing size. Default
  is `FALSE`.

- legend:

  Logical; if `TRUE`, a legend is included. Default is `FALSE`.

- fnames:

  Character. Column name for the x-axis of cluster function plots. If
  `NULL`, taken from the result's stored args.

- ...:

  Additional arguments passed to underlying plot functions.

## Value

A composed `patchwork` object.
