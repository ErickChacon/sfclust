# Plot function for `sfclust` objects

This function visualizes the estimated clusters from an `sfclust`
object. It can display: (1) a map of regions colored by their assigned
cluster, (2) the functional shapes of the linear predictors for each
cluster, and (3) a traceplot of the log marginal likelihood. A
conditional legend is added if the number of clusters is less than 10.

## Usage

``` r
# S3 method for class 'sfclust'
plot(
  x,
  sample = x$clust$id,
  which = 1:3,
  clusters = NULL,
  sort = FALSE,
  legend = FALSE,
  geom_before = NULL,
  ...
)
```

## Arguments

- x:

  An `sfclust` object containing the clustering results, including the
  cluster assignments and model parameters.

- sample:

  Integer specifying the clustering sample number to summarize. Defaults
  to the last sample.

- which:

  Integer vector indicating which plot to display. Options are: - 1: Map
  of regions colored by cluster assignment. - 2: Functional shapes of
  the linear predictors for each cluster. - 3: Traceplot of the log
  marginal likelihood.

- clusters:

  Optional vector specifying which clusters to plot. If `NULL`, all
  clusters are included.

- sort:

  Logical value indicating whether clusters should be relabeled based on
  the number of elements. Default is `FALSE`.

- legend:

  Logical value indicating whether a legend should be included in the
  plot. Default is `FALSE`.

- geom_before:

  An optional `ggplot2` geom layer to add before the cluster map layer
  (e.g., a background tile). Default is `NULL`.

- ...:

  Additional arguments passed to the underlying plotting functions.

## Value

A composed `patchwork` object displaying the selected subgraphs as
specified by `which`.
