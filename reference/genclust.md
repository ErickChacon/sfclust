# Generate initial cluster assignments

Creates an undirected graph from spatial data or a weighted adjacency
matrix, computes its minimum spanning tree (MST), and partitions it into
`nclust` clusters. Accepts weighted `matrix` or `Matrix` objects, and
`stars` objects with either vector geometry (contiguity via
[`sf::st_touches()`](https://r-spatial.github.io/sf/reference/geos_binary_pred.html))
or raster dimensions (4-connected grid adjacency for non-NA pixels).

## Usage

``` r
genclust(x, ...)

# Default S3 method
genclust(x, ...)

# S3 method for class 'matrix'
genclust(x, nclust = 10, weights = NULL, ...)

# S3 method for class 'Matrix'
genclust(x, nclust = 10, weights = NULL, ...)

# S3 method for class 'stars'
genclust(x, nclust = 10, spnames = NULL, response = NULL, weights = NULL, ...)
```

## Arguments

- x:

  Spatial data or adjacency matrix. Accepted classes: `matrix`,
  `Matrix`, or `stars`.

- ...:

  Not used; required for S3 method consistency.

- nclust:

  Integer. Number of initial clusters (default `10`).

- weights:

  Optional numeric vector of edge weights with `n^2` elements, where `n`
  is the total number of spatial units (before filtering). If `NULL`,
  random weights are assigned.

- spnames:

  Character vector with the names of the spatial dimensions. If `NULL`,
  auto-detected as dimensions with a non-`NA` regular `delta` in
  [`stars::st_dimensions()`](https://r-spatial.github.io/stars/reference/st_dimensions.html).
  Currently supports 1D (vector geometry) and 2D raster grids; 3D
  spatial grids (e.g. `x`, `y`, `z`) are not yet supported.

- response:

  Character. Name of the attribute in `x` to use for determining valid
  spatial cells (cells where all observations are NA are excluded). If
  `NULL` (default), all spatial cells are treated as valid.

## Value

A list with:

- `graph`: undirected igraph object representing spatial contiguity.

- `mst`: minimum spanning tree of `graph`.

- `membership`: integer vector of cluster assignments (length = number
  of valid spatial units).

- `valid_ids`: integer vector of flat spatial positions included in the
  graph (only present for `stars` input; `NULL` for matrix/Matrix
  input).

## Examples

``` r

library(sfclust)
library(sf)
library(stars)

# stars object with vector geometry
geom <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))
x <- st_as_stars(st_sf(z = 1:6, geometry = geom))
clust <- genclust(x, nclust = 3)
plot(st_sf(geom, cluster = factor(clust$membership)))


# stars raster input
x <- st_as_stars(cluster = matrix(1:35, 5))
clust <- genclust(x, nclust = 4)
x$cluster <- clust$membership
plot(x, col = rainbow(4))


# matrix input
x <- matrix(c(0,1,0,1, 1,0,1,0, 0,1,0,1, 1,0,1,0), nrow = 4)
clust <- genclust(x, nclust = 2)
clust$membership
#> [1] 1 2 1 1
```
