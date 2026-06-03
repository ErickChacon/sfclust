# Generate clusters for spatial clustering

Creates an undirected graph from spatial polygonal data, computes its
minimum spanning tree (MST), and generates `nclust` clusters. This
function is used to initialize cluster membership in a clustering
algorithm, such as `sfclust`.

## Usage

``` r
genclust(x, nclust = 10, weights = NULL)
```

## Arguments

- x:

  An `sf` or `sfc` object representing spatial polygonal data. It can
  also be a `matrix` or `Matrix` object with non-zero values
  representing weighted connectivity between units.

- nclust:

  Integer, specifying the initial number of clusters.

- weights:

  Optional `numeric` vector or `matrix` of weights between units in `x`.
  It should have dimensions `n^2`, where `n` is the number of units in
  `x`. If NULL, random weights are assigned.

## Value

A list with three elements:

- `graph`: The undirected graph object representing spatial contiguity.

- `mst`: The minimum spanning tree.

- `membership`: The cluster membership for elements in `x`.

## Examples

``` r

library(sfclust)
library(sf)

x <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))

# using distance between geometries
clust <- genclust(x, nclust = 3, weights = st_distance(st_centroid(x)))
print(clust)
#> $graph
#> IGRAPH 2e6aeab UNW- 6 11 -- 
#> + attr: name (v/c), weight (e/n)
#> + edges from 2e6aeab (vertex names):
#>  [1] 1--2 1--4 1--5 2--3 2--4 2--5 2--6 3--5 3--6 4--5 5--6
#> 
#> $mst
#> IGRAPH 274750c UNW- 6 5 -- 
#> + attr: name (v/c), vid (v/n), weight (e/n)
#> + edges from 274750c (vertex names):
#> [1] 1--4 2--3 3--6 4--5 5--6
#> 
#> $membership
#> 1 2 3 4 5 6 
#> 1 2 3 3 3 3 
#> 
plot(st_sf(x, cluster = factor(clust$membership)))


# using increasing weights
cluster_ini <- genclust(x, nclust = 3, weights = 1:36)
print(cluster_ini)
#> $graph
#> IGRAPH f61ed68 U-W- 6 11 -- 
#> + attr: weight (e/n)
#> + edges from f61ed68:
#>  [1] 1--2 1--4 1--5 2--3 2--4 2--5 2--6 3--5 3--6 4--5 5--6
#> 
#> $mst
#> IGRAPH 4ca6b67 U-W- 6 5 -- 
#> + attr: vid (v/n), weight (e/n)
#> + edges from 4ca6b67:
#> [1] 1--2 1--4 1--5 2--3 2--6
#> 
#> $membership
#> [1] 1 1 1 1 2 3
#> 
plot(st_sf(x, cluster = factor(cluster_ini$membership)))


# using on random weights
cluster_ini <- genclust(x, nclust = 3, weights = runif(36))
print(cluster_ini)
#> $graph
#> IGRAPH d6deba4 U-W- 6 11 -- 
#> + attr: weight (e/n)
#> + edges from d6deba4:
#>  [1] 1--2 1--4 1--5 2--3 2--4 2--5 2--6 3--5 3--6 4--5 5--6
#> 
#> $mst
#> IGRAPH b4c2201 U-W- 6 5 -- 
#> + attr: vid (v/n), weight (e/n)
#> + edges from b4c2201:
#> [1] 1--4 2--4 2--6 3--5 3--6
#> 
#> $membership
#> [1] 1 2 3 1 3 3
#> 
plot(st_sf(x, cluster = factor(cluster_ini$membership)))

```
