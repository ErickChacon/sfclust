# Prepare data for a cluster

Subset a spatio-temporal dataset for a cluster and convert it to a long
format with indices for time and spatial location.

## Usage

``` r
data_each(k, membership, stdata, stnames = c("geometry", "time"))
```

## Arguments

- k:

  The cluster number to subset.

- membership:

  A vector defining the cluster membership for each region.

- stdata:

  A stars object containing spatial-temporal dimensions defined in
  `stnames`.

- stnames:

  The names of the `spatial` and `temporal` dimensions.

## Value

A long-format data frame with ids for each observation and for spatial
and time indexing.

## Examples

``` r

library(sfclust)
library(stars)

dims <- st_dimensions(
  geometry = st_sfc(lapply(1:5, function(i) st_point(c(i, i)))),
  time = seq(as.Date("2024-01-01"), by = "1 day", length.out = 3)
)
stdata <- st_as_stars(cases = array(1:15, dim = c(5, 3)), dimensions = dims)

data_each(k = 2, membership = c(1, 1, 1, 2, 2), stdata)
#>   id ids idt       time cases
#> 1  4   4   1 2024-01-01     4
#> 2  5   5   1 2024-01-01     5
#> 3  9   4   2 2024-01-02     9
#> 4 10   5   2 2024-01-02    10
#> 5 14   4   3 2024-01-03    14
#> 6 15   5   3 2024-01-03    15
```
