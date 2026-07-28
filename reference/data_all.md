# Prepare data in long format

Convert spatio-temporal data to long format with a flat spatial index
(`ids`) and ordered observation indices (`id_<dimname>` for each
non-spatial dimension). This is a pure converter: all rows are returned
including cells that are all-NA. No filtering is applied. Use
[`filter_df()`](filter_df.md) downstream to restrict to valid cells.

## Usage

``` r
data_all(x, spnames = "geometry")
```

## Arguments

- x:

  A `stars` object containing the spatial and functional dimensions.

- spnames:

  Character vector with the names of the spatial dimensions of `x`. Use
  a single name (e.g. `"geometry"`) for vector geometry data, or two
  names (e.g. `c("x", "y")`) for raster data. Functional dimensions are
  derived automatically as all remaining dimensions.

## Value

A long-format data frame with columns `id` (flat array position in the
original stars object, column-major over all dimensions), `ids` (flat
spatial index, column-major over `spnames`, `spnames[1]` varies
fastest), `id_<dimname>` for each functional dimension, and all
variables from `x`. All rows are returned; `ids` is not remapped to
`1..n_valid`.

## Examples

``` r

library(sfclust)
library(stars)
#> Loading required package: abind
#> Loading required package: sf
#> Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.4.0; sf_use_s2() is TRUE

dims <- st_dimensions(
  geometry = st_sfc(lapply(1:5, function(i) st_point(c(i, i)))),
  time = seq(as.Date("2024-01-01"), by = "1 day", length.out = 3)
)
stdata <- st_as_stars(cases = array(1:15, dim = c(5, 3)), dimensions = dims)

data_all(stdata)
#>    id ids id_time       time cases
#> 1   1   1       1 2024-01-01     1
#> 2   2   2       1 2024-01-01     2
#> 3   3   3       1 2024-01-01     3
#> 4   4   4       1 2024-01-01     4
#> 5   5   5       1 2024-01-01     5
#> 6   6   1       2 2024-01-02     6
#> 7   7   2       2 2024-01-02     7
#> 8   8   3       2 2024-01-02     8
#> 9   9   4       2 2024-01-02     9
#> 10 10   5       2 2024-01-02    10
#> 11 11   1       3 2024-01-03    11
#> 12 12   2       3 2024-01-03    12
#> 13 13   3       3 2024-01-03    13
#> 14 14   4       3 2024-01-03    14
#> 15 15   5       3 2024-01-03    15
```
