# Spatio-temporal NDWI2 raster data from El Chaparrillo

A `stars` raster object containing the Normalized Difference Water Index
2 (NDWI2) derived from Sentinel-2 satellite imagery over El Chaparrillo,
an agricultural research center in Spain, across 86 time points.

## Usage

``` r
data(chapa)
```

## Format

A `stars` object with:

- ndwi2:

  NDWI2 values (numeric) derived from Sentinel-2 bands. Non-NA values
  define the study area mask.

- dimensions:

  Three dimensions: `x` (26 columns), `y` (27 rows), `time` (86 dates
  from 2025-05-15 to 2025-12-01). CRS: WGS 84 / UTM zone 30N.

## Examples

``` r

library(sfclust)

data(chapa)
chapa
#> stars object with 3 dimensions and 1 attribute
#> attribute(s):
#>               Min.   1st Qu.    Median     Mean   3rd Qu.      Max.   NAs
#> ndwi2  -0.03494883 0.1527428 0.1872349 0.193792 0.2270235 0.5227525 31562
#> dimension(s):
#>      from to  offset delta                refsys point
#> x       1 26  418080    10 WGS 84 / UTM zone 30N FALSE
#> y       1 27 4320790   -10 WGS 84 / UTM zone 30N FALSE
#> time    1 86      NA    NA                  Date    NA
#>                         values x/y
#> x                         NULL [x]
#> y                         NULL [y]
#> time 2025-05-15,...,2025-12-01    
plot(chapa["ndwi2"])

```
