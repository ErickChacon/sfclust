# Fit models and compute the log marginal likelihood for all clusters

Fit the specified INLA model to each cluster and compute the log
marginal likelihood for each cluster specified in the membership vector.

## Usage

``` r
log_mlik_all(
  membership,
  stdata,
  stnames = c("geometry", "time"),
  correction = TRUE,
  detailed = FALSE,
  ...
)
```

## Arguments

- membership:

  Integer, character or factor vector indicating the cluster membership
  for each spatial unit.

- stdata:

  A stars object with spatial-temporal dimensions defined in `stnames`,
  and including predictors and response variables.

- stnames:

  The names of the `spatial` and `temporal` dimensions of the `stdata`
  object.

- correction:

  Logical value indicating whether a correction for dispersion.

- detailed:

  Logical value indicating whether to return the INLA model instead of
  the log marginal likelihood. The argument `correction` is not applied
  in this case.

- ...:

  Arguments passed to the `inla` function (eg. `family`, `formula` and
  `E`).

## Value

A numeric vector containing the log marginal likelihood for each cluster
or the the fitted INLA model for each cluster when `detailed = TRUE`.

## Examples

``` r


# \donttest{
library(sfclust)
library(stars)

dims <- st_dimensions(
  geometry = st_sfc(lapply(1:5, function(i) st_point(c(i, i)))),
  time = seq(as.Date("2024-01-01"), by = "1 day", length.out = 3)
)
stdata <- st_as_stars(
  cases = array(rpois(15, 100 * exp(1)), dim = c(5, 3)),
  temperature = array(runif(15, 15, 20), dim = c(5, 3)),
  expected = array(100, dim = c(5, 3)),
  dimensions = dims
)

log_mlik_all(c(1, 1, 1, 2, 2), stdata,
  formula = cases ~ temperature, family = "poisson", E = expected)
#> [1] -46.58291 -33.48540

models = log_mlik_all(c(1, 1, 1, 2, 2), stdata, detailed = TRUE,
  formula = cases ~ temperature, family = "poisson", E = expected)
lapply(models, summary)
#> [[1]]
#> Time used:
#>     Pre = 0.198, Running = 0.138, Post = 0.00989, Total = 0.345 
#> Fixed effects:
#>               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
#> (Intercept)  1.174 0.280      0.625    1.174      1.723  1.174   0
#> temperature -0.009 0.016     -0.040   -0.009      0.022 -0.009   0
#> 
#> Marginal log-Likelihood:  -46.58 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
#> 
#> 
#> [[2]]
#> Time used:
#>     Pre = 0.242, Running = 0.138, Post = 0.0102, Total = 0.389 
#> Fixed effects:
#>               mean    sd 0.025quant 0.5quant 0.975quant   mode kld
#> (Intercept)  1.092 0.404      0.299    1.092      1.884  1.092   0
#> temperature -0.006 0.023     -0.052   -0.006      0.040 -0.006   0
#> 
#> Marginal log-Likelihood:  -33.48 
#>  is computed 
#> Posterior summaries for the linear predictor and the fitted values are computed
#> (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')
#> 
#> 
# }
```
