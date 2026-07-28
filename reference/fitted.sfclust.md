# Fitted Values for `sfclust` Objects

This function calculates the fitted values for a specific clustering
sample in an `sfclust` object, based on the estimated models for each
cluster. The fitted values are computed using the membership assignments
and model parameters associated with the selected clustering sample.

## Usage

``` r
# S3 method for class 'sfclust'
fitted(object, sample = object$clust$id, sort = FALSE, aggregate = FALSE, ...)
```

## Arguments

- object:

  An object of class 'sfclust', containing clustering results and
  models.

- sample:

  An integer specifying the clustering sample number for which the
  fitted values should be computed. The default is the `id` of the
  current clustering. The value must be between 1 and the total number
  of clustering (membership) samples.

- sort:

  Logical value indicating if clusters should be relabel based on number
  of elements.

- aggregate:

  Logical value indicating if fitted values are desired at cluster
  level. Only supported for `sfclust_stars` results.

- ...:

  Additional arguments, currently not used.

## Value

A data frame with fitted values and cluster assignments, keyed by `id`.
For `sfclust_stars` objects, a `stars` object is returned instead.

## Examples

``` r

# \donttest{
if (requireNamespace("INLA", quietly = TRUE)) {
library(sfclust)

data(stgaus)
result <- sfclust(stgaus, formula = y ~ f(id_time, model = "rw1"), niter = 10,
  nmessage = 1)

# Estimated values ordering clusters by size
df_est <- fitted(result, sort = TRUE)

# Estimated values aggregated by cluster
df_est <- fitted(result, aggregate = TRUE)

# Estimated values using a particular clustering sample
df_est <- fitted(result, sample = 3)
}
#> Iteration 1: clusters = 10, births = 0, deaths = 0, changes = 0, hypers = 0, log_mlike = -674.186842125447
#> Iteration 2: clusters = 11, births = 1, deaths = 0, changes = 0, hypers = 0, log_mlike = -512.007403492225
#> Iteration 3: clusters = 11, births = 1, deaths = 0, changes = 0, hypers = 0, log_mlike = -512.007403492225
#> Iteration 4: clusters = 11, births = 1, deaths = 0, changes = 0, hypers = 0, log_mlike = -512.007403492225
#> Iteration 5: clusters = 11, births = 1, deaths = 0, changes = 0, hypers = 0, log_mlike = -512.007403492225
#> Iteration 6: clusters = 12, births = 2, deaths = 0, changes = 0, hypers = 0, log_mlike = -444.59932393423
#> Iteration 7: clusters = 12, births = 2, deaths = 0, changes = 0, hypers = 0, log_mlike = -444.59932393423
#> Iteration 8: clusters = 12, births = 2, deaths = 0, changes = 0, hypers = 0, log_mlike = -444.59932393423
#> Iteration 9: clusters = 12, births = 2, deaths = 0, changes = 1, hypers = 0, log_mlike = -281.109178282004
#> Iteration 10: clusters = 12, births = 2, deaths = 0, changes = 1, hypers = 0, log_mlike = -281.109178282004
# }
```
