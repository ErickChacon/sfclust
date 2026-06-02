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
  level.

- ...:

  Additional arguments, currently not used.

## Value

A `stars` object with linear predictor fitted values at regions levels.
In case `aggregate = TRUE`, the `output` produces an `stars` objecto at
cluster levels.

## Details

The function first checks if the provided `sample` value is valid (i.e.,
it is within the range of available clustering samples). If the
specified `sample` does not match the current clustering `id`, the
`sfclust` object is updated accordingly. It then retrieves the
membership assignments and cluster models for the selected sample,
calculates the linear predictions for each cluster, and combines them
into a matrix of fitted values.

## Examples

``` r

# \donttest{
library(sfclust)

data(stgaus)
result <- sfclust(stgaus, formula = y ~ f(idt, model = "rw1"), niter = 10,
  nmessage = 1)
#> Iteration 1: clusters = 10, births = 0, deaths = 0, changes = 0, hypers = 0, log_mlike = -1867.54278987303
#> Iteration 2: clusters = 11, births = 1, deaths = 0, changes = 0, hypers = 0, log_mlike = -1818.05514376303
#> Iteration 3: clusters = 11, births = 1, deaths = 0, changes = 0, hypers = 0, log_mlike = -1818.05514376303
#> Iteration 4: clusters = 12, births = 2, deaths = 0, changes = 0, hypers = 0, log_mlike = -1503.31662432528
#> Iteration 5: clusters = 12, births = 2, deaths = 0, changes = 0, hypers = 0, log_mlike = -1503.31662432528
#> Iteration 6: clusters = 12, births = 2, deaths = 0, changes = 0, hypers = 0, log_mlike = -1503.31662432528
#> Iteration 7: clusters = 12, births = 2, deaths = 0, changes = 0, hypers = 0, log_mlike = -1503.31662432528
#> Iteration 8: clusters = 13, births = 3, deaths = 0, changes = 0, hypers = 0, log_mlike = -840.715686749306
#> Iteration 9: clusters = 14, births = 4, deaths = 0, changes = 0, hypers = 0, log_mlike = -426.665204489668
#> Iteration 10: clusters = 14, births = 4, deaths = 0, changes = 0, hypers = 0, log_mlike = -426.665204489668

# Estimated values ordering clusters by size
df_est <- fitted(result, sort = TRUE)

# Estimated values aggregated by cluster
df_est <- fitted(result, aggregate = TRUE)

# Estimated values using a particular clustering sample
df_est <- fitted(result, sample = 3)
# }
```
