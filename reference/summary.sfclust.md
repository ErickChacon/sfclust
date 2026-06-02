# Summary method for sfclust objects

This function summarizes the cluster assignments from the desired
clustering `sample`.

## Usage

``` r
# S3 method for class 'sfclust'
summary(object, sample = object$clust$id, sort = FALSE, ...)
```

## Arguments

- object:

  An object of class 'sfclust'.

- sample:

  An integer specifying the clustering sample number to be summarized
  (default is the last sample).

- sort:

  Logical value indicating if clusters should be relabel based on number
  of elements.

- ...:

  Additional arguments passed to `print.default`.

## Value

Invisibly returns a table with the number of regions in each cluster for
the selected sample. The function also prints a summary that includes:

- the within-cluster model formula,

- the total number of MCMC clustering samples,

- the cluster membership counts for the specified sample (optionally
  sorted),

- and the log marginal likelihood of the selected clustering sample.
