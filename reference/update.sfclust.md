# Update MCMC Clustering Procedure

This function continues the MCMC sampling of a `sfclust` object based on
previous results or update the model fitting for a specified sample
clustering if the argument `sample` is provided.

## Usage

``` r
# S3 method for class 'sfclust'
update(
  object,
  niter = NULL,
  burnin = 0,
  thin = 1,
  nmessage = 10,
  sample = NULL,
  path_save = NULL,
  nsave = nmessage,
  ...
)
```

## Arguments

- object:

  A `sfclust` object.

- niter:

  An integer specifying the number of additional MCMC iterations to
  perform. If `NULL` (default), no additional iterations are run and the
  within-cluster INLA models are refitted for the current clustering
  sample.

- burnin:

  An integer specifying the number of burn-in iterations to discard.

- thin:

  An integer specifying the thinning interval for recording results.

- nmessage:

  An integer specifying the number of messages to display during the
  process.

- sample:

  An integer specifying the clustering sample number to be executed. The
  default is the last sample (i.e., `nrow(x$samples$membership)`).

- path_save:

  A character string specifying the file path to save the results. If
  `NULL`, results are not saved.

- nsave:

  An integer specifying how often to save results. Defaults to
  `nmessage`.

- ...:

  Additional arguments (currently not used).

## Value

An updated `sfclust` object with (i) new clustering samples if `sample`
is not specified, or (ii) updated within-cluster model results if
`sample` is given.

## Details

This function takes the last state of the Markov chain from a previous
`sfclust_fit` / `sfclust` execution and uses it as the starting point
for additional MCMC iterations. If `niter = NULL` (the default), no
additional iterations are run; instead, the within-cluster INLA models
are fitted for the clustering given by `sample` (or the current
`clust$id` if `sample` is `NULL`). This is useful when loading a
checkpoint saved during a run, where the models were not yet stored:
`result <- update(result)`.
