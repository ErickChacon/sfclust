# Fit models and compute the log marginal likelihood for all clusters

Fit the specified INLA model to each cluster and compute the log
marginal likelihood for each cluster specified in the membership vector.

## Usage

``` r
log_mlik_all(
  membership,
  data,
  correction = TRUE,
  detailed = FALSE,
  inla_args = NULL
)
```

## Arguments

- membership:

  Integer, character or factor vector indicating the cluster membership
  for each spatial unit.

- data:

  A long-format data frame as returned by [`filter_df()`](filter_df.md),
  with columns `id` (flat array position), `ids` (flat spatial
  position), `sid` (valid-cell rank, 1 to n_valid, matching `membership`
  indices), observation index columns (`id_<dimname>`), and all
  response/covariate variables.

- correction:

  Logical value indicating whether a correction for dispersion.

- detailed:

  Logical value indicating whether to return the INLA model instead of
  the log marginal likelihood. The argument `correction` is not applied
  in this case.

- inla_args:

  A named list or pairlist of arguments passed to `inla()` (e.g.
  `formula`, `family`, `E`).

## Value

A numeric vector containing the log marginal likelihood for each cluster
or the the fitted INLA model for each cluster when `detailed = TRUE`.
