All notable changes to this project will be documented in this file.

# sfclust (development version)

## Added

* Added `sfclust_within_model()` for INLA within-cluster model specifications.
* Added a `within_model` argument to `sfclust()` for custom INLA models.
* Added support for joint INLA models using INLA stacks  while preserving the original univariate `formula` and `family` interface.
* Extended fitted value and plotting methods to support joint models when a compatible `fitted_fun` is provided.
* Added dense graph construction in `genclust()` for `matrix`, `Matrix`, and `stars` input through `graph_mode = "dense"`.

## Changed

* Kept the original `sfclust()` interface for univariate INLA models while allowing advanced users to define custom within-cluster model.

## Fixed

* Handled repeated INLA model fitting failures during proposal evaluation with a clear stopping rule controlled by `max_inla_failures`.

# sfclust 1.1.0

## New features

* Added support for raster `stars` input (`x`/`y` grid dimensions), alongside
  the existing vector-geometry input, including a new `raster_adjacency()`
  internal helper and dispatch in `genclust()`/`sfclust()`.
* Added `chapa`, a real NDWI2 vegetation-index raster dataset from El
  Chaparrillo, Spain, and a new article vignette demonstrating raster
  clustering (`vg13-ndwi-chaparrillo.Rmd`).
* Cells with missing/`NA` response values (e.g. raster cells outside a study
  region) are now automatically excluded from clustering via `response` /
  `valid_ids`, instead of requiring the user to pre-filter their data.
* Added `update()` methods for `sfclust` objects: `niter` continues MCMC
  sampling from where a previous fit left off, and `sample` refits the full
  INLA models for a specific stored sample without further sampling.
* Generalized spatial/functional dimension handling (`spnames`/`fnames`) so
  spatial dimension names no longer need to match a fixed convention.
* `plot_clusters_fitted()`/`fitted()` gained `inv_link`/aggregation options,
  returning cluster-level (`mean_cluster`) and inverse-link
  (`mean_cluster_inv`) summaries.

## Changes

* Internal architecture reworked: `sfclust()` now dispatches to
  `sfclust.data.frame()` (core interface) and `sfclust.stars()` (spatial
  convenience wrapper), replacing the previous single-path implementation.
  `data_all()`/`filter_df()` replace the old `stnames`-based long-format
  conversion.

## Bug fixes

* Fixed corner cases with `NULL` fitted models and raster inputs.
* Fixed adjacency handling for custom/user-supplied adjacency matrices and
  point-geometry plotting.
* Fixed cluster-level aggregation in fitted values and plots.
* `fitted()`/`plot()` no longer require INLA to be installed when applied to
  an already-fitted `sfclust` object (e.g. one loaded from a saved `.rds`);
  the inverse-link transform for standard families (gaussian, binomial,
  poisson) is now computed natively instead of via `INLA::inla.link.inv*()`.

# sfclust 1.0.1

* Initial CRAN release.
