# sfclust 1.1.1

## Documentation fixes

* Rewrote the "Advanced features" vignette (`vg02-advanced-features.Rmd`) to
  use `nclust` directly instead of an explicit `genclust()` call, and to
  demonstrate the `path_save`/`nsave` checkpointing arguments; regenerated
  the corresponding cached vignette data (`gaussian-mcmc1.rds`,
  `gaussian-mcmc2.rds`).
* Fixed the "Standardized daily temperature at Canadian stations" vignette
  (`vg12-temp-canada.Rmd`): added `sort = TRUE` to `plot_clusters_series()` so
  its cluster numbering matches `plot_clusters_map()`/`plot_clusters_fitted()`,
  and corrected the resulting narrative description of cluster shapes.

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
