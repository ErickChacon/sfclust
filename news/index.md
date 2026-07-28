# Changelog

## sfclust 1.1.0

### New features

- Added support for raster `stars` input (`x`/`y` grid dimensions),
  alongside the existing vector-geometry input, including a new
  `raster_adjacency()` internal helper and dispatch in
  [`genclust()`](../reference/genclust.md)/[`sfclust()`](../reference/sfclust.md).
- Added `chapa`, a real NDWI2 vegetation-index raster dataset from El
  Chaparrillo, Spain, and a new article vignette demonstrating raster
  clustering (`vg13-ndwi-chaparrillo.Rmd`).
- Cells with missing/`NA` response values (e.g. raster cells outside a
  study region) are now automatically excluded from clustering via
  `response` / `valid_ids`, instead of requiring the user to pre-filter
  their data.
- Added [`update()`](https://rdrr.io/r/stats/update.html) methods for
  `sfclust` objects: `niter` continues MCMC sampling from where a
  previous fit left off, and `sample` refits the full INLA models for a
  specific stored sample without further sampling.
- Generalized spatial/functional dimension handling (`spnames`/`fnames`)
  so spatial dimension names no longer need to match a fixed convention.
- [`plot_clusters_fitted()`](../reference/plot_clusters_fitted.md)/[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)
  gained `inv_link`/aggregation options, returning cluster-level
  (`mean_cluster`) and inverse-link (`mean_cluster_inv`) summaries.

### Changes

- Internal architecture reworked: [`sfclust()`](../reference/sfclust.md)
  now dispatches to [`sfclust.data.frame()`](../reference/sfclust.md)
  (core interface) and [`sfclust.stars()`](../reference/sfclust.md)
  (spatial convenience wrapper), replacing the previous single-path
  implementation.
  [`data_all()`](../reference/data_all.md)/[`filter_df()`](../reference/filter_df.md)
  replace the old `stnames`-based long-format conversion.

### Bug fixes

- Fixed corner cases with `NULL` fitted models and raster inputs.
- Fixed adjacency handling for custom/user-supplied adjacency matrices
  and point-geometry plotting.
- Fixed cluster-level aggregation in fitted values and plots.

## sfclust 1.0.1

CRAN release: 2025-05-19

- Initial CRAN release.
