## Submission

This is a minor release (1.0.1 -> 1.1.0) adding raster `stars` support, a new
built-in dataset (`chapa`), `update()` methods for continuing/refitting MCMC
samples, and several bug fixes. See NEWS.md for full details.

## Test environments

* local: Linux (Arch), R 4.6.0
* GitHub Actions: R CMD check (Linux, macOS, Windows), via .github/workflows
* R-hub

## R CMD check results

0 errors | 0 warnings | 0 notes

* Suggested package 'INLA' is not available on CRAN. It is required for the
  package's core Bayesian clustering procedure (`sfclust()`), which calls
  INLA at every MCMC iteration; without it, only the standalone
  preprocessing utilities (`genclust()`, `data_all()`) are usable. All
  INLA-dependent code in examples, tests, and vignettes is guarded by
  `requireNamespace("INLA", quietly = TRUE)`, so R CMD check passes cleanly
  without it. 'INLA' is installed from
  <https://inla.r-inla-download.org/R/stable>, as noted in the
  `Additional_repositories` field of DESCRIPTION.

## Downstream dependencies

There are currently no downstream dependencies for this package.
