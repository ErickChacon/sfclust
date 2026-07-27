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

* Suggested package 'INLA' is not available on CRAN; it is used only in
  examples and vignettes guarded by `requireNamespace("INLA", quietly = TRUE)`,
  and the package is fully functional without it for non-INLA workflows.
  'INLA' is installed from <https://inla.r-inla-download.org/R/stable>, as
  noted in the `Additional_repositories` field of DESCRIPTION.

## Downstream dependencies

There are currently no downstream dependencies for this package.
