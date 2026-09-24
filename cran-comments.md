## Submission

This is a patch release (1.1.0 -> 1.1.1) fixing
two vignettes (documentation only, no code or
API changes). See NEWS.md for full details.

## Test environments

* local: Ubuntu 24.04 (Docker), R 4.6.0, via
  `R CMD check --as-cran`
* GitHub Actions: R CMD check (Linux, macOS,
  Windows; R release and oldrel-1, plus
  R-devel on Linux)
* R-hub

## R CMD check results

0 errors | 0 warnings | 0 notes

## Note for reviewers

Suggested package 'INLA' is not on CRAN; it is
available via `Additional_repositories`
(<https://inla.r-inla-download.org/R/stable>)
and used conditionally, as in previous
accepted versions.

## Downstream dependencies

There are currently no downstream dependencies
for this package.
