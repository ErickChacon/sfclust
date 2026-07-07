#' @title sfclust: Bayesian Spatial Functional Clustering
#'
#' @description
#' The \code{sfclust} package implements the clustering algorithm proposed in *"Bayesian
#' Spatial Functional Data Clustering: Applications in Disease Surveillance"*, available
#' at <https://doi.org/10.1002/sim.70597>. The package provides tools for performing
#' Bayesian spatial functional clustering, as well as methods for diagnostic analysis,
#' visualization, and summarization of results.
#'
#' @author
#' Ruiman Zhong \email{ruiman.zhong@kaust.edu.sa},
#' Erick A. Chacón-Montalván \email{erick.chaconmontalvan@kaust.edu.sa},
#' Paula Moraga \email{paula.moraga@kaust.edu.sa}
#'
#' @importFrom stats runif setNames
#' @importFrom sf st_touches st_geometry
#'
#' @docType package
#' @keywords internal
#' @aliases sfclust-package
"_PACKAGE"

#' Spatio-temporal Binomial data
#'
#' A simulated `stars` object containing binomial response data with a functional clustering
#' pattern defined by polynomial fixed effects. This dataset includes the variables `cases`
#' and `population` observed across 100 simulated spatial regions over 91 time points.
#'
#' @format A `stars` object with:
#' \describe{
#'   \item{cases}{Number of observed cases (integer)}
#'   \item{population}{Population at risk (integer)}
#'   \item{dimensions}{Two dimensions: \code{geometry} (spatial features) and \code{time} (daily observations)}
#' }
#'
#' @usage data(stbinom)
#'
#' @examples
#'
#' library(sfclust)
#'
#' data(stbinom)
#' stbinom
#' plot(stbinom["cases"])
#'
#' @name stbinom
"stbinom"

#' Spatio-temporal Gaussian data
#'
#' A simulated `stars` object containing Gaussian response data with a functional
#' clustering pattern uging random walk processes. This dataset includes the response
#' variable `y` observed across 100 simulated spatial regions over 91 time points.
#'
#' @format A `stars` object with:
#' \describe{
#'   \item{y}{Response variable}
#' }
#'
#' @usage data(stgaus)
#'
#' @examples
#'
#' library(sfclust)
#'
#' data(stgaus)
#' stgaus
#' plot(stgaus["y"])
#'
#' @name stgaus
"stgaus"

#' Spatio-temporal NDWI2 raster data from El Chaparrillo
#'
#' A `stars` raster object containing the Normalized Difference Water Index 2 (NDWI2)
#' derived from Sentinel-2 satellite imagery over El Chaparrillo, an agricultural
#' research center in Spain, across 86 time points.
#'
#' @format A `stars` object with:
#' \describe{
#'   \item{ndwi2}{NDWI2 values (numeric) derived from Sentinel-2 bands. Non-NA values define the study area mask.}
#'   \item{dimensions}{Three dimensions: \code{x} (26 columns), \code{y} (27 rows),
#'     \code{time} (86 dates from 2025-05-15 to 2025-12-01). CRS: WGS 84 / UTM zone 30N.}
#' }
#'
#' @usage data(chapa)
#'
#' @examples
#'
#' library(sfclust)
#'
#' data(chapa)
#' chapa
#' plot(chapa["ndwi2"])
#'
#' @name chapa
"chapa"

## ggplot variables used in NSE contexts
globalVariables(c("log_mlike", "time", "mean_cluster", "cluster"))
