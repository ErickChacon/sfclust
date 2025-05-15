#' @title sfclust: Bayesian Spatial Functional Clustering
#'
#' @description
#' The \code{sfclust} package implements the clustering algorithm proposed in *"Bayesian
#' Spatial Functional Data Clustering: Applications in Disease Surveillance"*, available
#' at <https://arxiv.org/abs/2407.12633>. The package provides tools for performing
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
