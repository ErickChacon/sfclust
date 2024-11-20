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
#' @docType package
#' @keywords internal
#' @aliases sfclust-package
"_PACKAGE"

#' US COVID-19 Data
#'
#' This dataset contains COVID-19 case counts and population data for different regions in the United States.
#'
#' @format A list with the following components:
#' \describe{
#'   \item{Y}{Matrix of COVID-19 case counts.}
#'   \item{N}{Matrix of population data.}
#'   \item{map}{Spatial map of the regions.}
#'   \item{E}{Matrix of expected cases data.}
#' }
#' @source Processed data from the US COVID-19 data
#' @name us_covid
#' @docType data
"us_covid"
