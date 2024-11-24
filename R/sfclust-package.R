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

#' Example Dataset: toy
#'
#' A synthetic dataset used to demonstrate the functionality of the `sfclust` package.
#' This dataset contains spatially distributed functional data and associated attributes 
#' for testing and illustrating spatial functional clustering methods.
#'
#' @format A list with 7 elements:
#' \describe{
#'   \item{Y}{A numeric matrix of observed functional responses. Each row represents a spatial location, and each column corresponds to a time point.}
#'   \item{X1}{A numeric vector of the first covariate, representing spatially distributed explanatory data.}
#'   \item{X2}{A numeric vector of the second covariate, representing additional spatial information.}
#'   \item{map}{An `sf` object containing spatial geometries (e.g., polygons or points) that correspond to the locations of the data.}
#'   \item{E}{A numeric vector representing the expected number of cases at each spatial location. This is used for modeling relative risk.}
#'   \item{cluster}{An integer vector indicating the true cluster labels for each spatial location. Useful for validation and benchmarking.}
#'   \item{eta}{A numeric vector representing the intensity parameter for `Y`. This characterizes the underlying process generating the functional responses.}
#' }
#'
#' @usage data(toy)
#'
#' @examples
#' data(toy)
#' head(toy$Y)
#' plot(st_sf(toy$map), main = "Spatial Map")
#'
#' @keywords datasets
#' @name toy
NULL

#' Simulated Spatiotemporal and Spatial Data
#'
#' A list containing two objects for demonstrating Bayesian spatial functional clustering:
#'
#' @format A list with two elements:
#' \describe{
#'   \item{\code{stdata}}{A `stars` object containing spatiotemporal observations.}
#'   \item{\code{geodata}}{A `sf` object containing Voronoi polygons with cluster memberships.}
#' }
#'
#' @details
#' - \code{stdata} contains simulated spatiotemporal data where rows are time points and columns represent spatial units.
#' - \code{geodata} includes spatial Voronoi polygons and their associated metadata, such as cluster memberships.
#'
#' @examples
#' data(simu_data)
#'
#' # Access spatiotemporal data
#' stdata <- simu_data$stdata
#' plot(stdata)
#'
#' # Access spatial polygons
#' geodata <- simu_data$geodata
#' library(ggplot2)
#' ggplot(geodata) +
#'   geom_sf(aes(fill = factor(membership))) +
#'   labs(title = "Spatial Clusters", fill = "Cluster") +
#'   theme_minimal()
#'
#' @keywords datasets
#' @name simu_data
#' @docType data
NULL
