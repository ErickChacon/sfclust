% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model-within.R
\name{data_each}
\alias{data_each}
\title{Prepare data for a cluster}
\usage{
data_each(k, membership, stdata, stnames = c("geometry", "time"))
}
\arguments{
\item{k}{The cluster number to subset.}

\item{membership}{A vector defining the cluster membership for each region.}

\item{stdata}{A stars object containing spatial-temporal dimensions defined in \code{stnames}.}

\item{stnames}{The names of the \code{spatial} and \code{temporal} dimensions.}
}
\value{
A long-format data frame with ids for each observation and  for spatial and
time indexing.
}
\description{
Subset a spatio-temporal dataset for a cluster and convert it to a long format with
indices for time and spatial location.
}
\examples{

library(sfclust)
library(stars)

dims <- st_dimensions(
  geometry = st_sfc(lapply(1:5, function(i) st_point(c(i, i)))),
  time = seq(as.Date("2024-01-01"), by = "1 day", length.out = 3)
)
stdata <- st_as_stars(cases = array(1:15, dim = c(5, 3)), dimensions = dims)

data_each(k = 2, membership = c(1, 1, 1, 2, 2), stdata)

}
