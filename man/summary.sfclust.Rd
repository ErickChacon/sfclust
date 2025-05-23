% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sfclust-methods.R
\name{summary.sfclust}
\alias{summary.sfclust}
\title{Summary method for sfclust objects}
\usage{
\method{summary}{sfclust}(object, sample = object$clust$id, sort = FALSE, ...)
}
\arguments{
\item{object}{An object of class 'sfclust'.}

\item{sample}{An integer specifying the clustering sample number to be summarized (default
is the last sample).}

\item{sort}{Logical value indicating if clusters should be relabel based on number of
elements.}

\item{...}{Additional arguments passed to \code{print.default}.}
}
\value{
Invisibly returns a table with the number of regions in each cluster for the selected
sample. The function also prints a summary that includes:
\itemize{
\item the within-cluster model formula,
\item the total number of MCMC clustering samples,
\item the cluster membership counts for the specified sample (optionally sorted),
\item and the log marginal likelihood of the selected clustering sample.
}
}
\description{
This function summarizes the cluster assignments from the desired clustering \code{sample}.
}
