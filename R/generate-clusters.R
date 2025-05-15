#' Generate clusters for spatial clustering
#'
#' Creates an undirected graph from spatial polygonal data, computes its minimum spanning
#' tree (MST), and generates `nclust` clusters. This function is used to initialize
#' cluster membership in a clustering algorithm, such as `sfclust`.
#'
#' @param x An `sf` or `sfc` object representing spatial polygonal data. It can also be
#'   a `matrix` or `Matrix` object with non-zero values representing weighted
#'   connectivity between units.
#' @param nclust Integer, specifying the initial number of clusters.
#' @param weights Optional `numeric` vector or `matrix` of weights between units in `x`.
#'   It should have dimensions `n^2`, where `n` is the number of units in `x`. If NULL,
#'   random weights are assigned.
#'
#' @return A list with three elements:
#'   - `graph`: The undirected graph object representing spatial contiguity.
#'   - `mst`: The minimum spanning tree.
#'   - `membership`: The cluster membership for elements in `x`.
#'
#' @examples
#'
#' library(sfclust)
#' library(sf)
#'
#' x <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))
#'
#' # using distance between geometries
#' clust <- genclust(x, nclust = 3, weights = st_distance(st_centroid(x)))
#' print(clust)
#' plot(st_sf(x, cluster = factor(clust$membership)))
#'
#' # using increasing weights
#' cluster_ini <- genclust(x, nclust = 3, weights = 1:36)
#' print(cluster_ini)
#' plot(st_sf(x, cluster = factor(cluster_ini$membership)))
#'
#' # using on random weights
#' cluster_ini <- genclust(x, nclust = 3, weights = runif(36))
#' print(cluster_ini)
#' plot(st_sf(x, cluster = factor(cluster_ini$membership)))
#'
#' @import igraph
#' @importFrom methods as
#' @importFrom sf st_touches
#' @export
genclust <- function(x, nclust = 10, weights = NULL){

  # create adjacency, initial checks and weights if required
  if (inherits(x, c("sf", "sfc"))) {
    x <- as(st_touches(st_geometry(x)), "matrix")
  } else if (!inherits(x, c("matrix", "Matrix"))) {
    stop("`x` must be of class `sf`, `sfc`, `matrix` or `Matrix`.")
  }

  if (nclust > dim(x)[1]) {
    stop("`nclust` must be smaller that number of regions.")
  }

  if (is.null(weights)) weights <- runif(length(x))

  # create weighted graph and minimum spanning tree
  graph <- graph_from_adjacency_matrix(x * weights, mode = "upper", weighted = TRUE)
  mstgraph <- mst(graph)
  V(mstgraph)$vid <- 1:vcount(mstgraph)

  # partition mst into nclust groups
  rmid <- order(E(mstgraph)$weight, decreasing = TRUE)[1:(nclust - 1)]
  partition <- components(delete_edges(mstgraph, rmid))

  return(list(graph = graph, mst = mstgraph, membership = partition$membership))
}
