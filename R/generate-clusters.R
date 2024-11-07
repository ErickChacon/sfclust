#' Generate clusters for spatial clustering
#'
#' Creates a undirected graph from spatial polygonal data, computes its minimum spanning
#' tree (MST) and generate `nclust` clusters.
#'
#' @param x An `sf` or `sfc` object representing the spatial polygonal data. Additionally,
#' it can also be a `matrix` or `Matrix` object with non-zero values representing the
#' weighted connectivity between units.
#' @param nclust Integer, the initial number of clusters.
#' @param weights Optional `numeric` vector or `matrix` of weights between units in `x`.
#' It should be of dimension `n^2`, where `n` is the number of units in `x`. If NULL,
#' random weights are assigned.
#'
#' @return A list with three elements: 'graph' containing the initial undirected graph
#' object, 'mst' containing the initial minimum spanning tree, and 'membership' containing
#' the membership of `x`.
#'
#' @examples
#'
#'   librar(sf)
#'   x <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))
#'   cluster_ini <- genclust(x, nclust = 3, weights = 1:6))
#'   print(cluster_ini)
#'
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
