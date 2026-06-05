#' Generate clusters for spatial clustering
#'
#' Creates an undirected graph from spatial polygonal data, computes its minimum spanning
#' tree (MST), and generates `nclust` clusters. This function is used to initialize
#' cluster membership in a clustering algorithm, such as `sfclust`.
#'
#' @param x An `sf`, `sfc`, or `stars` object representing spatial data. It can also be
#'   a `matrix` or `Matrix` object with non-zero values representing weighted
#'   connectivity between units.
#' @param nclust Integer, specifying the initial number of clusters.
#' @param weights Optional `numeric` vector or `matrix` of weights between units in `x`.
#'   It should have dimensions `n^2`, where `n` is the number of units in `x`. If NULL,
#'   random weights are assigned.
#' @param sp_dims Character vector with the names of the spatial dimensions when `x` is a
#'   `stars` object. If `NULL`, spatial dimensions are auto-detected as those with a
#'   non-`NA` regular `delta` in `st_dimensions(x)`.
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
#' @importFrom stars st_dimensions
#' @export
genclust <- function(x, nclust = 10, weights = NULL, sp_dims = NULL){

  # create adjacency, initial checks and weights if required
  valid_ids <- NULL
  if (inherits(x, "stars")) {
    dims <- st_dimensions(x)
    if (is.null(sp_dims)) {
      sp_dims <- names(dims)[!is.na(sapply(dims, function(d) d$delta))]
      if (length(sp_dims) != 2)
        stop("Could not auto-detect 2 spatial dimensions from `stars` object. Provide `sp_dims`.")
    }
    sp_margins <- match(sp_dims, names(dim(x)))
    valid_ids <- which(apply(x[[1]], sp_margins, function(v) !all(is.na(v))))
    A <- raster_adjacency(dim(x)[sp_dims[1]], dim(x)[sp_dims[2]])
    x <- A[valid_ids, valid_ids, drop = FALSE]
  } else if (inherits(x, c("sf", "sfc"))) {
    x <- as(st_touches(st_geometry(x)), "matrix")
  } else if (!inherits(x, c("matrix", "Matrix"))) {
    stop("`x` must be of class `stars`, `sf`, `sfc`, `matrix` or `Matrix`.")
  }

  if (!is.numeric(nclust) || length(nclust) != 1 || nclust < 1) {
    stop("`nclust` must be a positive integer.")
  }
  if (nclust > dim(x)[1]) {
    stop("`nclust` must be smaller than number of regions.")
  }

  if (is.null(weights)) weights <- runif(length(x))
  result <- genclust_adj(x * weights, nclust = nclust)
  result$valid_ids <- valid_ids
  return(result)
}

# Core MST-based cluster initialization from a pre-weighted adjacency matrix.
# Returns list(graph, mst, membership).
genclust_adj <- function(x, nclust = 10) {
  graph <- graph_from_adjacency_matrix(x, mode = "upper", weighted = TRUE)
  mstgraph <- mst(graph)
  V(mstgraph)$vid <- 1:vcount(mstgraph)
  rmid <- sample.int(ecount(mstgraph), nclust - 1)
  partition <- components(delete_edges(mstgraph, rmid))
  list(graph = graph, mst = mstgraph, membership = partition$membership)
}

#' Build a 4-connected (Rook) adjacency matrix for a regular grid
#'
#' Constructs the sparse symmetric adjacency matrix for an `nx` by `ny` regular grid
#' under 4-connectivity (each cell connected to its left, right, top, and bottom
#' neighbors). The flat cell index is `(iy - 1) * nx + ix`.
#'
#' @param nx Integer. Number of cells in the x direction (columns).
#' @param ny Integer. Number of cells in the y direction (rows).
#'
#' @return A sparse symmetric `Matrix` of size `(nx * ny) x (nx * ny)` with 1s at
#'         adjacent cell pairs.
#'
#' @examples
#' library(sfclust)
#' A <- raster_adjacency(3, 2)
#' print(A)
#'
#' @importFrom Matrix sparseMatrix
#' @export
raster_adjacency <- function(nx, ny) {
  # horizontal neighbors: (ix, iy) -- (ix+1, iy)
  ix <- rep(seq_len(nx - 1L), ny)
  iy <- rep(seq_len(ny), each = nx - 1L)
  from_h <- (iy - 1L) * nx + ix
  to_h   <- from_h + 1L

  # vertical neighbors: (ix, iy) -- (ix, iy+1)
  ix <- rep(seq_len(nx), ny - 1L)
  iy <- rep(seq_len(ny - 1L), each = nx)
  from_v <- (iy - 1L) * nx + ix
  to_v   <- from_v + nx

  from <- c(from_h, from_v)
  to   <- c(to_h,   to_v)
  n    <- nx * ny

  Matrix::sparseMatrix(i = c(from, to), j = c(to, from), x = 1L, dims = c(n, n))
}
