#' Generate initial cluster assignments
#'
#' Creates an undirected graph from spatial data or a weighted adjacency matrix,
#' computes its minimum spanning tree (MST), and partitions it into `nclust` clusters.
#' Accepts `sf` and `sfc` objects (contiguity via [sf::st_touches()]), weighted
#' `matrix` or `Matrix` objects, and `stars` rasters (4-connected grid adjacency for
#' non-NA pixels).
#'
#' @param x Spatial data or adjacency matrix. Accepted classes: `matrix`, `Matrix`,
#'   `sf`, `sfc` (vector geometry), or `stars` (raster).
#' @param nclust Integer. Number of initial clusters (default `10`).
#' @param weights Optional numeric vector or matrix of edge weights with `n^2` elements,
#'   where `n` is the number of spatial units. Smaller values indicate stronger similarity
#'   between units. If `NULL`, random weights are assigned. Not applicable for `stars` input.
#'
#' @return A list with:
#'   - `graph`: undirected igraph object representing spatial contiguity.
#'   - `mst`: minimum spanning tree of `graph`.
#'   - `membership`: integer vector of cluster assignments (length = number of spatial units).
#'   - `valid_ids`: integer vector of non-NA pixel indices for `stars` input; `NULL` otherwise.
#'
#' @examples
#'
#' library(sfclust)
#' library(sf)
#'
#' x <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))
#'
#' # sfc input with custom weights
#' clust <- genclust(x, nclust = 3, weights = as.numeric(st_distance(st_centroid(x))))
#' plot(st_sf(x, cluster = factor(clust$membership)))
#'
#' # matrix input (weighted adjacency)
#' adj <- as(sf::st_touches(x), "matrix") * runif(36)
#' clust <- genclust(adj, nclust = 3)
#' plot(st_sf(x, cluster = factor(clust$membership)))
#'
#' @export
genclust <- function(x, ...) UseMethod("genclust")

#' @rdname genclust
#' @export
genclust.default <- function(x, ...) {
  stop("No genclust method for class '", paste(class(x), collapse = "/"), "'. ",
       "Provide a matrix, Matrix, sf, sfc, or stars object.")
}

#' @rdname genclust
#' @export
genclust.matrix <- function(x, nclust = 10, weights = NULL) {
  validate_nclust(nclust, dim(x)[1])
  if (is.null(weights)) weights <- runif(length(x))
  result <- genclust_adj(x * weights, nclust = nclust)
  result$valid_ids <- NULL
  result
}

#' @rdname genclust
#' @export
genclust.Matrix <- function(x, nclust = 10, weights = NULL) {
  genclust.matrix(x, nclust = nclust, weights = weights)
}

#' @rdname genclust
#' @importFrom methods as
#' @importFrom sf st_touches
#' @export
genclust.sfc <- function(x, nclust = 10, weights = NULL) {
  adj <- as(st_touches(x), "matrix") * 1L
  genclust.matrix(adj, nclust = nclust, weights = weights)
}

#' @rdname genclust
#' @importFrom sf st_geometry
#' @export
genclust.sf <- function(x, nclust = 10, weights = NULL) {
  genclust.sfc(st_geometry(x), nclust = nclust, weights = weights)
}

#' @rdname genclust
#' @param nclust Integer. Number of initial clusters (default `10`).
#' @param sp_dims Character vector with the names of the two spatial dimensions. If
#'   `NULL`, auto-detected as dimensions with a non-`NA` regular `delta` in
#'   [stars::st_dimensions()].
#' @importFrom stars st_dimensions
#' @export
genclust.stars <- function(x, nclust = 10, sp_dims = NULL) {
  dims <- st_dimensions(x)
  if (is.null(sp_dims)) {
    sp_dims <- names(dims)[!is.na(sapply(dims, function(d) d$delta))]
    if (length(sp_dims) != 2)
      stop("Could not auto-detect 2 spatial dimensions from `stars` object. Provide `sp_dims`.")
  }
  sp_margins <- match(sp_dims, names(dim(x)))
  valid_ids  <- which(apply(x[[1]], sp_margins, function(v) !all(is.na(v))))
  A   <- raster_adjacency(dim(x)[sp_dims[1]], dim(x)[sp_dims[2]])
  adj <- A[valid_ids, valid_ids, drop = FALSE]
  result <- genclust_adj(adj * runif(length(adj)), nclust = nclust)
  result$valid_ids <- valid_ids
  result
}

validate_nclust <- function(nclust, n) {
  if (!is.numeric(nclust) || length(nclust) != 1 || nclust < 1)
    stop("`nclust` must be a positive integer.")
  if (nclust > n)
    stop("`nclust` must be smaller than number of regions.")
}

#' @importFrom igraph graph_from_adjacency_matrix mst V vcount ecount delete_edges components
#' @importFrom Matrix sparseMatrix
NULL

genclust_adj <- function(x, nclust = 10) {
  graph <- graph_from_adjacency_matrix(x, mode = "upper", weighted = TRUE)
  mstgraph <- mst(graph)
  V(mstgraph)$vid <- 1:vcount(mstgraph)
  rmid <- sample.int(ecount(mstgraph), nclust - 1)
  partition <- components(delete_edges(mstgraph, rmid))
  list(graph = graph, mst = mstgraph, membership = partition$membership)
}

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
