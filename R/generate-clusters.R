#' Generate initial cluster assignments
#'
#' Creates an undirected graph from spatial data or a weighted adjacency matrix,
#' computes its minimum spanning tree (MST), and partitions it into `nclust` clusters.
#' Accepts weighted `matrix` or `Matrix` objects, and `stars` objects with either
#' vector geometry (contiguity via [sf::st_touches()]) or raster dimensions
#' (4-connected grid adjacency for non-NA pixels).
#'
#' @param x Spatial data or adjacency matrix. Accepted classes: `matrix`, `Matrix`,
#'   or `stars`.
#' @param nclust Integer. Number of initial clusters (default `10`).
#' @param weights Optional numeric vector or matrix of edge weights with `n^2` elements,
#'   where `n` is the number of spatial units. Smaller values indicate stronger similarity
#'   between units. If `NULL`, random weights are assigned.
#' @param ... Not used; required for S3 method consistency.
#'
#' @return A list with:
#'   - `graph`: undirected igraph object representing spatial contiguity.
#'   - `mst`: minimum spanning tree of `graph`.
#'   - `membership`: integer vector of cluster assignments (length = number of valid spatial units).
#'   - `valid_ids`: integer vector of flat spatial positions included in the graph
#'     (only present for `stars` input; `NULL` for matrix/Matrix input).
#'
#' @examples
#'
#' library(sfclust)
#' library(sf)
#' library(stars)
#'
#' # stars object with vector geometry
#' geom <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))
#' x <- st_as_stars(st_sf(z = 1:6, geometry = geom))
#' clust <- genclust(x, nclust = 3)
#' plot(st_sf(geom, cluster = factor(clust$membership)))
#'
#' # stars raster input
#' x <- st_as_stars(cluster = matrix(1:35, 5))
#' clust <- genclust(x, nclust = 4)
#' x$cluster <- clust$membership
#' plot(x, col = rainbow(4))
#'
#' # matrix input
#' x <- matrix(c(0,1,0,1, 1,0,1,0, 0,1,0,1, 1,0,1,0), nrow = 4)
#' clust <- genclust(x, nclust = 2)
#' clust$membership
#'
#' @export
genclust <- function(x, ...) UseMethod("genclust")

#' @rdname genclust
#' @export
genclust.default <- function(x, ...) {
  stop("No genclust method for class '", paste(class(x), collapse = "/"), "'. ",
       "Provide a matrix, Matrix, or stars object.")
}

#' @rdname genclust
#' @export
genclust.matrix <- function(x, nclust = 10, weights = NULL, ...) {
  validate_nclust(nclust, nrow(x))
  if (is.null(weights)) weights <- runif(length(x))
  genclust_adj(x * weights, nclust = nclust)
}

#' @rdname genclust
#' @export
genclust.Matrix <- function(x, nclust = 10, weights = NULL, ...) {
  genclust.matrix(x, nclust = nclust, weights = weights)
}

#' @rdname genclust
#' @param nclust Integer. Number of initial clusters (default `10`).
#' @param spnames Character vector with the names of the spatial dimensions. If
#'   `NULL`, auto-detected as dimensions with a non-`NA` regular `delta` in
#'   [stars::st_dimensions()].
#' @param response Character. Name of the attribute in `x` to use for determining
#'   valid spatial cells (cells where all observations are NA are excluded).
#'   If `NULL` (default), all spatial cells are treated as valid.
#' @param weights Optional numeric vector of edge weights with `n^2` elements,
#'   where `n` is the total number of spatial units (before filtering). If `NULL`,
#'   random weights are assigned.
#' @importFrom stars st_dimensions
#' @export
genclust.stars <- function(x, nclust = 10, spnames = NULL, response = NULL, weights = NULL, ...) {
  spnames   <- detect_spnames(x, spnames)
  domain    <- create_domain(x, spnames, response)
  valid_ids <- which(domain[[1]])
  adj       <- create_adj(domain, weights, valid_ids)
  validate_nclust(nclust, nrow(adj))
  result <- genclust_adj(adj, nclust = min(nclust, nrow(adj)))
  result$valid_ids <- valid_ids
  result
}

validate_nclust <- function(nclust, n) {
  if (!is.numeric(nclust) || length(nclust) != 1 || nclust < 1)
    stop("`nclust` must be a positive integer.")
  if (nclust > n)
    stop("`nclust` must be smaller than number of valid spatial units.")
}

#' @importFrom stars st_apply
create_domain <- function(x, spnames, response = NULL) {
  if (is.null(response)) {
    st_apply(x[1], spnames, function(v) TRUE)
  } else {
    st_apply(x[response], spnames, function(v) any(!is.na(v)))
  }
}

#' @importFrom stars st_dimensions st_get_dimension_values
#' @importFrom methods as
#' @importFrom sf st_touches
create_adj <- function(domain, weights = NULL, valid_ids = which(domain[[1]])) {
  spnames <- dimnames(domain)
  if (length(spnames) == 1) {
    geom <- st_get_dimension_values(domain, spnames)
    adj <- as(st_touches(geom), "sparseMatrix")
  } else if (length(spnames) == 2) {
    adj <- raster_adjacency(dim(domain)[[1]], dim(domain)[[2]])
  } else {
    stop("create_adj only supports 1 (geometry) or 2 (raster x/y) spatial dimensions.")
  }
  if (is.null(weights)) weights <- runif(length(adj))
  adj <- as(adj * weights, "CsparseMatrix")
  adj[valid_ids, valid_ids, drop = FALSE]
}

#' @importFrom igraph graph_from_adjacency_matrix mst V "V<-" vcount ecount delete_edges components
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
