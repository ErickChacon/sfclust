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
#' @param graph_mode Character. Use `"spatial"` to build the graph from spatial
#'   contiguity, or `"dense"` to allow every spatial unit to be connected. The
#'   default is `"spatial"`.
#' @param k Optional integer. When using `graph_mode = "dense"`, if `k` is provided, only the `k`
#'   nearest neighbors of each spatial unit are kept. This requires `weights`.
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
genclust.matrix <- function(x, nclust = 10, weights = NULL, graph_mode = c("spatial", "dense"), k = NULL, ...) {
  graph_mode <- match.arg(graph_mode)
  validate_nclust(nclust, nrow(x))

  if (is.null(weights)) weights <- matrix(runif(length(x)), nrow(x), ncol(x))
  if (!is.matrix(weights) && !inherits(weights, "Matrix")) {
    weights <- matrix(weights, nrow(x), ncol(x))
  }

  if (graph_mode == "dense") {
    adj <- weights
    if (!is.null(k)) {
      if (is.null(weights)) stop("'weights' must be provided to select 'k' nearest neighbors.")
      if (!is.numeric(k) || length(k) != 1 || is.na(k) || k != as.integer(k) || k < 1 || k > nrow(adj)) {
        stop("'k' must be a single integer between 1 and the number of regions.")
      }

      diag(adj) <- NA # skip self matches
      knn_mat <- matrix(0, nrow(adj), ncol(adj))
      for (i in seq_len(nrow(adj))) {
        knn <- order(adj[i, ], decreasing = FALSE)[seq_len(k)]
        knn_mat[i, knn] <- 1
        knn_mat[knn, i] <- 1
      }
      adj <- knn_mat * adj
    }
  } else {
    adj <- x * weights
  }
  diag(adj) <- 0
  genclust_adj(adj, nclust = nclust)
}

#' @rdname genclust
#' @export
genclust.Matrix <- function(x, nclust = 10, weights = NULL, graph_mode = c("spatial", "dense"), k = NULL, ...) {
  graph_mode <- match.arg(graph_mode)
  genclust.matrix(x, nclust = nclust, weights = weights, graph_mode = graph_mode, k = k)
}

#' @rdname genclust
#' @param nclust Integer. Number of initial clusters (default `10`).
#' @param spnames Character vector with the names of the spatial dimensions. If
#'   `NULL`, auto-detected as dimensions with a non-`NA` regular `delta` in
#'   [stars::st_dimensions()]. Currently supports 1D (vector geometry) and 2D
#'   raster grids; 3D spatial grids (e.g. `x`, `y`, `z`) are not yet supported.
#' @param response Character. Name of the attribute in `x` to use for determining
#'   valid spatial cells (cells where all observations are NA are excluded).
#'   If `NULL` (default), all spatial cells are treated as valid.
#' @param weights Optional numeric vector of edge weights with `n^2` elements,
#'   where `n` is the total number of spatial units (before filtering). If `NULL`,
#'   random weights are assigned.
#' @importFrom stars st_dimensions
#' @export
genclust.stars <- function(x, nclust = 10, spnames = NULL, response = NULL, weights = NULL, graph_mode = c("spatial", "dense"), k = NULL, ...) {
  graph_mode <- match.arg(graph_mode)
  spnames   <- detect_spnames(x, spnames)
  domain    <- create_domain(x, spnames, response)
  valid_ids <- which(domain[[1]])
  adj       <- create_adj(domain, weights, valid_ids, graph_mode, k)
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
create_adj <- function(domain, weights = NULL, valid_ids = which(domain[[1]]), graph_mode = c("spatial", "dense"), k = NULL) {
  graph_mode <- match.arg(graph_mode)
  spnames <- dimnames(domain)

  if (graph_mode == "dense") {
    n <- length(domain[[1]])
    adj <- matrix(1, n, n)
  } else {
    if (length(spnames) == 1) {
      geom <- st_get_dimension_values(domain, spnames)
      adj <- as(st_touches(geom), "sparseMatrix")
    } else if (length(spnames) == 2) {
      adj <- raster_adjacency(dim(domain)[[1]], dim(domain)[[2]])
    } else if (length(spnames) == 3) {
      stop("3D raster support (e.g. x, y, z brain voxels) is not yet implemented.")
    } else {
      stop("create_adj only supports 1 (geometry) or 2 (raster) spatial dimensions.")
    }
  }

  if (is.null(weights)) weights <- matrix(runif(length(adj)), nrow(adj), ncol(adj))
  if (!is.matrix(weights) && !inherits(weights, "Matrix")) {
    weights <- matrix(weights, nrow(adj), ncol(adj))
  }

  if (graph_mode == "dense") {
    adj <- weights
    if (!is.null(k)) {
      if (is.null(weights)) stop("'weights' must be provided to select 'k' nearest neighbors.")
      if (!is.numeric(k) || length(k) != 1 || is.na(k) || k != as.integer(k) || k < 1 || k > nrow(adj)) {
        stop("'k' must be a single integer between 1 and the number of regions.")
      }

      diag(adj) <- NA
      knn_mat <- matrix(0, nrow(adj), ncol(adj))
      for (i in seq_len(nrow(adj))) {
        knn <- order(adj[i, ], decreasing = FALSE)[seq_len(k)]
        knn_mat[i, knn] <- 1
        knn_mat[knn, i] <- 1
      }
      adj <- knn_mat * adj
    }
  } else {
    adj <- adj * weights
  }

  diag(adj) <- 0
  adj <- as(adj, "CsparseMatrix")
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

# Builds a 4-connected adjacency matrix for a 2D grid of size nx * ny.
# Flat position: (ix, iy) -> (iy - 1) * nx + ix  (first dim varies fastest).
# Extension to 3D (e.g. brain voxels) would follow the same pattern with
# 6-connected neighbors and flat position (iz-1)*nx*ny + (iy-1)*nx + ix.
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
