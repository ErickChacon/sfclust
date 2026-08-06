library(sf)
library(stars)
library(igraph)
library(Matrix)

# --- raster_adjacency -------------------------------------------------

test_that("raster_adjacency: 4-connected grid", {
  # 1x1: no neighbors
  x <- raster_adjacency(1, 1)
  A <- sparseMatrix(i = integer(0), j = integer(0), x = integer(0), dims = c(1, 1))
  expect_equal(x, A)

  # 4x1: chain (max 2 neighbors)
  x <- raster_adjacency(4, 1)
  A <- sparseMatrix(
    i = c(1, 2, 3), j = c(2, 3, 4),
    x = 1L, dims = c(4, 4), symmetric = TRUE
  )
  expect_equal(x, as(A, "generalMatrix"))

  # 3x2: corner and edge cells (max 3 neighbors)
  x <- raster_adjacency(3, 2)
  A <- sparseMatrix(
    i = c(1, 2, 4, 5, 1, 2, 3),
    j = c(2, 3, 5, 6, 4, 5, 6),
    x = 1L, dims = c(6, 6), symmetric = TRUE
  )
  expect_equal(x, as(A, "generalMatrix"))

  # 2x3: swapped dimensions
  x <- raster_adjacency(2, 3)
  A <- sparseMatrix(
    i = c(1, 1, 2, 3, 3, 4, 5),
    j = c(2, 3, 4, 4, 5, 6, 6),
    x = 1L, dims = c(6, 6), symmetric = TRUE
  )
  expect_equal(x, as(A, "generalMatrix"))

  # 3x4: interior cells with all 4 neighbors
  x <- raster_adjacency(3, 4)
  A <- sparseMatrix(
    i = c(1, 1, 2, 2, 3, 4, 4, 5, 5, 6,  7,  7,  8,  8,  9, 10, 11),
    j = c(2, 4, 3, 5, 6, 5, 7, 6, 8, 9,  8, 10,  9, 11, 12, 11, 12),
    x = 1L, dims = c(12, 12), symmetric = TRUE
  )
  expect_equal(x, as(A, "generalMatrix"))
})

# --- genclust: matrix -------------------------------------------------

test_that("genclust: input validation", {
  x <- sparseMatrix(i = 1:5, j = 2:6, x = 1, dims = c(6, 6), symmetric = TRUE)

  expect_error(genclust("x", nclust = 5), "No genclust method")
  expect_error(genclust(x, nclust = 0), "`nclust` must be a positive integer.")
  expect_error(genclust(x, nclust = -1), "`nclust` must be a positive integer.")
  expect_error(genclust(x, nclust = "3"), "`nclust` must be a positive integer.")
  expect_error(genclust(x, nclust = 10), "`nclust` must be smaller than number of valid spatial units.")
})

test_that("genclust: matrix input", {
  x <- sparseMatrix(i = 1:5, j = 2:6, x = 1, dims = c(6, 6), symmetric = TRUE)
  set.seed(42)
  clust <- genclust(x, nclust = 3, weights = 1:length(x))

  expect_equal(as_adjacency_matrix(clust$graph), as(x, "generalMatrix"))
  expect_equal(as_adjacency_matrix(clust$mst), as(x, "generalMatrix"))
  expect_equal(clust$membership, c(1, 2, 2, 2, 2, 3))
})

test_that("genclust: dense matrix input", {
  x <- sparseMatrix(i = 1:3, j = 2:4, x = 1, dims = c(4, 4), symmetric = TRUE)
  weights <- as.matrix(dist(seq_len(4), diag = TRUE, upper = TRUE))

  set.seed(1)
  clust <- genclust(x, nclust = 2, weights = weights, graph_mode = "dense", k = 1)

  A <- sparseMatrix(
    i = 1:3, j = 2:4,
    x = 1, dims = c(4, 4), symmetric = TRUE
  )
  graph_adj <- as_adjacency_matrix(clust$graph)
  dimnames(graph_adj) <- NULL
  expect_equal(graph_adj, as(A, "generalMatrix"))
  expect_equal(length(clust$membership), 4L)

  expect_error(
    genclust(x, nclust = 2, graph_mode = "dense", k = 1),
    "Weights must be provided to select k nearest neighbors.",
    fixed = TRUE
  )
})

test_that("genclust: dense stars input", {
  geom <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(4, 1))
  x <- st_as_stars(st_sf(val = 1:4, geometry = geom))

  set.seed(1)
  clust <- genclust(x, nclust = 2, graph_mode = "dense")

  expect_equal(ecount(clust$graph), choose(length(geom), 2))
  expect_equal(length(clust$membership), length(geom))
  expect_equal(clust$valid_ids, seq_along(geom))
})

# --- genclust: stars --------------------------------------------------

test_that("genclust: stars with vector geometry", {
  geom <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))
  x <- st_as_stars(st_sf(val = 1:6, geometry = geom))

  # weights based on distance
  set.seed(42)
  clust <- genclust(x, nclust = 3, weights = st_distance(st_centroid(geom)))

  A <- sparseMatrix(
    i = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5),
    j = c(2, 4, 5, 3, 4, 5, 6, 5, 6, 5, 6),
    x = 1, dims = c(6, 6), symmetric = TRUE
  )
  expect_equal(as_adjacency_matrix(clust$graph), as(A, "generalMatrix"))
  A <- sparseMatrix(
    i = c(1, 2, 3, 4, 5),
    j = c(4, 3, 6, 5, 6),
    x = 1, dims = c(6, 6), symmetric = TRUE
  )
  expect_equal(as_adjacency_matrix(clust$mst), as(A, "generalMatrix"))
  expect_equal(clust$membership, c(1, 2, 3, 1, 1, 3))
  expect_equal(clust$valid_ids, 1:6)

  # weights as sequence
  set.seed(42)
  clust <- genclust(x, nclust = 3, weights = 1:6^2)

  A <- sparseMatrix(
    i = c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5),
    j = c(2, 4, 5, 3, 4, 5, 6, 5, 6, 5, 6),
    x = 1, dims = c(6, 6), symmetric = TRUE
  )
  expect_equal(as_adjacency_matrix(clust$graph), as(A, "generalMatrix"))
  A <- sparseMatrix(
    i = c(1, 1, 1, 2, 2),
    j = c(2, 4, 5, 3, 6),
    x = 1, dims = c(6, 6), symmetric = TRUE
  )
  expect_equal(as_adjacency_matrix(clust$mst), as(A, "generalMatrix"))
  expect_equal(clust$membership, c(1, 2, 2, 1, 1, 3))
  expect_equal(clust$valid_ids, 1:6)
})

test_that("genclust: stars with vector geometry and NA region", {
  geom <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(4, 3))
  vals <- 1:12; vals[c(3, 7, 10)] <- NA
  x <- st_as_stars(st_sf(val = vals, geometry = geom))

  set.seed(21)
  clust <- genclust(x, nclust = 3, response = "val")

  expect_equal(clust$valid_ids, c(1, 2, 4, 5, 6, 8, 9, 11, 12))
  expect_equal(clust$membership, c(1, 1, 2, 1, 3, 2, 1, 3, 3))
  A <- sparseMatrix(
    i = c(1, 1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6, 8),
    j = c(2, 4, 5, 4, 5, 6, 5, 7, 7, 8, 8, 9, 9),
    x = 1, dims = c(9, 9), symmetric = TRUE
  )
  expect_equal(as_adjacency_matrix(clust$graph), as(A, "generalMatrix"))
  A <- sparseMatrix(
    i = c(1, 1, 1, 3, 4, 5, 6, 8),
    j = c(2, 4, 5, 6, 7, 8, 9, 9),
    x = 1, dims = c(9, 9), symmetric = TRUE
  )
  expect_equal(as_adjacency_matrix(clust$mst), as(A, "generalMatrix"))
})

test_that("genclust: stars with raster", {
  x <- st_as_stars(matrix(1:12, 4, 3))

  set.seed(7)
  clust <- genclust(x, nclust = 3)

  A <- sparseMatrix(
    i = c(1, 1, 2, 2, 3, 3, 4, 5, 5, 6, 6,  7,  7,  8,  9, 10, 11),
    j = c(2, 5, 3, 6, 4, 7, 8, 6, 9, 7, 10, 8, 11, 12, 10, 11, 12),
    x = 1, dims = c(12, 12), symmetric = TRUE
  )
  expect_equal(as_adjacency_matrix(clust$graph), as(A, "generalMatrix"))
  A <- sparseMatrix(
    i = c(1, 2, 3, 4, 5, 5, 6,  7,  9, 10, 11),
    j = c(2, 3, 7, 8, 6, 9, 7,  8, 10, 11, 12),
    x = 1, dims = c(12, 12), symmetric = TRUE
  )
  expect_equal(as_adjacency_matrix(clust$mst), as(A, "generalMatrix"))
  expect_equal(clust$membership, c(1, 1, 2, 2, 3, 2, 2, 2, 3, 3, 3, 3))
  expect_equal(clust$valid_ids, 1:12)
})

test_that("genclust: stars with raster and NA cells", {
  x <- st_as_stars(val = matrix(1:20, 5, 4))
  x$val[1, c(1, 2, 4)] <- NA
  x$val[2, 1] <- NA
  x$val[4, c(1,2)]       <- NA
  x$val[5, c(1, 2)] <- NA

  set.seed(42)
  clust <- genclust(x, nclust = 3, response = "val")

  expect_equal(clust$valid_ids, c(3, 7, 8, 11, 12, 13, 14, 15, 17, 18, 19, 20))
  expect_equal(clust$membership, c(1, 1, 1, 2, 2, 2, 3, 3, 2, 3, 3, 3))
  A <- sparseMatrix(
    i = c(1, 2, 2, 3, 4, 5, 5,  6,  6,  7,  7,  8,  9, 10, 11),
    j = c(3, 3, 5, 6, 5, 6, 9,  7, 10,  8, 11, 12, 10, 11, 12),
    x = 1, dims = c(12, 12), symmetric = TRUE
  )
  expect_equal(as_adjacency_matrix(clust$graph), as(A, "generalMatrix"))
  A <- sparseMatrix(
    i = c(1, 2, 3, 4, 5, 5,  6,  7,  7, 10, 11),
    j = c(3, 3, 6, 5, 6, 9, 10,  8, 11, 11, 12),
    x = 1, dims = c(12, 12), symmetric = TRUE
  )
  expect_equal(as_adjacency_matrix(clust$mst), as(A, "generalMatrix"))
})
