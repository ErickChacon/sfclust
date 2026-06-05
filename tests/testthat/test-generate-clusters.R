library(sf)
library(igraph)
library(Matrix)

test_that("generate clusters", {
  # sfc objects
  x <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))

  ## weights based in distance
  set.seed(42)
  cluster_ini <- genclust(x, nclust = 3, weights = st_distance(st_centroid(x)))

  i <- c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5)
  j <- c(2, 4, 5, 3, 4, 5, 6, 5, 6, 5, 6)
  A <- sparseMatrix(i = i, j = j, x = 1, dims = c(6, 6), symmetric = TRUE)
  expect_equal(unname(as_adjacency_matrix(cluster_ini$graph)), as(A, "generalMatrix"))

  i <- c(1, 2, 3, 4, 5)
  j <- c(4, 3, 6, 5, 6)
  A <- sparseMatrix(i = i, j = j, x = 1, dims = c(6, 6), symmetric = TRUE)
  expect_equal(unname(as_adjacency_matrix(cluster_ini$mst)), as(A, "generalMatrix"))

  expect_equal(unname(cluster_ini$membership), c(1, 2, 2, 3, 3, 2))

  ## weights as sequence
  set.seed(42)
  cluster_ini <- genclust(x, nclust = 3, weights = 1:length(x)^2)

  i <- c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5)
  j <- c(2, 4, 5, 3, 4, 5, 6, 5, 6, 5, 6)
  A <- sparseMatrix(i = i, j = j, x = 1, dims = c(6, 6), symmetric = TRUE)
  expect_equal(unname(as_adjacency_matrix(cluster_ini$graph)), as(A, "generalMatrix"))

  i <- c(1, 1, 1, 2, 2)
  j <- c(2, 4, 5, 3, 6)
  A <- sparseMatrix(i = i, j = j, x = 1, dims = c(6, 6), symmetric = TRUE)
  expect_equal(unname(as_adjacency_matrix(cluster_ini$mst)), as(A, "generalMatrix"))

  expect_equal(unname(cluster_ini$membership), c(1, 2, 2, 1, 1, 3))

  # matrices
  x <- sparseMatrix(i = 1:5, j = 2:6, x = 1, dims = c(6, 6), symmetric = TRUE)

  ## weights as sequence
  set.seed(42)
  cluster_ini <- genclust(x, nclust = 3, weights = 1:length(x))

  expect_equal(unname(as_adjacency_matrix(cluster_ini$graph)), as(x, "generalMatrix"))
  expect_equal(unname(as_adjacency_matrix(cluster_ini$mst)), as(x, "generalMatrix"))
  expect_equal(unname(cluster_ini$membership), c(1, 2, 2, 2, 2, 3))

  # missspecified x
  expect_error(genclust("x", nclust = 5),
    "`x` must be of class `stars`, `sf`, `sfc`, `matrix` or `Matrix`.")

  # invalid nclust values
  expect_error(genclust(x, nclust = 0), "`nclust` must be a positive integer.")
  expect_error(genclust(x, nclust = -1), "`nclust` must be a positive integer.")
  expect_error(genclust(x, nclust = "3"), "`nclust` must be a positive integer.")
})

test_that("raster_adjacency builds correct 4-connected grid", {
  # 3x2 grid: 6 cells
  A <- raster_adjacency(3, 2)
  expect_equal(dim(A), c(6, 6))
  expect_true(Matrix::isSymmetric(A))

  # cell ids: (ix, iy) -> id = (iy-1)*3 + ix
  # (1,1)=1, (2,1)=2, (3,1)=3, (1,2)=4, (2,2)=5, (3,2)=6
  # horizontal: 1-2, 2-3, 4-5, 5-6
  # vertical:   1-4, 2-5, 3-6
  expected_edges <- list(c(1,2), c(2,3), c(4,5), c(5,6), c(1,4), c(2,5), c(3,6))
  for (e in expected_edges) {
    expect_equal(A[e[1], e[2]], 1)
    expect_equal(A[e[2], e[1]], 1)
  }
  # no diagonal connections
  expect_equal(A[1, 5], 0)
  expect_equal(A[2, 4], 0)

  # 1x1 grid: no neighbors
  A1 <- raster_adjacency(1, 1)
  expect_equal(dim(A1), c(1, 1))
  expect_equal(A1[1, 1], 0)
})
