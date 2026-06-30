library(sf)
library(stars)
library(igraph)
library(Matrix)

test_that("genclust input validation", {
  x <- sparseMatrix(i = 1:5, j = 2:6, x = 1, dims = c(6, 6), symmetric = TRUE)

  expect_error(genclust("x", nclust = 5), "No genclust method")
  expect_error(genclust(x, nclust = 0), "`nclust` must be a positive integer.")
  expect_error(genclust(x, nclust = -1), "`nclust` must be a positive integer.")
  expect_error(genclust(x, nclust = "3"), "`nclust` must be a positive integer.")
  expect_error(genclust(x, nclust = 10), "`nclust` must be smaller than number of regions.")
})

test_that("genclust for matrix input", {
  x <- sparseMatrix(i = 1:5, j = 2:6, x = 1, dims = c(6, 6), symmetric = TRUE)
  set.seed(42)
  cluster_ini <- genclust(x, nclust = 3, weights = 1:length(x))

  expect_equal(unname(as_adjacency_matrix(cluster_ini$graph)), as(x, "generalMatrix"))
  expect_equal(unname(as_adjacency_matrix(cluster_ini$mst)), as(x, "generalMatrix"))
  expect_equal(unname(cluster_ini$membership), c(1, 2, 2, 2, 2, 3))
})

test_that("genclust for sfc objects", {
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
})

test_that("genclust for stars raster input", {
  x <- st_as_stars(cluster = matrix(1:35, 5))

  set.seed(42)
  clust <- genclust(x, nclust = 4)

  expect_equal(length(clust$membership), 35)
  expect_false(is.null(clust$valid_ids))
  expect_equal(length(clust$valid_ids), 35)
  expect_equal(length(unique(clust$membership)), 4)
  expect_true(igraph::is_igraph(clust$graph))
  expect_true(igraph::is_igraph(clust$mst))
})

test_that("raster_adjacency builds correct 4-connected grid", {
  x <- raster_adjacency(3, 2)
  expected <- sparseMatrix(
    i = c(1, 2, 4, 5, 1, 2, 3),
    j = c(2, 3, 5, 6, 4, 5, 6),
    x = 1L, dims = c(6, 6), symmetric = TRUE
  )
  expect_equal(x, as(expected, "generalMatrix"))

  # 1x1 grid: no neighbors
  x <- raster_adjacency(1, 1)
  expected <- sparseMatrix(i = integer(0), j = integer(0), x = integer(0), dims = c(1, 1))
  expect_equal(x, expected)
})
