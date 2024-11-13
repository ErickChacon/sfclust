library(sf)
library(igraph)
library(Matrix)

test_that("generate clusters", {
  # sfc objects
  x <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))

  ## weights based in distance
  cluster_ini <- genclust(x, nclust = 3, weights = st_distance(st_centroid(x)))

  i <- c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5)
  j <- c(2, 4, 5, 3, 4, 5, 6, 5, 6, 5, 6)
  A <- sparseMatrix(i = i, j = j, x = 1, dims = c(6, 6), symmetric = TRUE)
  expect_equal(unname(as_adj(cluster_ini$graph)), as(A, "generalMatrix"))

  i <- c(1, 2, 3, 4, 5)
  j <- c(4, 3, 6, 5, 6)
  A <- sparseMatrix(i = i, j = j, x = 1, dims = c(6, 6), symmetric = TRUE)
  expect_equal(unname(as_adj(cluster_ini$mst)), as(A, "generalMatrix"))

  expect_equal(unname(cluster_ini$membership), c(1, 2, 3, 3, 3, 3))

  ## weights as sequence
  cluster_ini <- genclust(x, nclust = 3, weights = 1:length(x))

  i <- c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5)
  j <- c(2, 4, 5, 3, 4, 5, 6, 5, 6, 5, 6)
  A <- sparseMatrix(i = i, j = j, x = 1, dims = c(6, 6), symmetric = TRUE)
  expect_equal(unname(as_adj(cluster_ini$graph)), as(A, "generalMatrix"))

  i <- c(1, 1, 1, 2, 2)
  j <- c(2, 4, 5, 3, 6)
  A <- sparseMatrix(i = i, j = j, x = 1, dims = c(6, 6), symmetric = TRUE)
  expect_equal(unname(as_adj(cluster_ini$mst)), as(A, "generalMatrix"))

  expect_equal(unname(cluster_ini$membership), c(1, 1, 2, 1, 1, 3))

  # matrices
  x <- sparseMatrix(i = 1:5, j = 2:6, x = 1, dims = c(6, 6), symmetric = TRUE)

  ## weights as sequence
  cluster_ini <- genclust(x, nclust = 3, weights = 1:length(x))

  expect_equal(unname(as_adj(cluster_ini$graph)), as(x, "generalMatrix"))
  expect_equal(unname(as_adj(cluster_ini$mst)), as(x, "generalMatrix"))
  expect_equal(unname(cluster_ini$membership), c(1, 1, 1, 1, 2, 3))

  # missspecified x
  expect_error(genclust("x", nclust = 5),
    "`x` must be of class `sf`, `sfc`, `matrix` or `Matrix`.")

})

test_that("membership order", {
    mem <- c(4, 2, 4, 1, 3, 3, 1, 3)
    mem_new <- sort_membership(mem)
    expect_equal(as.integer(mem_new), c(3, 4, 3, 2, 1, 1, 2, 1))
    expect_equal((attr(mem_new, "order")), c(3, 1, 4, 2))

    mem <- c(4, 1, 1, 3, 2, 1, 2, 3)
    mem_new <- sort_membership(mem)
    expect_equal(as.integer(mem_new), mem)
    expect_equal((attr(mem_new, "order")), 1:4)
})
