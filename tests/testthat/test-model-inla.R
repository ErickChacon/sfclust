library(stars)
library(Matrix)

# --- unique_clusters ---------------------------------------------------

test_that("unique_clusters: integer membership", {
  membership <- c(2, 4, 4, 3, 1)
  expect_equal(unique_clusters(membership), 1:4)

  membership <- c(3, 7, 7, 4, 2)
  expect_error(unique_clusters(membership))
})

test_that("unique_clusters: factor and character membership", {
  membership <- factor(c(3, 7, 7, 4, 2))
  aux <- as.character(c(2, 3, 4, 7))
  expect_equal(unique_clusters(membership), setNames(aux, aux))

  membership <- factor(c("c", "e", "e", "d", "b"))
  aux <- c("b", "c", "d", "e")
  expect_equal(unique_clusters(membership), setNames(aux, aux))

  membership <- c("c", "e", "e", "d", "b")
  aux <- c("b", "c", "d", "e")
  expect_equal(unique_clusters(membership), setNames(aux, aux))
})

# --- correction_required -----------------------------------------------

test_that("correction_required: identifies rw1/rw2/crw1 terms", {
  expect_equal(correction_required(y ~ x + z), character())
  expect_equal(correction_required(y ~ x + f(z, model = "rw1")), "z")
  expect_equal(correction_required(y ~ f(z, model = "rw2") + x), "z")

  formula <- y ~ x +
    f(z, model = "crw1", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.01)))) +
    f(w, model = "rw2") +
    f(v, model = "ar")
  expect_equal(correction_required(formula), c("z", "w"))
})

# --- log_mlik_correction and get_structure_matrix ----------------------

test_that("log_mlik_correction: rw1, rw2, and combined", {
  skip_on_cran()

  n <- 10
  data <- data.frame(y = rnorm(n), time = 1:n, time2 = 1:n)

  ## rw1
  formula <- y ~ f(time, model = "rw1")
  model <- INLA::inla(formula, data = data, control.compute = list(config = TRUE))

  i <- c(1:n, 1:(n-1))
  j <- c(1:n, 2:n)
  vals <- c(c(1, rep(2, n-2), 1) + 0.0001, rep(-1, n-1))
  expect_equal(
    get_structure_matrix(model)$time,
    sparseMatrix(i = i, j = j, x = vals)
  )
  expect_equal(log_mlik_correction(model), -3.45305292)

  ## rw2
  formula <- y ~ f(time, model = "rw2")
  model <- INLA::inla(formula, data = data, control.compute = list(config = TRUE))

  i <- c(1:n, 1:(n-1), 1:(n-2))
  j <- c(1:n, 2:n, 3:n)
  vals <- c(c(1, 5, rep(6, n-4), 5, 1) + 0.0001, c(-2, rep(-4, n-3), -2), rep(1, n-2))
  expect_equal(
    get_structure_matrix(model)$time,
    sparseMatrix(i = i, j = j, x = vals)
  )
  expect_equal(log_mlik_correction(model), -5.8514497)

  ## rw1 and rw2
  formula <- y ~ f(time, model = "rw1") + f(time2, model = "rw2")
  model <- INLA::inla(formula, data = data, control.compute = list(config = TRUE))
  expect_equal(log_mlik_correction(model), -5.8514497 - 3.45305292)
})
