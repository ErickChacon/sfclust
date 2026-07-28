library(stars)
library(Matrix)

# --- unique_clusters ---------------------------------------------------

test_that("unique_clusters: integer membership", {
  expect_equal(unique_clusters(c(1)), 1L)
  expect_equal(unique_clusters(c(1, 1, 1)), 1L)
  expect_equal(unique_clusters(c(2, 4, 4, 3, 1)), 1:4)
  expect_error(unique_clusters(c(3, 7, 7, 4, 2)))
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
  skip_if_not_installed("INLA")

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

# --- log_mlik_each INLA failure ----------------------------------------

test_that("log_mlik_each returns -Inf and warns when INLA fails", {
  skip_if_not_installed("INLA")
  data <- data.frame(y = rnorm(5), id_time = 1:5, ids = 1:5, id = 1:5, sid = 1:5)
  inla_args <- list(formula = y ~ f(nonexistent, model = "rw1"), family = "gaussian")
  expect_warning(
    result <- log_mlik_each(1, rep(1L, 5), data, FALSE, FALSE, inla_args),
    "INLA failed"
  )
  expect_equal(result, -Inf)
})

# --- log_mlik_all and log_mlik_each ------------------------------------

test_that("log_mlik_all: gaussian, binomial, poisson", {
  skip_if_not_installed("INLA")

  ns <- 3; nt <- 5
  membership <- c(1, 1, 2)
  df <- data.frame(
    sid     = rep(1:ns, nt),
    ids     = rep(1:ns, nt),
    id      = seq_len(ns * nt),
    id_time = rep(1:nt, each = ns),
    n       = c(8L, 10L, 12L, 9L, 11L, 10L, 8L, 9L, 12L, 10L, 11L, 8L, 9L, 10L, 12L),
    y_gaus  = c(1.37, -0.56, 0.36, 0.63, 0.40,  -0.11,  1.51, -0.09,
                2.02, -0.06, 1.30, 2.29, -1.39, -0.28, -0.13),
    y_bin   = c(4L, 4L, 3L, 4L, 0L, 4L, 0L, 2L, 5L, 3L, 2L, 3L, 1L, 6L, 3L),
    y_pois  = c(9L, 8L, 6L, 10L, 6L, 4L, 4L, 4L, 7L, 1L, 6L, 6L, 3L, 3L, 5L)
  )

  # Gaussian
  inla_args <- list(formula = y_gaus ~ f(id_time, model = "rw1"), family = "gaussian")
  res <- log_mlik_all(membership, df, correction = FALSE, inla_args = inla_args)
  expect_equal(res, c(-24.24213, -18.05626), tolerance = 1e-3)

  # Binomial
  inla_args <- list(
    formula = y_bin ~ f(id_time, model = "rw1"),
    family = "binomial", Ntrials = quote(n)
  )
  res <- log_mlik_all(membership, df, correction = FALSE, inla_args = inla_args)
  expect_equal(res, c(-22.44368, -8.355529), tolerance = 1e-3)

  # Poisson
  inla_args <- list(
    formula = y_pois ~ f(id_time, model = "rw1"),
    family = "poisson", E = quote(n)
  )
  res <- log_mlik_all(membership, df, correction = FALSE, inla_args = inla_args)
  expect_equal(res, c(-26.86092, -10.86669), tolerance = 1e-3)
})

