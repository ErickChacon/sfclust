library(stars)

test_that('convert stars object to long format', {

  space <- st_sfc(lapply(1:10, function(i) st_point(c(i, i))))
  time <- seq(as.Date("2024-10-01"), by = "1 day", length.out = 3)
  ns <- length(space)
  nt <- length(time)
  membership <- rep(1:3, length = ns)
  k <- 1
  nk <- sum(membership == k)

  # space dimension on rows and time dimension on columns
  stdata <- st_as_stars(
    cases = array(1:(ns * nt), dim = c(ns, nt)),
    dimensions = st_dimensions(geometry = space, time = time)
  )

  stdata_k <- data_each(k = k, membership, stdata)
  expect_equal(nrow(stdata_k), nk * nt)
  expect_equal(stdata_k$id, as.numeric(outer(which(membership == k), ns * (1:nt - 1), `+`)))
  expect_equal(stdata_k$ids, rep(which(membership == k), nt))
  expect_equal(stdata_k$idt, rep(1:nt, each = nk))
  expect_equal(stdata_k$cases, as.numeric(outer(which(membership == k), ns * (1:nt - 1), `+`)))

  stdata_long <- data_all(stdata)
  expect_equal(stdata_long$id, 1:(ns*nt))
  expect_equal(stdata_long$ids, rep(1:ns, nt))
  expect_equal(stdata_long$idt, rep(1:nt, each = ns))
  expect_equal(stdata_long$time, rep(time, each = ns))
  expect_equal(stdata_long$cases, 1:(ns*nt))

  ## additional third dimension
  stdata <- st_as_stars(
    cases = array(1:(ns * nt), dim = c(ns, nt, 1)),
    dimensions = st_dimensions(geometry = space, time = time, band = 1, point = TRUE)
  )

  xk_aux <- data_each(k = k, membership, stdata)
  expect_equal(stdata_k, subset(xk_aux, select = - band))
  expect_equal(stdata_long, subset(data_all(stdata), select = - band))

  # time dimension on rows and space dimension on columns
  stdata <- st_as_stars(
    cases = t(array(1:(ns * nt), dim = c(ns, nt))),
    dimensions = st_dimensions(time = time, geometry = space)
  )

  stdata_k <- data_each(k = k, membership, stdata)
  expect_equal(nrow(stdata_k), nk * nt)
  expect_equal(stdata_k$id, as.numeric(outer(1:nt, nt * (which(membership == k) - 1), `+`)))
  expect_equal(stdata_k$ids, rep(which(membership == k), each = nt))
  expect_equal(stdata_k$idt, rep(1:nt, nk))
  expect_equal(stdata_k$cases, as.numeric(outer(ns * (1:nt-1), which(membership == k), `+`)))

  stdata_long <- data_all(stdata)
  expect_equal(stdata_long$id, 1:(ns*nt))
  expect_equal(stdata_long$ids, rep(1:ns, each = nt))
  expect_equal(stdata_long$idt, rep(1:nt, ns))
  expect_equal(stdata_long$time, rep(time, ns))
  expect_equal(stdata_long$cases, as.numeric(outer(ns * (1:nt-1), 1:ns, `+`)))
})

test_that('compute log marginal correction', {
  skip_on_cran()

  # terms that require correction
  formula <- y ~ x + z
  expect_equal(correction_required(formula), character())

  formula <- y ~ x + f(z, model = "rw1")
  expect_equal(correction_required(formula), "z")

  formula <- y ~ f(z, model="rw2") + x
  expect_equal(correction_required(formula), "z")

  formula <- y ~ x +
      f(z, model = "crw1", hyper = list(prec = list(prior = "loggamma", param = c(1, 0.01)))) +
      f(w, model = "rw2") +
      f(v, model = "ar")
  expect_equal(correction_required(formula), c("z", "w"))

  # structure matrix and log marginal likelihood correction
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
  formula <- y ~ f(time, model = "rw1") +  f(time2, model = "rw2")
  model <- INLA::inla(formula, data = data, control.compute = list(config = TRUE))
  expect_equal(log_mlik_correction(model), - 5.8514497 - 3.45305292)
})

test_that('obtain unique clusters from membership', {
  membership <- c(2, 4, 4, 3, 1)
  expect_equal(unique_clusters(membership), 1:4)

  membership <- c(3, 7, 7, 4, 2)
  expect_error(unique_clusters(membership))

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
