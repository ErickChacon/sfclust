library(stars)

test_that('filter stars object and convert to long format', {
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
  expect_equal(stdata_k$ids, rep(1:nk, nt))
  expect_equal(stdata_k$idt, rep(1:nt, each = nk))
  expect_equal(stdata_k$cases, as.numeric(outer(c(1, 4, 7, 10), 10 * (0:2), `+`)))

  ## additional third dimension
  stdata <- st_as_stars(
    cases = array(1:(ns * nt), dim = c(ns, nt, 1)),
    dimensions = st_dimensions(geometry = space, time = time, band = 1, point = TRUE)
  )

  xk_aux <- data_each(k = k, membership, stdata)
  expect_equal(stdata_k, subset(xk_aux, select = - band))

  # time dimension on rows and space dimension on columns
  stdata <- st_as_stars(
    cases = t(array(1:(ns * nt), dim = c(ns, nt))),
    dimensions = st_dimensions(time = time, geometry = space)
  )

  stdata_k <- data_each(k = k, membership, stdata)
  expect_equal(nrow(stdata_k), nk * nt)
  expect_equal(stdata_k$ids, rep(1:nk, each = nt))
  expect_equal(stdata_k$idt, rep(1:nt, nk))
  expect_equal(stdata_k$cases, as.numeric(outer(10 * (0:2), c(1, 4, 7, 10), `+`)))
})

test_that('get structure matrix and log marginal correction', {
})

test_that('log marginal likelihood', {
})



test_that('get structure matrix', {
})
