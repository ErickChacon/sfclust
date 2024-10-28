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
  x <- st_as_stars(
    cases = array(1:(ns * nt), dim = c(ns, nt)),
    dimensions = st_dimensions(geometry = space, time = time)
  )

  xk <- data_each(x, k = k, membership)
  expect_equal(nrow(xk), nk * nt)
  expect_equal(xk$ids, rep(1:nk, nt))
  expect_equal(xk$idt, rep(1:nt, each = nk))
  expect_equal(xk$cases, as.numeric(outer(c(1, 4, 7, 10), 10 * (0:2), `+`)))

  ## additional third dimension
  x <- st_as_stars(
    cases = array(1:(ns * nt), dim = c(ns, nt, 1)),
    dimensions = st_dimensions(geometry = space, time = time, band = 1, point = TRUE)
  )

  xk_aux <- data_each(x, k = k, membership)
  expect_equal(xk, subset(xk_aux, select = - band))

  # time dimension on rows and space dimension on columns
  x <- st_as_stars(
    cases = t(array(1:(ns * nt), dim = c(ns, nt))),
    dimensions = st_dimensions(time = time, geometry = space)
  )

  xk <- data_each(x, k = k, membership)
  expect_equal(nrow(xk), nk * nt)
  expect_equal(xk$ids, rep(1:nk, each = nt))
  expect_equal(xk$idt, rep(1:nt, nk))
  expect_equal(xk$cases, as.numeric(outer(10 * (0:2), c(1, 4, 7, 10), `+`)))
})
