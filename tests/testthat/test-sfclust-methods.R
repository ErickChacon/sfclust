library(sf)
library(stars)

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

# --- shared fixtures ---------------------------------------------------

make_sfclust_vec <- function() {
  geom <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))
  time <- seq(as.Date("2024-01-01"), by = "1 day", length.out = 4)
  ns <- length(geom); nt <- length(time)
  x <- st_as_stars(
    y = array(rnorm(ns * nt), dim = c(ns, nt)),
    dimensions = st_dimensions(geometry = geom, time = time)
  )
  set.seed(42)
  gd <- genclust(x, nclust = 2)
  suppressWarnings(
    sfclust(x, graphdata = gd, formula = y ~ 1, family = "gaussian",
            niter = 3, burnin = 0, thin = 1, nmessage = 1)
  )
}

make_sfclust_vec_binom <- function() {
  geom <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(3, 2))
  time <- seq(as.Date("2024-01-01"), by = "1 day", length.out = 4)
  ns <- length(geom); nt <- length(time)
  population <- array(rep(100L, ns * nt), dim = c(ns, nt))
  x <- st_as_stars(
    cases      = array(rbinom(ns * nt, 100L, 0.3), dim = c(ns, nt)),
    population = population,
    dimensions = st_dimensions(geometry = geom, time = time)
  )
  set.seed(42)
  gd <- genclust(x, nclust = 2)
  suppressWarnings(
    sfclust(x, graphdata = gd, formula = cases ~ 1, family = "binomial",
            Ntrials = population, niter = 3, burnin = 0, thin = 1, nmessage = 1)
  )
}

make_sfclust_raster_na <- function() {
  nx <- 3; ny <- 2; nt <- 4
  vals <- array(rnorm(nx * ny * nt), dim = c(nx, ny, nt))
  vals[2, 1, ] <- NA
  x <- st_as_stars(
    z = vals,
    dimensions = st_dimensions(
      x    = 1:nx, y = 1:ny,
      time = seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
    )
  )
  set.seed(42)
  gd <- genclust(x, spnames = c("x", "y"), response = "z", nclust = 2)
  suppressWarnings(
    sfclust(x, graphdata = gd, spnames = c("x", "y"),
            formula = z ~ 1, family = "gaussian",
            niter = 3, burnin = 0, thin = 1, nmessage = 1)
  )
}

# --- print / summary ---------------------------------------------------

test_that("print.sfclust returns object invisibly and shows formula", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  out <- capture.output(ret <- print(result))
  expect_identical(ret, result)
  expect_true(any(grepl("formula", out, ignore.case = TRUE)))
})

test_that("summary.sfclust returns table summing to ns", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  s <- summary(result)
  expect_s3_class(s, "table")
  expect_equal(sum(s), 6L)
})

test_that("summary.sfclust: sort relabels by size", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  s <- summary(result, sort = TRUE)
  expect_true(s[[1]] >= s[[length(s)]])
})

# --- fitted ------------------------------------------------------------

test_that("fitted.sfclust returns data frame with cluster and mean columns", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  df <- sfclust:::fitted.sfclust(result)
  expect_true(is.data.frame(df))
  expect_true(all(c("cluster", "mean") %in% names(df)))
  expect_equal(nrow(df), 6L * 4L)
})

test_that("fitted.sfclust_stars returns stars object covering full grid", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  s <- fitted(result)
  expect_s3_class(s, "stars")
  expect_equal(prod(dim(s)), 6L * 4L)
})

test_that("fitted.sfclust_stars: raster with NA cells has NA for invalid cell", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_raster_na()
  s <- fitted(result)
  expect_s3_class(s, "stars")
  expect_equal(prod(dim(s)), 3L * 2L * 4L)
  df <- as.data.frame(s)
  expect_true(any(is.na(df$mean)))   # flat position 2 (NA cell) is NA in output
})

# --- column-reference INLA args (Ntrials, E) ---------------------------

test_that("sfclust with Ntrials column reference fits without error", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec_binom()
  expect_s3_class(result, "sfclust_stars")
  expect_equal(nrow(result$samples$membership), 3L)
})

test_that("update.sfclust with Ntrials column reference continues chain", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec_binom()
  result2 <- suppressWarnings(update(result, niter = 2, nmessage = 1))
  expect_s3_class(result2, "sfclust_stars")
  expect_equal(nrow(result2$samples$membership), 2L)
})

# --- inla_args stored --------------------------------------------------

test_that("inla_args stores formula and family from sfclust call", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec_binom()
  ia <- attr(result, "inla_args")
  expect_equal(eval(ia$formula), cases ~ 1)
  expect_equal(eval(ia$family), "binomial")
})

# --- update ------------------------------------------------------------

test_that("update.sfclust continues chain from last state", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  result2 <- suppressWarnings(update(result, niter = 2, nmessage = 1))
  expect_s3_class(result2, "sfclust_stars")
  expect_equal(nrow(result2$samples$membership), 2L)
})

test_that("update.sfclust with sample refits within-cluster models", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  result2 <- update(result, sample = 1)
  expect_equal(result2$clust$id, 1L)
  expect_equal(length(result2$clust$models), max(result$samples$membership[1, ]))
})

# --- plots -------------------------------------------------------------

test_that("plot_log_mlik returns ggplot", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  gg <- plot_log_mlik(result)
  expect_s3_class(gg, "ggplot")
})

test_that("plot_clusters_fitted returns ggplot", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  gg <- plot_clusters_fitted(result)
  expect_s3_class(gg, "ggplot")
})

test_that("plot_clusters_series returns ggplot", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  gg <- plot_clusters_series(result, var = y)
  expect_s3_class(gg, "ggplot")
})

test_that("plot_clusters_map returns ggplot", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_vec()
  gg <- plot_clusters_map(result)
  expect_s3_class(gg, "ggplot")
})

test_that("plot_clusters_map: raster with NAs returns ggplot", {
  skip_if_not_installed("INLA")
  result <- make_sfclust_raster_na()
  gg <- plot_clusters_map(result)
  expect_s3_class(gg, "ggplot")
})
