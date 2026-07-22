library(sf)
library(stars)
library(igraph)
library(Matrix)

# --- sort_membership ---------------------------------------------------

test_that("sort_membership: relabels by cluster size", {
  mem <- c(4, 2, 4, 1, 3, 3, 1, 3)
  mem_new <- sort_membership(mem)
  expect_equal(as.integer(mem_new), c(3, 4, 3, 2, 1, 1, 2, 1))
  expect_equal(attr(mem_new, "order"), c(3, 1, 4, 2))
})

test_that("sort_membership: identity when already sorted by size", {
  mem <- c(4, 1, 1, 3, 2, 1, 2, 3)
  mem_new <- sort_membership(mem)
  expect_equal(as.integer(mem_new), mem)
  expect_equal(attr(mem_new, "order"), 1:4)
})

# --- sfclust: vector geometry, gaussian --------------------------------

test_that("sfclust_stars: vector geometry, gaussian", {
  skip_if_not_installed("INLA")
  nx <- 3; ny <- 2; ns <- nx * ny; nt <- 4
  geom <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(nx, ny))
  time <- seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
  x <- st_as_stars(
    y = array(rnorm(ns * nt), dim = c(ns, nt)),
    dimensions = st_dimensions(geometry = geom, time = time)
  )

  set.seed(1)
  result <- suppressWarnings(
    sfclust(x, nclust = 2, formula = y ~ 1, family = "gaussian",
            niter = 3, burnin = 0, thin = 1, nmessage = 1)
  )

  test_that("structure", {
    expect_s3_class(result, c("sfclust_stars", "sfclust"))
    expect_named(result, c("samples", "clust"))
    expect_equal(dim(result$samples$membership), c(3L, 6L))
    expect_equal(length(result$samples$log_mlike), 3L)
    expect_named(result$samples$move_counts, c("births", "deaths", "changes", "hypers"))
    expect_equal(length(result$clust$membership), 6L)
    expect_equal(length(result$clust$models), max(result$clust$membership))
    input_args <- attr(result, "input_args")
    expect_s3_class(input_args$stars, "stars")
    expect_equal(input_args$spnames, "geometry")
    expect_equal(input_args$fnames, "time")
    fit_args <- attr(result, "fit_args")
    expect_equal(fit_args$niter, 3L)
    expect_equal(fit_args$burnin, 0L)
    expect_equal(fit_args$thin, 1L)
    expect_equal(fit_args$correction, FALSE)
  })

  test_that("print and summary", {
    out <- capture.output(ret <- print(result))
    expect_identical(ret, result)
    s <- summary(result)
    expect_s3_class(s, "table")
    expect_equal(sum(s), 6L)
    s_sorted <- summary(result, sort = TRUE)
    expect_true(s_sorted[[1]] >= s_sorted[[length(s_sorted)]])
  })

  test_that("fitted", {
    s <- fitted(result)
    expect_s3_class(s, "stars")
    expect_equal(dim(s), c(geometry = 6L, time = 4L))
  })

  test_that("update", {
    result2 <- suppressWarnings(update(result, niter = 2, nmessage = 1))
    expect_s3_class(result2, "sfclust_stars")
    expect_equal(nrow(result2$samples$membership), 2L)
    result3 <- update(result, sample = 1)
    expect_equal(result3$clust$id, 1L)
    expect_equal(length(result3$clust$models), max(result$samples$membership[1, ]))
  })

  test_that("plots", {
    expect_s3_class(plot_log_mlik(result), "ggplot")
    expect_s3_class(plot_clusters_fitted(result), "ggplot")
    expect_s3_class(plot_clusters_series(result, var = y), "ggplot")
    expect_s3_class(plot_clusters_map(result), "ggplot")
  })
})

# --- sfclust: vector geometry, binomial --------------------------------

test_that("sfclust_stars: vector geometry, binomial", {
  skip_if_not_installed("INLA")
  nx <- 3; ny <- 2; ns <- nx * ny; nt <- 4
  geom <- st_make_grid(cellsize = c(1, 1), offset = c(0, 0), n = c(nx, ny))
  time <- seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
  x <- st_as_stars(
    cases      = array(rbinom(ns * nt, 100L, 0.3), dim = c(ns, nt)),
    population =array(rep(100L, ns * nt), dim = c(ns, nt)),
    dimensions = st_dimensions(geometry = geom, time = time)
  )
  set.seed(2)
  result <- suppressWarnings(
    sfclust(x, nclust = 3, formula = cases ~ 1, family = "binomial",
            Ntrials = population, niter = 3, burnin = 0, thin = 1, nmessage = 1)
  )

  test_that("structure and inla_args", {
    expect_s3_class(result, c("sfclust_stars", "sfclust"))
    expect_equal(dim(result$samples$membership), c(3L, 6L))
    inla_args <- attr(result, "inla_args")
    expect_equal(eval(inla_args$formula), cases ~ 1)
    expect_equal(eval(inla_args$family), "binomial")
    expect_equal(inla_args$Ntrials, quote(population))
  })

  test_that("print, summary, fitted", {
    expect_no_error(print(result))
    expect_no_error(summary(result))
    s <- fitted(result)
    expect_s3_class(s, "stars")
    expect_equal(dim(s), c(geometry = ns, time = nt))
  })

  test_that("update", {
    result2 <- suppressWarnings(update(result, niter = 2, nmessage = 1))
    expect_s3_class(result2, "sfclust_stars")
    expect_equal(nrow(result2$samples$membership), 2L)
  })

  test_that("plots", {
    expect_s3_class(plot_log_mlik(result), "ggplot")
    expect_s3_class(plot_clusters_fitted(result), "ggplot")
    expect_s3_class(plot_clusters_series(result, var = cases), "ggplot")
    expect_s3_class(plot_clusters_map(result), "ggplot")
  })
})

# --- sfclust: raster, poisson ------------------------------------------

test_that("sfclust_stars: raster, poisson", {
  skip_if_not_installed("INLA")
  nx <- 3; ny <- 2; ns <- nx * ny; nt <- 4
  x <- st_as_stars(
    counts   = array(rpois(nx * ny * nt, lambda = 5), dim = c(nx, ny, nt)),
    expected = array(rep(4L, nx * ny * nt), dim = c(nx, ny, nt)),
    dimensions = st_dimensions(
      x    = 1:nx, y = 1:ny,
      time = seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
    )
  )

  set.seed(5)
  result <- suppressWarnings(
    sfclust(x, nclust = 3, spnames = c("x", "y"),
            formula = counts ~ 1, family = "poisson", E = expected,
            niter = 3, burnin = 0, thin = 1, nmessage = 1)
  )

  test_that("structure and inla_args", {
    expect_s3_class(result, c("sfclust_stars", "sfclust"))
    expect_equal(dim(result$samples$membership), c(3L, 6L))
    expect_equal(attr(result, "input_args")$spnames, c("x", "y"))
    inla_args <- attr(result, "inla_args")
    expect_equal(eval(inla_args$formula), counts ~ 1)
    expect_equal(eval(inla_args$family), "poisson")
    expect_equal(inla_args$E, quote(expected))
  })

  test_that("print, summary, fitted", {
    expect_no_error(print(result))
    expect_no_error(summary(result))
    s <- fitted(result)
    expect_s3_class(s, "stars")
    expect_equal(dim(s), c(x = nx, y = ny, time = nt))
  })

  test_that("update and plots", {
    result2 <- suppressWarnings(update(result, niter = 2, nmessage = 1))
    expect_s3_class(result2, "sfclust_stars")
    expect_equal(nrow(result2$samples$membership), 2L)
    expect_s3_class(plot_log_mlik(result), "ggplot")
    expect_s3_class(plot_clusters_fitted(result), "ggplot")
    expect_s3_class(plot_clusters_series(result, var = counts), "ggplot")
    expect_s3_class(plot_clusters_map(result), "ggplot")
  })
})

# --- sfclust: raster with NA cells -------------------------------------

test_that("sfclust_stars: raster with NA cells", {
  skip_if_not_installed("INLA")
  nx <- 5; ny <- 4; nt <- 4
  vals <- array(rnorm(nx * ny * nt), dim = c(nx, ny, nt))
  vals[1, c(1, 2, 4), ] <- NA
  vals[2, 1, ]          <- NA
  vals[4, c(1, 2), ]    <- NA
  vals[5, c(1, 2), ]    <- NA
  x <- st_as_stars(
    z = vals,
    dimensions = st_dimensions(
      x    = 1:nx, y = 1:ny,
      time = seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
    )
  )

  set.seed(3)
  result <- suppressWarnings(
    sfclust(x, nclust = 4, spnames = c("x", "y"),
            formula = z ~ 1, family = "gaussian",
            niter = 3, burnin = 0, thin = 1, nmessage = 1)
  )

  test_that("structure", {
    expect_s3_class(result, c("sfclust_stars", "sfclust"))
    expect_equal(dim(result$samples$membership), c(3L, 12L))
    input_args <- attr(result, "input_args")
    expect_equal(input_args$spnames, c("x", "y"))
    expect_equal(input_args$fnames, "time")
  })

  test_that("print, summary, fitted", {
    expect_no_error(print(result))
    s <- summary(result)
    expect_s3_class(s, "table")
    expect_equal(sum(s), 12L)
    fitted_s <- fitted(result)
    expect_s3_class(fitted_s, "stars")
    expect_equal(dim(fitted_s), c(x = 5L, y = 4L, time = 4L))
    expect_equal(is.na(fitted_s[["mean"]]), is.na(x[["z"]]))
  })

  test_that("update and plots", {
    result2 <- suppressWarnings(update(result, niter = 2, nmessage = 1))
    expect_s3_class(result2, "sfclust_stars")
    expect_equal(nrow(result2$samples$membership), 2L)
    expect_s3_class(plot_log_mlik(result), "ggplot")
    expect_s3_class(plot_clusters_fitted(result), "ggplot")
    expect_s3_class(plot_clusters_series(result, var = z), "ggplot")
    expect_s3_class(plot_clusters_map(result), "ggplot")
  })
})

# --- sfclust: data.frame -----------------------------------------------

test_that("sfclust: data.frame", {
  skip_if_not_installed("INLA")
  ns <- 6L; nt <- 4L
  set.seed(4)
  df  <- data.frame(
    id      = seq_len(ns * nt),
    ids     = rep(seq_len(ns), nt),
    id_time = rep(seq_len(nt), each = ns),
    y       = rnorm(ns * nt, mean = rep(c(-1, 1), each = ns / 2, nt))
  )
  adj <- sparseMatrix(
    i = c(1, 2, 4, 5, 1, 2, 3), j = c(2, 3, 5, 6, 4, 5, 6),
    x = 1L, dims = c(ns, ns), symmetric = TRUE
  )

  result <- suppressWarnings(
    sfclust(df, adjacency = adj, nclust = 4, fnames = "id_time",
            formula = y ~ 1, family = "gaussian",
            niter = 3, burnin = 0, thin = 1, nmessage = 1)
  )

  test_that("structure", {
    expect_s3_class(result, "sfclust")
    expect_false(inherits(result, "sfclust_stars"))
    expect_named(result, c("samples", "clust"))
    expect_equal(dim(result$samples$membership), c(3L, 6L))
    expect_equal(length(result$samples$log_mlike), 3L)
    expect_named(result$samples$move_counts, c("births", "deaths", "changes", "hypers"))
    expect_equal(length(result$clust$membership), 6L)
    expect_equal(attr(result, "input_args")$fnames, "id_time")
  })

  test_that("print, summary, fitted", {
    expect_no_error(print(result))
    expect_no_error(summary(result))
    df_fitted <- fitted(result)
    expect_true(is.data.frame(df_fitted))
    expect_true(all(c("cluster", "mean") %in% names(df_fitted)))
    expect_equal(nrow(df_fitted), 6L * 4L)
  })

  test_that("update and plots", {
    result2 <- suppressWarnings(update(result, niter = 2, nmessage = 1))
    expect_s3_class(result2, "sfclust")
    expect_equal(nrow(result2$samples$membership), 2L)
    expect_s3_class(plot_log_mlik(result), "ggplot")
    expect_s3_class(plot_clusters_fitted(result), "ggplot")
    expect_s3_class(plot_clusters_series(result, var = y), "ggplot")
  })
})
