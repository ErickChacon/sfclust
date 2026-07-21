library(stars)

# --- helpers -----------------------------------------------------------

make_geom_time <- function(ns, nt, dim_order = c("geometry", "time")) {
  space <- st_sfc(lapply(seq_len(ns), function(i) st_point(c(i, i))))
  time  <- seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
  vals  <- array(seq_len(ns * nt), dim = c(ns, nt))
  dims  <- if (identical(dim_order, c("geometry", "time"))) {
    st_dimensions(geometry = space, time = time)
  } else {
    vals <- t(vals)
    st_dimensions(time = time, geometry = space)
  }
  st_as_stars(val = vals, dimensions = dims)
}

make_raster <- function(nx, ny, nt, dim_order = c("x", "y", "time")) {
  time <- seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
  vals <- array(seq_len(nx * ny * nt), dim = c(nx, ny, nt))
  perm <- match(dim_order, c("x", "y", "time"))
  vals <- aperm(vals, perm)
  do.call(st_as_stars, c(
    list(val = vals),
    list(dimensions = do.call(st_dimensions, setNames(
      list(1:nx, 1:ny, time)[perm],
      dim_order
    )))
  ))
}

# --- data_all: vector geometry -----------------------------------------

test_that("data_all: geometry, time", {
  ns <- 10; nt <- 3
  stdata <- make_geom_time(ns, nt)
  df <- data_all(stdata)

  # dim order geometry(fast), time(slow)
  expect_equal(nrow(df),   ns * nt)
  expect_equal(df$id,      seq_len(ns * nt))
  expect_equal(df$ids,     rep(seq_len(ns), nt))
  expect_equal(df$id_time, rep(seq_len(nt), each = ns))
  expect_equal(df$val,     seq_len(ns * nt))
})

test_that("data_all: time, geometry", {
  ns <- 10; nt <- 3
  stdata <- make_geom_time(ns, nt, dim_order = c("time", "geometry"))
  df <- data_all(stdata)

  # dim order time(fast), geometry(slow)
  expect_equal(nrow(df),   ns * nt)
  expect_equal(df$id,      seq_len(ns * nt))
  expect_equal(df$ids,     rep(seq_len(ns), each = nt))
  expect_equal(df$id_time, rep(seq_len(nt), ns))
  expect_equal(df$val,     as.vector(t(array(seq_len(ns * nt), dim = c(ns, nt)))))
})

test_that("data_all: wavelength, time, geometry", {
  ns <- 5; nt <- 3; nw <- 2
  space <- st_sfc(lapply(seq_len(ns), function(i) st_point(c(i, i))))
  time  <- seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
  stdata <- st_as_stars(
    val = array(seq_len(nw * nt * ns), dim = c(nw, nt, ns)),
    dimensions = st_dimensions(wavelength = seq_len(nw), time = time, geometry = space)
  )
  df <- data_all(stdata, spnames = "geometry")

  # dim order wavelength(fast), time, geometry(slow)
  expect_equal(nrow(df),         ns * nt * nw)
  expect_equal(df$id,            seq_len(ns * nt * nw))
  expect_equal(df$ids,           rep(seq_len(ns), each = nt * nw))
  expect_equal(df$id_time,       rep(rep(seq_len(nt), each = nw), ns))
  expect_equal(df$id_wavelength, rep(seq_len(nw), ns * nt))
  expect_equal(df$val,           seq_len(ns * nt * nw))
})

test_that("data_all: vector geometry-only stars", {
  ns <- 5
  space <- st_sfc(lapply(seq_len(ns), function(i) st_point(c(i, i))))
  stdata <- st_as_stars(
    val = array(1:ns, dim = ns),
    dimensions = st_dimensions(geometry = space)
  )
  df <- data_all(stdata)

  expect_equal(nrow(df), ns)
  expect_equal(df$id,  seq_len(ns))
  expect_equal(df$ids, seq_len(ns))
  expect_equal(df$val, 1:ns)
  expect_equal(ncol(df[, grepl("^id_", names(df)), drop = FALSE]), 0L)
})

test_that("data_all: vector geometry with NA values", {
  ns <- 4; nt <- 3
  stdata <- make_geom_time(ns, nt)
  stdata$val[2, ] <- NA
  df <- data_all(stdata)

  # dim order geometry(fast), time(slow)
  expect_equal(nrow(df),   ns * nt)
  expect_equal(df$id,      seq_len(ns * nt))
  expect_equal(df$ids,     rep(seq_len(ns), nt))
  expect_equal(df$id_time, rep(seq_len(nt), each = ns))
  expect_equal(df$val,     as.vector(stdata$val))
})

# --- data_all: raster --------------------------------------------------

test_that("data_all: x, y, time", {
  nx <- 3; ny <- 2; nt <- 4
  stdata <- make_raster(nx, ny, nt)
  df <- data_all(stdata, spnames = c("x", "y"))
  df_yx <- data_all(stdata, spnames = c("y", "x"))

  # dim order x(fast), y, time(slow)
  expect_equal(nrow(df),    nx * ny * nt)
  expect_equal(df$id,       seq_len(nx * ny * nt))
  expect_equal(df$ids,      rep(seq_len(nx * ny), nt))
  expect_equal(df$id_time,  rep(seq_len(nt), each = nx * ny))
  expect_equal(df$val,      seq_len(nx * ny * nt))
  expect_equal(df,          df_yx)
})

test_that("data_all: x, y, time, wavelength", {
  nx <- 3; ny <- 2; nt <- 4; nw <- 2
  time <- seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
  vals <- array(seq_len(nx * ny * nt * nw), dim = c(nx, ny, nt, nw))
  stdata <- st_as_stars(
    val = vals,
    dimensions = st_dimensions(x = 1:nx, y = 1:ny, time = time, wavelength = seq_len(nw))
  )
  df <- data_all(stdata, spnames = c("x", "y"))
  df_yx <- data_all(stdata, spnames = c("y", "x"))

  # dim order x(fast), y, time, wavelength(slow)
  expect_equal(nrow(df),          nx * ny * nt * nw)
  expect_equal(df$id,             seq_len(nx * ny * nt * nw))
  expect_equal(df$ids,            rep(seq_len(nx * ny), nt * nw))
  expect_equal(df$id_time,        rep(rep(seq_len(nt), each = nx * ny), nw))
  expect_equal(df$id_wavelength,  rep(seq_len(nw), each = nx * ny * nt))
  expect_equal(df$val,            seq_len(nx * ny * nt * nw))
  expect_equal(df,                df_yx)
})

test_that("data_all: time, x, y", {
  nx <- 3; ny <- 2; nt <- 4
  stdata <- make_raster(nx, ny, nt, dim_order = c("time", "x", "y"))
  df <- data_all(stdata, spnames = c("x", "y"))
  df_yx <- data_all(stdata, spnames = c("y", "x"))

  # dim order time(fast), x, y(slow)
  expect_equal(nrow(df),   nx * ny * nt)
  expect_equal(df$id,      seq_len(nx * ny * nt))
  expect_equal(df$ids,     rep(seq_len(nx * ny), each = nt))
  expect_equal(df$id_time, rep(seq_len(nt), nx * ny))
  val <- as.vector(aperm(array(seq_len(nx * ny * nt), dim = c(nx, ny, nt)), c(3, 1, 2)))
  expect_equal(df$val,     val)
  expect_equal(df,         df_yx)
})

test_that("data_all: y, x, time", {
  nx <- 3; ny <- 2; nt <- 4
  stdata <- make_raster(nx, ny, nt, dim_order = c("y", "x", "time"))
  df <- data_all(stdata, spnames = c("y", "x"))
  df_xy <- data_all(stdata, spnames = c("x", "y"))

  # dim order y(fast), x, time(slow)
  expect_equal(nrow(df),   nx * ny * nt)
  expect_equal(df$id,      seq_len(nx * ny * nt))
  expect_equal(df$ids,     rep(seq_len(nx * ny), nt))
  expect_equal(df$id_time, rep(seq_len(nt), each = nx * ny))
  val <- as.vector(aperm(array(seq_len(nx * ny * nt), dim = c(nx, ny, nt)), c(2, 1, 3)))
  expect_equal(df$val,     val)
  expect_equal(df,         df_xy)
})

test_that("data_all: raster-only stars", {
  nx <- 3; ny <- 2
  stdata <- st_as_stars(
    val = matrix(seq_len(nx * ny), nx, ny),
    dimensions = st_dimensions(x = 1:nx, y = 1:ny)
  )
  df <- data_all(stdata, spnames = c("x", "y"))

  expect_equal(nrow(df), nx * ny)
  expect_equal(df$id,  seq_len(nx * ny))
  expect_equal(df$ids, seq_len(nx * ny))
  expect_equal(df$val, seq_len(nx * ny))
  expect_equal(ncol(df[, grepl("^id_", names(df)), drop = FALSE]), 0L)
})

test_that("data_all: raster NA values — interior cell at flat pos 2", {
  nx <- 3; ny <- 2; nt <- 2
  stdata <- make_raster(nx, ny, nt)
  stdata$val[2, 1, ] <- NA
  df <- data_all(stdata, spnames = c("x", "y"))

  # dim order x(fast), y, time(slow)
  expect_equal(nrow(df),   nx * ny * nt)
  expect_equal(df$id,      seq_len(nx * ny * nt))
  expect_equal(df$ids,     rep(seq_len(nx * ny), nt))
  expect_equal(df$id_time, rep(seq_len(nt), each = nx * ny))
  expect_equal(df$val,     as.vector(stdata$val))
})

# --- filter_df: vector geometry  ---------------------------------------

test_that("filter_df: vector subset a cluster", {
  ns <- 6; nt <- 3
  membership <- c(1, 1, 2, 2, 2, 3)
  stdata <- make_geom_time(ns, nt)

  # check domain
  domain <- create_domain(stdata, "geometry")
  expect_equal(as.vector(domain[[1]]), rep(TRUE, ns))

  # check filter
  df <- filter_df(data_all(stdata), which(domain[[1]]))
  cluster_units <- which(membership == 2)
  inla_data <- df[df$sid %in% cluster_units, , drop = FALSE]
  expect_equal(inla_data$sid,     rep(cluster_units, nt))
  expect_equal(inla_data$ids,     rep(cluster_units, nt))
  expect_equal(inla_data$id_time, rep(seq_len(nt), each = length(cluster_units)))
  expect_equal(inla_data$val,     as.vector(stdata$val[cluster_units, ]))
})

test_that("filter_df: vector geometry with NA region", {
  ns <- 4; nt <- 3
  stdata <- make_geom_time(ns, nt)
  stdata$val[2, ] <- NA  # region 2 is all NA

  # check domain: region 2 excluded
  domain <- create_domain(stdata, "geometry", "val")
  expect_equal(as.vector(domain[[1]]), c(TRUE, FALSE, TRUE, TRUE))

  # check filter: valid_ids = c(1, 3, 4), sid remapped to 1..3
  valid_ids <- which(domain[[1]])
  df <- filter_df(data_all(stdata), valid_ids)
  full_val <- as.vector(stdata$val)
  expect_equal(df$sid,     rep(seq_len(length(valid_ids)), nt))
  expect_equal(df$ids,     rep(valid_ids, nt))
  expect_equal(df$id_time, rep(seq_len(nt), each = length(valid_ids)))
  expect_equal(df$val,     full_val[!is.na(full_val)])
})

# --- filter_df: raster -------------------------------------------------

test_that("filter_df: raster subset a cluster", {
  nx <- 3; ny <- 2; nt <- 3
  membership <- c(1, 1, 2, 2, 2, 1)
  x <- make_raster(nx, ny, nt)

  # check domain: all cells valid
  domain <- create_domain(x, c("x", "y"))
  expect_equal(as.vector(domain[[1]]), rep(TRUE, nx * ny))

  # check filter
  df <- filter_df(data_all(x, c("x", "y")), which(domain[[1]]))
  cluster_units <- which(membership == 2)
  inla_data <- df[df$sid %in% cluster_units, , drop = FALSE]
  expect_equal(inla_data$sid,     rep(cluster_units, nt))
  expect_equal(inla_data$ids,     rep(cluster_units, nt))
  expect_equal(inla_data$id_time, rep(seq_len(nt), each = length(cluster_units)))
  expect_equal(inla_data$val,     df$val[df$sid %in% cluster_units])
})

test_that("filter_df: raster with NA cells", {
  nx <- 3; ny <- 2; nt <- 2
  membership <- c(1, 2, 1, 2, 1)
  x <- make_raster(nx, ny, nt)
  x$val[2, 1, ] <- NA

  # check domain: flat pos 2 (x=2, y=1) excluded
  domain <- create_domain(x, c("x", "y"), "val")
  expect_equal(as.vector(domain[[1]]), c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))

  # check filter
  valid_ids <- which(domain[[1]])
  df <- filter_df(data_all(x, c("x", "y")), valid_ids)
  full_val <- as.vector(x$val)
  expect_equal(df$sid,     rep(seq_len(length(valid_ids)), nt))
  expect_equal(df$ids,     rep(valid_ids, nt))
  expect_equal(df$id_time, rep(seq_len(nt), each = length(valid_ids)))
  expect_equal(df$val,     full_val[!is.na(full_val)])

  # sid remapping: cluster 1 sids c(1,3,5) map to flat ids c(1,4,6)
  cluster_units <- which(membership == 1)
  correct_flat_ids <- valid_ids[cluster_units]
  expect_equal(df$ids[df$sid %in% cluster_units], rep(correct_flat_ids, nt))
})
