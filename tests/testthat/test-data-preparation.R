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
  st_as_stars(cases = vals, dimensions = dims)
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

test_that("data_all: geometry, time — ids, id_time, id", {
  ns <- 10; nt <- 3
  stdata <- make_geom_time(ns, nt)
  df <- data_all(stdata)

  expect_equal(nrow(df), ns * nt)
  expect_equal(df$id,       seq_len(ns * nt))
  expect_equal(df$ids,      rep(seq_len(ns), nt))
  expect_equal(df$id_time, rep(seq_len(nt), each = ns))
  expect_equal(df$cases,    seq_len(ns * nt))
})

test_that("data_all: time, geometry — ids, id_time, id differ from geometry-first", {
  ns <- 10; nt <- 3
  st_gt <- make_geom_time(ns, nt)
  st_tg <- make_geom_time(ns, nt, dim_order = c("time", "geometry"))

  df_gt <- data_all(st_gt)
  df_tg <- data_all(st_tg)

  # ids and id_time should match when sorted by value (geometry index)
  expect_equal(df_gt$ids[order(df_gt$cases)], df_tg$ids[order(df_tg$cases)])
  expect_equal(df_gt$id_time[order(df_gt$cases)], df_tg$id_time[order(df_tg$cases)])

  # row order differs between layouts (different dim varies fastest)
  expect_false(identical(df_gt$cases, df_tg$cases))
})

test_that("data_all: wavelength, time, geometry — obs_dims derived automatically", {
  ns <- 5; nt <- 3; nw <- 2
  space <- st_sfc(lapply(seq_len(ns), function(i) st_point(c(i, i))))
  time  <- seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
  vals  <- array(seq_len(nw * nt * ns), dim = c(nw, nt, ns))
  stdata <- st_as_stars(
    val = vals,
    dimensions = st_dimensions(wavelength = seq_len(nw), time = time, geometry = space)
  )

  df <- data_all(stdata, spnames = "geometry")

  expect_equal(nrow(df), ns * nt * nw)
  expect_equal(sort(unique(df$ids)),        seq_len(ns))
  expect_equal(sort(unique(df$id_time)),   seq_len(nt))
  expect_equal(sort(unique(df$id_wavelength)), seq_len(nw))
  # id is the original flat array position — must span 1..n
  expect_equal(sort(df$id), seq_len(ns * nt * nw))
})

# --- data_all: raster --------------------------------------------------

test_that("data_all: x, y, time — ids column-major over x then y", {
  nx <- 3; ny <- 2; nt <- 4
  stdata <- make_raster(nx, ny, nt)
  df <- data_all(stdata, spnames = c("x", "y"))

  expect_equal(nrow(df), nx * ny * nt)
  expected_ids <- rep((rep(seq_len(ny), each = nx) - 1L) * nx + rep(seq_len(nx), ny), nt)
  expect_equal(df$ids, expected_ids)
  expect_equal(sort(unique(df$id_time)), seq_len(nt))
  expect_equal(sort(df$id), seq_len(nx * ny * nt))
})

test_that("data_all: x, y, time, wavelength — two obs_dims", {
  nx <- 3; ny <- 2; nt <- 4; nw <- 2
  time <- seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
  vals <- array(seq_len(nx * ny * nt * nw), dim = c(nx, ny, nt, nw))
  stdata <- st_as_stars(
    val = vals,
    dimensions = st_dimensions(x = 1:nx, y = 1:ny, time = time, wavelength = seq_len(nw))
  )

  df <- data_all(stdata, spnames = c("x", "y"))

  expect_equal(nrow(df), nx * ny * nt * nw)
  expect_equal(sort(unique(df$ids)),            seq_len(nx * ny))
  expect_equal(sort(unique(df$id_time)),        seq_len(nt))
  expect_equal(sort(unique(df$id_wavelength)),  seq_len(nw))
  expect_equal(sort(df$id), seq_len(nx * ny * nt * nw))
})

test_that("data_all: ids consistent across dimension orderings (x before y in stars dims)", {
  nx <- 3; ny <- 2; nt <- 4
  # x, y, time / time, x, y / x, time, y — all have x before y in stars dims
  st_xyt <- make_raster(nx, ny, nt, dim_order = c("x", "y", "time"))
  st_txy <- make_raster(nx, ny, nt, dim_order = c("time", "x", "y"))
  st_xty <- make_raster(nx, ny, nt, dim_order = c("x", "time", "y"))

  df_xyt <- data_all(st_xyt, spnames = c("x", "y"))
  df_txy <- data_all(st_txy, spnames = c("x", "y"))
  df_xty <- data_all(st_xty, spnames = c("x", "y"))

  # sort by val (unique 1..n) to get canonical order for comparison
  expect_equal(df_xyt$ids[order(df_xyt$val)],     df_txy$ids[order(df_txy$val)])
  expect_equal(df_xyt$ids[order(df_xyt$val)],     df_xty$ids[order(df_xty$val)])
  expect_equal(df_xyt$id_time[order(df_xyt$val)], df_txy$id_time[order(df_txy$val)])
  expect_equal(df_xyt$id_time[order(df_xyt$val)], df_xty$id_time[order(df_xty$val)])

  # row order differs between layouts (different dim varies fastest)
  expect_false(identical(df_xyt$val, df_txy$val))
  expect_false(identical(df_xyt$val, df_xty$val))
})

test_that("data_all: time, y, x — spnames order-invariant within that object", {
  nx <- 3; ny <- 2; nt <- 4
  # y appears before x in stars dims -> canonical spnames normalizes to c("y","x")
  st_tyx <- make_raster(nx, ny, nt, dim_order = c("time", "y", "x"))

  df_xy <- data_all(st_tyx, spnames = c("x", "y"))
  df_yx <- data_all(st_tyx, spnames = c("y", "x"))

  expect_equal(df_xy, df_yx)
})

test_that("data_all: spnames order-invariant — c('x','y') same as c('y','x')", {
  nx <- 3; ny <- 2; nt <- 4
  stdata <- make_raster(nx, ny, nt)

  df_xy <- data_all(stdata, spnames = c("x", "y"))
  df_yx <- data_all(stdata, spnames = c("y", "x"))

  expect_equal(df_xy, df_yx)
})

test_that("data_all: spatial-only stars (no functional dimensions)", {
  ns <- 5
  space <- st_sfc(lapply(seq_len(ns), function(i) st_point(c(i, i))))
  stdata <- st_as_stars(list(val = array(1:ns, dim = ns)),
                        dimensions = st_dimensions(geometry = space))
  df <- data_all(stdata)

  expect_equal(nrow(df), ns)
  expect_equal(df$id,  seq_len(ns))
  expect_equal(df$ids, seq_len(ns))
  expect_equal(df$val, 1:ns)
  expect_equal(ncol(df[, grepl("^id_", names(df)), drop = FALSE]), 0L)
})

test_that("data_all: raster spatial-only stars (no functional dimensions)", {
  nx <- 3; ny <- 2
  stdata <- st_as_stars(
    val = matrix(seq_len(nx * ny), nx, ny),
    dimensions = st_dimensions(x = 1:nx, y = 1:ny)
  )
  df <- data_all(stdata, spnames = c("x", "y"))

  expect_equal(nrow(df), nx * ny)
  expect_equal(sort(unique(df$ids)), seq_len(nx * ny))
  expect_equal(ncol(df[, grepl("^id_", names(df)), drop = FALSE]), 0L)
})

test_that("data_all: no filtering — NA cells are kept", {
  nx <- 3; ny <- 2; nt <- 2
  vals <- array(seq_len(nx * ny * nt), dim = c(nx, ny, nt))
  vals[1, 1, ] <- NA
  stdata <- st_as_stars(
    val = vals,
    dimensions = st_dimensions(x = 1:nx, y = 1:ny,
      time = seq(as.Date("2024-01-01"), by = "1 day", length.out = nt))
  )

  df <- data_all(stdata, spnames = c("x", "y"))

  # all rows returned including NA cells
  expect_equal(nrow(df), nx * ny * nt)
  expect_true(any(is.na(df$val)))
  # ids still span full 1..nx*ny
  expect_equal(sort(unique(df$ids)), seq_len(nx * ny))
})

# --- filter_df ---------------------------------------------------------

test_that("filter_df: sid correctly selects cluster rows for vector geometry", {
  ns <- 6; nt <- 3
  space <- st_sfc(lapply(seq_len(ns), function(i) st_point(c(i, i))))
  time  <- seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
  stdata <- st_as_stars(
    cases = array(rpois(ns * nt, 5), dim = c(ns, nt)),
    dimensions = st_dimensions(geometry = space, time = time)
  )
  domain <- create_domain(stdata, "geometry")
  df <- filter_df(data_all(stdata), which(domain[[1]]))
  membership <- c(1, 1, 2, 2, 2, 3)

  # no NAs: sid == ids == 1..ns; cluster 2 has sid 3, 4, 5
  cluster_units <- which(membership == 2)
  inla_data <- df[df$sid %in% cluster_units, , drop = FALSE]

  expect_equal(sort(unique(inla_data$sid)), 3:5)
  expect_equal(sort(unique(inla_data$ids)), 3:5)
  expect_equal(nrow(inla_data), 3 * nt)
})

test_that("filter_df: sid remapped to 1..n_valid for raster with NA cells", {
  # 3x2 raster where flat position 2 (cell x=2, y=1) is NA
  # valid flat positions: 1, 3, 4, 5, 6  (n_valid = 5)
  nx <- 3; ny <- 2; nt <- 2
  vals <- array(seq_len(nx * ny * nt), dim = c(nx, ny, nt))
  vals[2, 1, ] <- NA
  x <- st_as_stars(
    val = vals,
    dimensions = st_dimensions(
      x = 1:nx, y = 1:ny,
      time = seq(as.Date("2024-01-01"), by = "1 day", length.out = nt)
    )
  )

  domain <- create_domain(x, c("x", "y"), "val")
  df <- filter_df(data_all(x, c("x", "y")), which(domain[[1]]))

  # valid_ids maps valid unit index (1..5) -> flat position
  valid_ids <- which(domain[[1]])   # c(1, 3, 4, 5, 6)

  # membership over n_valid = 5 units
  membership <- c(1, 2, 1, 2, 1)

  # correct flat positions for cluster 1: valid units 1, 3, 5 -> flat 1, 4, 6
  correct_flat_ids <- valid_ids[which(membership == 1)]  # c(1, 4, 6)

  # sid = 1..n_valid; which(membership == 1) = c(1, 3, 5) matches sid correctly
  cluster_units <- which(membership == 1)
  selected_ids <- sort(unique(df$ids[df$sid %in% cluster_units]))

  expect_equal(selected_ids, sort(correct_flat_ids))
})
