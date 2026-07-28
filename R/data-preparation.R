#' Prepare data in long format
#'
#' Convert spatio-temporal data to long format with a flat spatial index (`ids`) and
#' ordered observation indices (`id_<dimname>` for each non-spatial dimension).
#' This is a pure converter: all rows are returned including cells that are all-NA.
#' No filtering is applied. Use `filter_df()` downstream to restrict to valid cells.
#'
#' @param x A `stars` object containing the spatial and functional dimensions.
#' @param spnames Character vector with the names of the spatial dimensions of `x`.
#'        Use a single name (e.g. `"geometry"`) for vector geometry data, or two names
#'        (e.g. `c("x", "y")`) for raster data. Functional dimensions are derived
#'        automatically as all remaining dimensions.
#'
#' @return A long-format data frame with columns `id` (flat array position in the original
#'         stars object, column-major over all dimensions), `ids` (flat spatial index,
#'         column-major over `spnames`, `spnames[1]` varies fastest), `id_<dimname>` for
#'         each functional dimension, and all variables from `x`. All rows are
#'         returned; `ids` is not remapped to `1..n_valid`.
#'
#' @examples
#'
#' library(sfclust)
#' library(stars)
#'
#' dims <- st_dimensions(
#'   geometry = st_sfc(lapply(1:5, function(i) st_point(c(i, i)))),
#'   time = seq(as.Date("2024-01-01"), by = "1 day", length.out = 3)
#' )
#' stdata <- st_as_stars(cases = array(1:15, dim = c(5, 3)), dimensions = dims)
#'
#' data_all(stdata)
#'
#' @importFrom stars expand_dimensions
#' @export
data_all <- function(x, spnames = "geometry") {
  validate_stdata_input(x, spnames)
  spnames <- spnames[order(match(spnames, dimnames(x)))]
  x[["id"]] <- seq_len(prod(dim(x)))

  # dimensions in long format
  dims0 <- expand_dimensions(x)
  spsizes <- dim(x)[spnames]
  dims <- list()
  for (d in names(dims0)) {
    if (d %in% spnames) dims[[d]] <- seq_len(spsizes[[d]])
    else dims[[paste0("id_", d)]] <- order(dims0[[d]])
  }
  dims <- do.call(expand.grid, dims)
  dims <- cbind(ids = spatial_index(dims, spsizes), dims[!names(dims) %in% spnames])

  # merge dimensions and dataframe
  df <- as.data.frame(x)
  cbind(df["id"], dims, df[, !names(df) %in% c("id", spnames), drop = FALSE])
}

# Verigy if spnames is part of the stars object
validate_stdata_input <- function(x, spnames) {
  if (!inherits(x, "stars")) {
    stop("Argument `x` must be a `stars` object.")
  }
  if (any(!(spnames %in% dimnames(x)))) {
    stop("Dimension names in `spnames` not found in stars object.")
  }
}

# Try to detect the spnames of the stars object respecting the order
detect_spnames <- function(x, spnames = NULL) {
  if (is.null(spnames)) {
    dims <- st_dimensions(x)
    spnames <- names(dims)[!is.na(sapply(dims, function(d) d$delta))]
    if (length(spnames) != 2) {
      geom_dims <- names(dims)[sapply(dims, function(d) inherits(d$values, "sfc"))]
      if (length(geom_dims) == 1) spnames <- geom_dims
      else stop("Could not auto-detect spatial dimensions from `stars` object. Provide `spnames`.")
    }
  }
  spnames[order(match(spnames, dimnames(x)))]
}

# Create unique spatial id
spatial_index <- function(df, sp_sizes) {
  if (length(sp_sizes) == 1L) {
    df[[names(sp_sizes)]]
  } else {
    strides <- cumprod(c(1L, unname(sp_sizes[-length(sp_sizes)])))
    as.integer(as.matrix(df[names(sp_sizes)] - 1L) %*% strides) + 1L
  }
}


#' Filter a long-format data frame to valid spatial cells
#'
#' Keeps only rows whose spatial index (`ids`) appears in `valid_ids` and adds
#' a `sid` column (1..n_valid) giving each valid cell's rank among the valid
#' cells. `sid` matches the indexing of `membership` vectors returned by
#' [genclust()]. If `valid_ids` is `NULL`, all rows are kept and `sid` is
#' assigned by rank of `ids`.
#'
#' @param df A long-format data frame as returned by [data_all()].
#' @param valid_ids Integer vector of valid flat spatial positions, as returned
#'   in `genclust(...)$valid_ids`. If `NULL`, all spatial cells are treated as valid.
#'
#' @return A filtered data frame with rows corresponding to valid spatial cells only,
#'   with an additional `sid` column (integer, 1..n_valid).
filter_df <- function(df, valid_ids = NULL) {
  if (is.null(valid_ids)) valid_ids <- sort(unique(df$ids))
  df_filtered <- df[df$ids %in% valid_ids, , drop = FALSE]
  df_filtered$sid <- match(df_filtered$ids, valid_ids)
  df_filtered
}
