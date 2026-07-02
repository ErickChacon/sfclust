#' Fit models and compute the log marginal likelihood for all clusters
#'
#' Fit the specified INLA model to each cluster and compute the log marginal likelihood
#' for each cluster specified in the membership vector.
#'
#' @param membership Integer, character or factor vector indicating the cluster membership
#'        for each spatial unit.
#' @param data A long-format data frame as returned by [filter_df()], with columns `id`
#'        (flat array position), `ids` (flat spatial position), `sid` (valid-cell rank,
#'        1 to n_valid, matching `membership` indices), observation index columns
#'        (`id_<dimname>`), and all response/covariate variables.
#' @param correction Logical value indicating whether a correction for dispersion.
#' @param detailed Logical value indicating whether to return the INLA model instead of
#'        the log marginal likelihood. The argument `correction` is not applied in this
#'        case.
#' @param ... Arguments passed to the `inla` function (eg. `family`, `formula` and `E`).
#'
#' @return A numeric vector containing the log marginal likelihood for each cluster or the
#'         the fitted INLA model for each cluster when `detailed = TRUE`.
#'
#' @examples
#'
#'
#' \donttest{
#' library(sfclust)
#' library(stars)
#'
#' dims <- st_dimensions(
#'   geometry = st_sfc(lapply(1:5, function(i) st_point(c(i, i)))),
#'   time = seq(as.Date("2024-01-01"), by = "1 day", length.out = 3)
#' )
#' stdata <- st_as_stars(
#'   cases = array(rpois(15, 100 * exp(1)), dim = c(5, 3)),
#'   temperature = array(runif(15, 15, 20), dim = c(5, 3)),
#'   expected = array(100, dim = c(5, 3)),
#'   dimensions = dims
#' )
#'
#' df <- data_all(stdata)
#' log_mlik_all(c(1, 1, 1, 2, 2), df,
#'   formula = cases ~ temperature, family = "poisson", E = expected)
#'
#' models <- log_mlik_all(c(1, 1, 1, 2, 2), df, detailed = TRUE,
#'   formula = cases ~ temperature, family = "poisson", E = expected)
#' lapply(models, summary)
#' }
#'
#' @export
log_mlik_all <- function(membership, data, correction = TRUE, detailed = FALSE, ...) {
  clusters <- unique_clusters(membership)

  if (detailed) {
    lapply(clusters, log_mlik_each, membership, data, correction, detailed, ...)
  } else {
    sapply(clusters, log_mlik_each, membership, data, correction, detailed, ...)
  }
}

get_data <- function(object) {
  attr(object, "data")
}

unique_clusters <- function (membership) {
  if (is.character(membership)) membership <- as.factor(membership)
  if (is.factor(membership)){
    clusters <- levels(membership)
    setNames(clusters, clusters)
  } else if (is.numeric(membership)) {
    if (all(1:max(membership) %in% membership)) {
      1:max(membership)
    } else {
      stop("`membership` vector does not contain all groups until `max(membership)`.")
    }
  } else {
    stop("`membership` vector does not contain all groups until `max(membership)`.")
  }
}

log_mlik_each <- function(k, membership, data, correction = TRUE, detailed = FALSE, ...) {
  cluster_units <- which(membership == k)
  inla_data <- data[data$sid %in% cluster_units, , drop = FALSE]
  model <- INLA::inla(
    data = inla_data,
    control.predictor = list(compute = TRUE),
    control.compute = list(config = correction),
    ...
  )

  if (detailed) {
    model
  } else if (!correction) {
    model[["mlik"]][[1]]
  } else {
    model[["mlik"]][[1]] + log_mlik_correction(model)
  }
}

#' @importFrom Matrix diag
log_mlik_correction <- function(model) {
  Slist <- get_structure_matrix(model)
  if (!length(Slist) == 0) {
    Slogdet <- sapply(Slist, function(x) 2 * sum(log(Matrix::diag(SparseM::chol(x)))))
    0.5 * sum(Slogdet)
  } else {
    # warning("No structure matrix found to apply correction.")
    0.0
  }
}

get_structure_matrix <- function(model) {
  # obtain settings from model
  prior_diagonal <- model[[".args"]][["control.compute"]][["control.gcpo"]][["prior.diagonal"]] # 0.0001
  formula <- model[[".args"]][["formula"]]
  model <- model[["misc"]][["configs"]]

  # effects dimension information
  x_info <- model[["contents"]]
  ef_start <- setNames(x_info$start[-1] - x_info$length[1], x_info$tag[-1])
  ef_end <- ef_start + x_info$length[-1] - 1

  # select effect that requires correction
  effs_to_correct <- correction_required(formula)

  # provide structure matrix for selected effects
  ind <- which.max(sapply(model[["config"]], function(x) x$log.posterior))

  out <- list()
  for (x in effs_to_correct) {
    i <- ef_start[x]
    j <- ef_end[x]
    Qaux <- model[["config"]][[ind]][["Qprior"]][i:j, i:j]
    Matrix::diag(Qaux) <- Matrix::diag(Qaux) - prior_diagonal
    Qaux <- Qaux /
      exp(model[["config"]][[ind]][["theta"]][paste0("Log precision for ", x)])
    Matrix::diag(Qaux) <- Matrix::diag(Qaux) + prior_diagonal
    out[[x]] <- Qaux
  }

  return(out)
}

#' @importFrom stats terms
correction_required <- function (formula) {
  effects <- as.list(attr(terms(formula), "variables"))[c(-1, -2)]
  need_correction <- grepl("model\\s*=\\s*\"c{0,1}rw", sapply(effects, deparse1))
  sapply(effects, all.vars)[need_correction]
}

# Core spatial index helper: given a long-format data frame of dimension indices
# and named sp_sizes, returns
# the flat column-major spatial index (spnames[1] varies fastest). Used by
# data_all() and create_domain().
spatial_index <- function(df, sp_sizes) {
  if (length(sp_sizes) == 1L) {
    df[[names(sp_sizes)]]
  } else {
    strides <- cumprod(c(1L, unname(sp_sizes[-length(sp_sizes)])))
    as.integer(as.matrix(df[names(sp_sizes)] - 1L) %*% strides) + 1L
  }
}

#' Prepare data in long format
#'
#' Convert spatio-temporal data to long format with a flat spatial index (`ids`) and
#' ordered observation indices (`idf_<dimname>` for each non-spatial dimension).
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

validate_stdata_input <- function(x, spnames) {
  if (!inherits(x, "stars")) {
    stop("Argument `x` must be a `stars` object.")
  }
  if (any(!(spnames %in% dimnames(x)))) {
    stop("Dimension names in `spnames` not found in stars object.")
  }
}


#' @importFrom stars st_apply
create_domain <- function(x, spnames, response = NULL) {
  if (is.null(response)) {
    st_apply(x[1], spnames, function(v) TRUE)
  } else {
    st_apply(x[response], spnames, function(v) any(!is.na(v)))
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
