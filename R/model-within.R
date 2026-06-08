#' Fit models and compute the log marginal likelihood for all clusters
#'
#' Fit the specified INLA model to each cluster and compute the log marginal likelihood
#' for each cluster specified in the membership vector.
#'
#' @param membership Integer, character or factor vector indicating the cluster membership
#'        for each spatial unit.
#' @param data A long-format data frame as returned by [data_all()], with columns `id`
#'        (flat index), `ids` (spatial unit index, 1 to ns), functional index columns
#'        (`idf_<dimname>`), and all response/covariate variables.
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
  if (inherits(object, "sfclust_stars")) {
    data_all(attr(object, "stdata"), attr(object, "sp_dims"), attr(object, "fun_dims"))
  } else {
    attr(object, "args")$data
  }
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
  inla_data <- data[data$ids %in% which(membership == k), , drop = FALSE]
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

#' Prepare data in long format
#'
#' Convert spatio-functional data to long format with a flat spatial index (`ids`) and
#' ordered functional indices (`idf_<dimname>` for each functional dimension).
#'
#' @param stdata A stars object containing the spatial and functional dimensions.
#' @param sp_dims Character vector with the names of the spatial dimensions of `stdata`.
#'        Use a single name (e.g. `"geometry"`) for vector geometry data, or two names
#'        (e.g. `c("x", "y")`) for raster data.
#' @param fun_dims Character vector with the names of the functional dimensions of `stdata`
#'        (e.g. `"time"`).
#'
#' @return A long-format data frame with columns `id` (flat array index), `ids` (flat
#'         spatial index, 1 to ns), `idf_<dimname>` for each functional dimension, and
#'         all variables from `stdata`.
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
data_all <- function(stdata, sp_dims = "geometry", fun_dims = "time") {
  validate_stdata_input(stdata, sp_dims, fun_dims)
  stdata[["id"]] <- 1:prod(dim(stdata))
  valid_ids <- attr(stdata, "valid_ids")

  dims <- expand_dimensions(stdata)

  # replace spatial dim values with sequential indices, functional dims with order
  sp_sizes <- setNames(sapply(sp_dims, function(d) length(dims[[d]])), sp_dims)
  for (d in sp_dims) dims[[d]] <- seq_len(sp_sizes[[d]])
  for (d in fun_dims) dims[[d]] <- order(dims[[d]])

  # full grid; row order matches as.data.frame(stdata)
  grid <- do.call(expand.grid, dims)

  # flat spatial index (column-major over sp_dims)
  if (length(sp_dims) == 1) {
    ids <- grid[[sp_dims]]
  } else {
    strides <- cumprod(c(1L, as.integer(sp_sizes[-length(sp_sizes)])))
    ids <- as.integer(as.matrix(grid[sp_dims] - 1L) %*% strides) + 1L
  }

  fun_cols <- setNames(as.data.frame(grid[fun_dims]), paste0("idf_", fun_dims))

  stdata_df <- as.data.frame(stdata)
  result <- cbind(stdata_df["id"], ids = ids, fun_cols,
                  stdata_df[, !names(stdata_df) %in% c("id", sp_dims), drop = FALSE])

  # for raster: filter to valid (non-all-NA) pixels and remap ids to analysis indices
  if (!is.null(valid_ids)) {
    result <- result[result$ids %in% valid_ids, ]
    result$ids <- match(result$ids, valid_ids)
    result$id  <- seq_len(nrow(result))
  }

  result
}

#' Prepare data for a cluster
#'
#' Subset a spatio-functional dataset to a single cluster and convert it to long format.
#'
#' @param k The cluster number to subset.
#' @param membership A vector defining the cluster membership for each spatial unit.
#' @param stdata A stars object containing the spatial and functional dimensions.
#' @param sp_dims Character vector with the names of the spatial dimensions of `stdata`.
#'        Use a single name (e.g. `"geometry"`) for vector geometry data, or two names
#'        (e.g. `c("x", "y")`) for raster data.
#' @param fun_dims Character vector with the names of the functional dimensions of `stdata`
#'        (e.g. `"time"`).
#'
#' @return A long-format data frame for cluster `k` with columns `id`, `ids`, `idf_<dimname>`
#'         for each functional dimension, and all variables from `stdata`.
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
#' data_each(k = 2, membership = c(1, 1, 1, 2, 2), stdata)
#'
#' @importFrom stars expand_dimensions
#' @importFrom dplyr filter
#' @export
data_each <- function(k, membership, stdata, sp_dims = "geometry", fun_dims = "time") {
  validate_stdata_input(stdata, sp_dims, fun_dims)
  stdata[["id"]] <- 1:prod(dim(stdata))

  cluster_sp_ids <- which(membership == k)

  if (length(sp_dims) == 1) {
    stdata <- filter(stdata, !!as.name(sp_dims) %in% cluster_sp_ids)

    dims <- expand_dimensions(stdata)
    dims[[sp_dims]] <- cluster_sp_ids
    for (d in fun_dims) dims[[d]] <- order(dims[[d]])

    grid <- do.call(expand.grid, dims)
    fun_cols <- setNames(as.data.frame(grid[fun_dims]), paste0("idf_", fun_dims))

    stdata_df <- as.data.frame(stdata)
    cbind(stdata_df["id"], ids = grid[[sp_dims]], fun_cols,
          stdata_df[, !names(stdata_df) %in% c("id", sp_dims), drop = FALSE])
  } else {
    all_data <- data_all(stdata, sp_dims, fun_dims)
    all_data[all_data$ids %in% cluster_sp_ids, , drop = FALSE]
  }
}

validate_stdata_input <- function(stdata, sp_dims, fun_dims) {
  if (!inherits(stdata, "stars")) {
    stop("Argument `stdata` must be a `stars` object.")
  }
  all_dims <- c(sp_dims, fun_dims)
  if (any(!(all_dims %in% dimnames(stdata)))) {
    stop("Dimension names in `sp_dims` and `fun_dims` not found in stars object `stdata`.")
  }
}
