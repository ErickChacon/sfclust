#' Fit models and compute the log marginal likelihood for all clusters
#'
#' Fit the specified INLA model to each cluster and compute the log marginal likelihood
#' for each cluster specified in the membership vector.
#'
#' @param membership Integer, character or factor vector indicating the cluster membership
#'        for each spatial unit.
#' @param stdata A stars object with spatial-temporal dimensions defined in `stnames`, and
#'        including predictors and response variables.
#' @param stnames The names of the `spatial` and `temporal` dimensions of the `stdata`
#'        object.
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
#' log_mlik_all(c(1, 1, 1, 2, 2), stdata,
#'   formula = cases ~ temperature, family = "poisson", E = expected)
#'
#' models = log_mlik_all(c(1, 1, 1, 2, 2), stdata, detailed = TRUE,
#'   formula = cases ~ temperature, family = "poisson", E = expected)
#' lapply(models, summary)
#' }
#'
#' @export
log_mlik_all <- function(membership, stdata, stnames = c("geometry", "time"),
                         correction = TRUE, detailed = FALSE, ...) {
  clusters <- unique_clusters(membership)

  if (detailed) {
    lapply(clusters, log_mlik_each, membership, stdata, stnames, correction, detailed, ...)
  } else {
    sapply(clusters, log_mlik_each, membership, stdata, stnames, correction, detailed, ...)
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

log_mlik_each <- function(k, membership, stdata, stnames = c("geometry", "time"),
                          correction = TRUE, detailed = FALSE, ...) {
  inla_data <- data_each(k, membership, stdata, stnames)
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
#' Convert spatio-temporal data to long format with indices for time and spatial location.
#'
#' @param stdata A stars object containing spatial-temporal dimensions defined in `stnames`.
#' @param stnames The names of the `spatial` and `temporal` dimensions.
#'
#' @return A long-format data frame with ids for each observation and  for spatial and
#'         time indexing.
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
#' @export
data_all <- function(stdata, stnames = c("geometry", "time")) {
  # check if the input is a stars object and dimension names
  if (!inherits(stdata, "stars")) {
    stop("Argument `stdata` must be a `stars` object.")
  }
  if (any(!(stnames %in% dimnames(stdata)))) {
    stop("Provided dimension names in `stnames` not found in stars object `stdata`.")
  }

  stdata[["id"]] <- 1:prod(dim(stdata))

  # spatio-temporal dimension in long format
  dims <- expand_dimensions(stdata)
  dims[[stnames[1]]] <- seq_along(dims[[stnames[1]]])
  dims[[stnames[2]]] <- order(dims[[stnames[2]]])
  dims <- setNames(do.call(expand.grid, dims)[stnames], c("ids", "idt"))

  # merge dimensions and dataframe
  stdata <- as.data.frame(stdata)
  cbind(stdata["id"], dims, stdata[, !names(stdata) %in% c("id", stnames[1])])
}

#' Prepare data for a cluster
#'
#' Subset a spatio-temporal dataset for a cluster and convert it to a long format with
#' indices for time and spatial location.
#'
#' @param k The cluster number to subset.
#' @param membership A vector defining the cluster membership for each region.
#' @param stdata A stars object containing spatial-temporal dimensions defined in `stnames`.
#' @param stnames The names of the `spatial` and `temporal` dimensions.
#'
#' @return A long-format data frame with ids for each observation and  for spatial and
#'         time indexing.
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
#' @import cubelyr stars
#' @importFrom dplyr filter
#'
#' @export
data_each <- function(k, membership, stdata, stnames = c("geometry", "time")) {
  # check if the input is a stars object and dimension names
  if (!inherits(stdata, "stars")) {
    stop("Argument `stdata` must be a `stars` object.")
  }
  if (any(!(stnames %in% dimnames(stdata)))) {
    stop("Provided dimension names in `stnames` not found in stars object `stdata`.")
  }

  stdata[["id"]] <- 1:prod(dim(stdata))

  # subset cluster k
  cluster_elements <- which(membership == k)
  stdata <- filter(stdata, !!as.name(stnames[1]) %in% cluster_elements)

  # spatio-temporal dimension in long format
  dims <- expand_dimensions(stdata)
  dims[[stnames[1]]] <- cluster_elements
  dims[[stnames[2]]] <- order(dims[[stnames[2]]])
  dims <- setNames(do.call(expand.grid, dims)[stnames], c("ids", "idt"))

  # merge dimensions and dataframe
  stdata <- as.data.frame(stdata)
  cbind(stdata["id"], dims, stdata[, !names(stdata) %in% c("id", stnames[1])])
}
