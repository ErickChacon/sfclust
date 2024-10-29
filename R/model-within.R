#' Fit models and compute the log marginal likelihood for all clusters
#'
#' Fit the specified INLA model to each cluster and compute the log marginal likelihood
#' for each cluster specified in the membership vector.
#'
#' @param membership Integer vector indicating the cluster membership for each spatial
#'        unit.
#' @param stdata A stars object with spatial-temporal dimensions defined in `stnames`, and
#'        including predictors and response variables.
#' @param stnames The names of the `spatial` and `temporal` dimensions of the `stdata`
#'        object.
#' @param correction Logical value indicating whether a correction for dispersion.
#' @param detailed Logical value indicating whether to return the INLA model instead of
#'        the log marginal likelihood.
#' @param ... Arguments passed to the `inla` function (eg. `family`, `formula` and `E`).
#'
#' @return A numeric vector containing the log marginal likelihood for each cluster or the
#'         the fitted INLA model for each cluster when `detailed = TRUE`.
#'
#' @examples
#' \dontrun{
#'
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
#'
#' }
#' @export
log_mlik_all <- function(membership, stdata, stnames = c("geometry", "time"),
                         correction = FALSE, detailed = FALSE, ...) {
  q <- max(membership)
  if (detailed) {
    lapply(1:q, log_mlik_each, membership, stdata, stnames, correction, detailed, ...)
  } else {
    sapply(1:q, log_mlik_each, membership, stdata, stnames, correction, detailed, ...)
  }
}

log_mlik_each <- function(k, membership, stdata, stnames = c("geometry", "time"),
                          correction = FALSE, detailed = FALSE, ...) {
  inla_data <- data_each(k, membership, stdata, stnames)
  model <- INLA::inla(
    data = inla_data,
    control.predictor = list(compute = TRUE),
    control.compute = list(config = TRUE),
    ...
  )

  if (detailed) {
    model
  } else if (!correction) {
    model[["mlik"]][[1]]
  } else {
    model[["mlik"]][[1]] + log_mlik_correction(model, formula)
  }
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
#' @return A long-format data frame with ids for spatial and time indexing.
#'
#' @examples
#'
#' library(stars)
#'
#' dims <- st_dimensions(
#'   geometry = st_sfc(lapply(1:5, function(i) st_point(c(i, i)))),
#'   time = seq(as.Date("2024-01-01"), by = "1 day", length.out = 3)
#' )
#' stdata <- st_as_stars(cases = array(1:15, dim = c(5, 3)), dimensions = dims)
#'
#' data_each(stdata, k = 2, membership = c(1, 1, 1, 2, 2))
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

  # subset cluster k
  cluster <- which(membership == k)
  stdata <- filter(stdata, !!as.name(stnames[1]) %in% cluster)

  # spatio-temporal dimenstion in long format
  dims <- expand_dimensions(stdata)
  dims[[stnames[1]]] <- seq_along(dims[[stnames[1]]])
  dims[[stnames[2]]] <- order(dims[[stnames[2]]])
  dims <- setNames(do.call(expand.grid, dims)[stnames], c("ids", "idt"))

  # merge dimensions and dataframe
  stdata <- as.data.frame(stdata)
  stdata[[stnames[1]]] <- NULL
  cbind(list(id = 1:nrow(dims)), dims, stdata)
}

log_mlik_correction <- function(model, formula) {
  Slist <- get_structure_matrix(model, formula)
  Slogdet <- sapply(Slist, function(x) 2 * sum(log(Matrix::diag(SparseM::chol(x)))))
  0.5 * sum(Slogdet)
}

get_structure_matrix <- function(model, formula) {
  # prior diagonal
  prior_diagonal <- model[[".args"]][["control.compute"]][["control.gcpo"]][["prior.diagonal"]]

  # model config
  model <- model[["misc"]][["configs"]]

  # effects dimension information
  x_info <- model[["contents"]]
  ef_start <- setNames(x_info$start[-1] - x_info$length[1], x_info$tag[-1])
  ef_end <- ef_start + x_info$length[-1] - 1

  # select effect that requires correction
  fs <- as.list(attr(terms(formula), "variables"))[c(-1, -2)]

  # modify the regular expression to handle cases with or without spaces around the assignment
  fs_rw <- grepl("model\\s*=\\s*\"rw", sapply(fs, deparse))

  fs_vars <- sapply(fs, all.vars)[fs_rw]

  # provide structure matrix for selected effects
  ind <- which.max(sapply(model[["config"]], function(x) x$log.posterior))

  out <- list()
  for (x in fs_vars) {
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
