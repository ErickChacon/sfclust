
#' Evaluate Log Marginal Likelihood for Each Cluster Using INLA
#'
#' This function computes the log marginal likelihood for a specific cluster `k` within a dataset.
#' It supports multiple family distributions and allows for the inclusion of exposure or trial sizes if applicable.
#' The function is flexible, accepting custom formulas and additional arguments for the INLA model fitting.
#'
#' @param k Integer, the cluster index for which to compute the log marginal likelihood.
#' @param stars_obj A stars object containing spatial data, including covariates and response variables.
#' @param membership Integer vector, indicating the cluster membership of each observation.
#' @param formula Object of class \code{\link[stats]{formula}}, specifying the model to be fitted.
#' @param family Character string specifying the family of the model, with support for "normal", "poisson",
#' "binomial", and "nbinomial" (negative binomial variant 1); defaults to "normal".
#' @param correction Logical, indicating whether a correction should be applied to the computed log marginal
#' likelihood; defaults to FALSE.
#' @param time_var The name of the time variable in the stars object.
#' @param N_var The name of the variable used as exposure or trials (for Poisson, binomial, etc.).
#' @param detailed Logical, specifying whether to return the full model object or just the log marginal
#' likelihood; defaults to FALSE.
#' @param ... Additional arguments passed to the INLA function.
#'
#' @return Depending on the `detailed` parameter:
#'   - If `detailed` is TRUE, returns the full `INLA` model object.
#'   - If `detailed` is FALSE, returns the log marginal likelihood (corrected if `correction` is TRUE).
#'
#' @examples
#' \dontrun{
#' stars_obj <- create_stars_object(...)
#' membership <- sample(1:2, 10, replace = TRUE)
#' formula <- Y ~ temperature + precipitation + f(ids, model = "iid")
#' result <- log_mlik_each(1, stars_obj, membership, formula = formula, family = "poisson", correction = TRUE)
#' }
#' @export
log_mlik_each <- function(k, stars_obj, membership, formula, family = "normal", correction = FALSE, time_var, N_var, detailed = FALSE, ...) {
  inla_data <- preprocess_data_each(stars_obj, k, membership, time_var = time_var)
  if(!is.null(N_var)){N <- get(N_var, envir = as.environment(stars_obj))}
  if (family == "poisson") {
    model <- INLA::inla(formula, family,
                        E = N,  # Use N from the stars object, if present
                        data = inla_data, control.predictor = list(compute = TRUE),
                        control.compute = list(config = TRUE), ...
    )
  } else if (family == "binomial") {
    model <- INLA::inla(formula, family,
                        Ntrials = N,  # Use N from the stars object, if present
                        data = inla_data, control.predictor = list(compute = TRUE),
                        control.compute = list(config = TRUE), ...
    )
  } else if (family == "nbinomial") {
    model <- INLA::inla(formula, family,
                        control.family = list(variant = 1), Ntrials = N,  # Use N from stars_obj
                        data = inla_data, control.predictor = list(compute = TRUE),
                        control.compute = list(config = TRUE), ...
    )
  } else if (family == "normal" | is.null(family)) {
    model <- INLA::inla(formula, family,
                        data = inla_data, control.predictor = list(compute = TRUE),
                        control.compute = list(config = TRUE), ...
    )
  }

  if (detailed) {
    return(model)
  } else {
    if (correction) {
      return(model[["mlik"]][[1]] + log_mlik_correction(model, formula))
    } else {
      return(model[["mlik"]][[1]])
    }
  }
}

#########################################################################################################
#' Calculate Log Marginal Likelihood for All Clusters
#'
#' This function computes the log marginal likelihood for each cluster specified in the membership vector.
#' It applies a specified model to each cluster to calculate the likelihood.
#'
#' @param stars_obj A stars object containing spatial data, including covariates and response variables.
#' @param membership Integer vector indicating the cluster membership for each observation.
#' @param formula An object of class \code{\link[stats]{formula}} specifying the model to be used in INLA.
#' @param family Character string specifying the family of distributions to use for the model. Defaults to "normal".
#' @param correction Logical indicating whether a correction for dispersion or other factors should be applied. Defaults to FALSE.
#' @param time_var The name of the time variable in the stars object.
#' @param N_var The name of the variable used as exposure or trials (for Poisson, binomial, etc.).
#' @param ... Additional arguments passed to the underlying `log_mlik_each` function.
#'
#' @return A numeric vector containing the log marginal likelihood for each cluster.
#'
#' @examples
#' \dontrun{
#' log_mlikelihoods <- log_mlik_all(stars_obj, membership, family = "poisson", formula = Y ~ temperature + precipitation, time_var = "time", N_var = "expected_cases")
#' }
#' @export
log_mlik_all <- function(stars_obj, membership, formula, family = "normal", correction = FALSE, time_var, N_var, ...) {
  k <- max(membership)
  sapply(1:k, log_mlik_each, stars_obj, membership, formula, family, correction, time_var, N_var, ...)
}

########################################################################################################################
#' This function preprocesses the data from a stars object by subsetting based on a cluster
#' and converting the data to a long format with proper indexing for time and spatial location.
#'
#' @param stars_obj A stars object containing spatial data.
#' @param k The cluster number to subset the data for.
#' @param membership A vector or matrix that defines cluster membership for each region.
#' @param time_var The name of the time variable (default is 'time').
#'
#' @return A long-format data frame with ids for spatial and time indexing, suitable for INLA.
#' @examples
#' preprocess_data_each(stars_obj, k = 1, membership, time_var = "time")
#'
#' @export
preprocess_data_each <- function(stars_obj, k, membership, time_var = "time") {
  # Check if the input is a stars object
  if (!inherits(stars_obj, "stars")) {
    stop("Input must be a 'stars' object.")
  }

  # Get the indices of the regions that belong to the kth cluster
  cluster_regions <- which(membership == k)
  dims <- st_dimensions(stars_obj)
  space_dim_pos <- which(names(dims) == time_var)
  if (is.vector(membership)) {

    # Subset the stars object to only include the regions in the kth cluster
    stars_obj_subset <- stars_obj[, ,cluster_regions, drop = FALSE]
  }

  # Convert the subsetted stars object to a long data frame
  data_long <- as.data.frame(stars_obj_subset, long = TRUE)

  # Calculate the number of time points (nt) and spatial locations (nk) in the cluster
  nt <- length(unique(data_long[[time_var]]))  # number of time points
  nk <- length(cluster_regions)  # number of regions in the kth cluster

  # Create id, idt, and ids ensuring correct indexing based on spatial locations and time
  data_long$id <- seq_len(nt * nk)  # Create unique IDs for each observation
  data_long$idt <- rep(seq_len(nt), time = nk)  # Create time index (idt) for each region
  data_long$ids <- rep(seq_len(nk), each = nt)  # Create spatial index (ids) for each time point

  # Remove the first column (if it contains geometry or spatial information)
  data_inla <- data_long[, -1]

  return(data_inla)  # Return the processed data for the kth cluster
}

##############################################################

## Auxiliary function to correct the marginal likelihood of INLA model
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

  # Modify the regular expression to handle cases with or without spaces around the assignment
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

log_mlik_corrected <- function(model, formula) {
  Slist <- get_structure_matrix(model, formula)
  Slogdet <- sapply(Slist, function(x) 2 * sum(log(Matrix::diag(SparseM::chol(x)))))
  model[["mlik"]][[1]] + 0.5 * sum(Slogdet)
}

log_mlik_correction <- function(model, formula) {
  Slist <- get_structure_matrix(model, formula)
  Slogdet <- sapply(Slist, function(x) 2 * sum(log(Matrix::diag(SparseM::chol(x)))))
  0.5 * sum(Slogdet)
}
