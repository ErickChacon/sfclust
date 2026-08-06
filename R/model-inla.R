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
#' @param inla_args A named list or pairlist of arguments passed to `inla()` (e.g.
#'        `formula`, `family`, `E`).
#' @param within_model Optional object created by `sfclust_within_model()`. When
#'        provided with a `stack_fun`, the function is called once per cluster to
#'        build the INLA stack used for that cluster.
#'
#' @details
#' When a within-model object has a `stack_fun`, it is called as `stack_fun(data)`.
#' The `data` argument is the long-format data for one cluster. It contains `id`,
#' `ids`, `sid`, functional index columns such as `id_time`, and all variables
#' from the original data.
#'
#' `stack_fun(data)` may return an `INLA::inla.stack()` object directly, or a
#' list with components `stack` and optional `formula`. Family and other INLA
#' arguments should be provided through `sfclust_within_model()` or `sfclust()`,
#' not returned by `stack_fun`.
#'
#' @return A numeric vector containing the log marginal likelihood for each cluster or the
#'         the fitted INLA model for each cluster when `detailed = TRUE`.
#'
log_mlik_all <- function(membership, data, correction = TRUE, detailed = FALSE, inla_args = NULL, within_model = NULL) {
  clusters <- unique_clusters(membership)

  if (detailed) {
    lapply(clusters, log_mlik_each, membership, data, correction, detailed, inla_args, within_model)
  } else {
    sapply(clusters, log_mlik_each, membership, data, correction, detailed, inla_args, within_model)
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

log_mlik_each <- function(k, membership, data, correction = TRUE, detailed = FALSE, inla_args = NULL, within_model = NULL) {
  cluster_units <- which(membership == k)
  inla_data <- data[data$sid %in% cluster_units, , drop = FALSE]
  inla_args <- as.list(inla_args)

  inla_call <- list(
    data = inla_data,
    control.predictor = list(compute = TRUE),
    control.compute = list(config = correction)
  )

  if (!is.null(within_model$stack_fun)) {
    inla_stack <- within_model$stack_fun(inla_data)
    stack_formula <- NULL
    if (!inherits(inla_stack, "inla.data.stack")) {
      if (!is.list(inla_stack)) {
        stop("`stack_fun` must return either an INLA stack or a list with `stack` and optional `formula` components.", call. = FALSE)
      }
      if (is.null(names(inla_stack)) || !"stack" %in% names(inla_stack)) {
        stop("A list returned by `stack_fun` must contain a `stack` component.", call. = FALSE)
      }
      extra_components <- setdiff(names(inla_stack), c("stack", "formula"))
      if (length(extra_components) > 0) {
        stop("A list returned by `stack_fun` may only contain `stack` and optional `formula` components.", call. = FALSE)
      }
      stack_formula <- inla_stack$formula
      inla_stack <- inla_stack$stack
    }
    if (!inherits(inla_stack, "inla.data.stack")) {
      stop("The `stack` component returned by `stack_fun` must be an INLA stack.", call. = FALSE)
    }
    if (!is.null(stack_formula) && !inherits(stack_formula, "formula")) {
      stop("The `formula` component returned by `stack_fun` must be a formula.", call. = FALSE)
    }
    if (!is.null(stack_formula)) inla_args$formula <- stack_formula

    inla_call$data <- INLA::inla.stack.data(inla_stack)
    inla_call$control.predictor$A <- INLA::inla.stack.A(inla_stack)
  }

  if (!is.null(inla_args$control.predictor)) {
    user_control_predictor <- inla_args$control.predictor
    user_control_predictor$compute <- TRUE
    if (!is.null(inla_call$control.predictor$A)) {
      user_control_predictor$A <- inla_call$control.predictor$A
    }
    inla_call$control.predictor <- user_control_predictor
  }

  if (!is.null(inla_args$control.compute)) {
    user_control_compute <- inla_args$control.compute
    user_control_compute$config <- correction
    inla_call$control.compute <- user_control_compute
  }

  inla_args <- inla_args[setdiff(names(inla_args), names(inla_call))]
  inla_call <- c(inla_call, inla_args)

  model <- tryCatch(
    do.call(INLA::inla, inla_call),
    error = function(e) {
      warning("INLA failed for cluster ", k, ": ", conditionMessage(e), immediate. = TRUE)
      NULL
    }
  )

  if (is.null(model)) return(if (detailed) NULL else -Inf)

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
