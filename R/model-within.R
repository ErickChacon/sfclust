#' Create a within-cluster model specification
#'
#' Defines the INLA model fitted inside each cluster. The usual formula and
#' family interface remains available directly in `sfclust()`. This object is
#' useful when the within-cluster model needs a data preparation step,
#' such as an INLA stack or a custom function to obtain fitted values.
#'
#' @param formula An INLA formula. It may be `NULL` only when `stack_fun`
#'   returns a cluster-specific formula.
#' @param family The INLA family specification.
#' @param stack_fun Optional function called once per cluster. It receives the
#'   long-format cluster data and returns either an `INLA::inla.stack()` object
#'   or a list with components `stack` and optional `formula`.
#' @param fitted_fun Optional function used by `fitted()` to extract fitted
#'   values from custom within-cluster models.
#' @details
#' `stack_fun` is called as `stack_fun(data)`, where `data` contains the rows
#' for one cluster. The data are in the same long format used internally by
#' `sfclust()`: `id` is the row position in the original data, `ids` is the
#' spatial unit id, `sid` is the valid spatial unit id used by the clustering
#' algorithm, and columns such as `id_time` identify functional dimensions.
#'
#' `stack_fun(data)` may return an INLA stack directly:
#'
#' \preformatted{
#' INLA::inla.stack(...)
#' }
#'
#' or a list:
#'
#' \preformatted{
#' list(
#'   stack = INLA::inla.stack(...),
#'   formula = optional_formula
#' )
#' }
#'
#' The `formula` component is optional. Family and other INLA arguments should
#' be provided through `sfclust_within_model()` or `sfclust()`, not returned by
#' `stack_fun`.
#'
#' @return An object of class `"sfclust_within_model"`.
#'
#' @examples
#' model <- sfclust_within_model(
#'   formula = y ~ 1,
#'   family = "gaussian"
#' )
#'
#' @export
sfclust_within_model <- function(formula = NULL, family,
                                 stack_fun = NULL, fitted_fun = NULL) {
  if (missing(family)) family <- NULL

  object <- list(
    formula = formula,
    family = family,
    stack_fun = stack_fun,
    fitted_fun = fitted_fun
  )
  class(object) <- "sfclust_within_model"

  validate_sfclust_within_model(object)
  object
}

validate_sfclust_within_model <- function(object) {
  if (!inherits(object, "sfclust_within_model")) {
    stop("`object` must inherit from class \"sfclust_within_model\".", call. = FALSE)
  }

  if (is.null(object$formula) && is.null(object$stack_fun)) {
    stop("`formula` may be `NULL` only when `stack_fun` is provided.", call. = FALSE)
  }

  if (!is.null(object$formula) && !inherits(object$formula, "formula")) {
    stop("`formula` must be a formula or `NULL`.", call. = FALSE)
  }

  if (is.null(object$family)) {
    stop("`family` is required.", call. = FALSE)
  }

  if (!is.null(object$stack_fun) && !is.function(object$stack_fun)) {
    stop("`stack_fun`, if provided, must be a function.", call. = FALSE)
  }

  if (!is.null(object$fitted_fun) && !is.function(object$fitted_fun)) {
    stop("`fitted_fun`, if provided, must be a function.", call. = FALSE)
  }

  invisible(object)
}

#' @export
print.sfclust_within_model <- function(x, ...) {
  validate_sfclust_within_model(x)

  cat("sfclust within-cluster model\n")

  cat("\nFormula:\n")
  if (is.null(x$formula)) {
    cat("<provided by stack_fun>\n")
  } else {
    print(x$formula, showEnv = FALSE, ...)
  }

  cat("\nFamily:\n")
  print(x$family, ...)

  cat("\nStack function:", if (!is.null(x$stack_fun)) "<function>" else "<none>", "\n")
  cat("Fitted function:", if (!is.null(x$fitted_fun)) "<function>" else "<none>", "\n")

  invisible(x)
}
