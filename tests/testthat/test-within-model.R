library(Matrix)

# --- sfclust_within_model ----------------------------------------------

test_that("sfclust_within_model: validates model objects", {
  stack_fun <- function(data){data}
  fitted_fun <- function(model, data, cluster, ...){data.frame(cluster = cluster)}

  model <- sfclust_within_model(
    formula = y ~ 1,
    family = "gaussian",
    stack_fun = stack_fun,
    fitted_fun = fitted_fun
  )

  expect_s3_class(model, "sfclust_within_model")
  expect_identical(model$formula, y ~ 1)
  expect_identical(model$family, "gaussian")
  expect_identical(model$stack_fun, stack_fun)
  expect_identical(model$fitted_fun, fitted_fun)

  expect_s3_class(
    sfclust_within_model(formula = NULL, family = "gaussian", stack_fun = stack_fun),
    "sfclust_within_model"
  )

  expect_error(
    sfclust_within_model(family = "gaussian"),
    "`formula` may be `NULL` only when `stack_fun` is provided"
  )

  expect_error(
    sfclust_within_model(formula = "y ~ 1", family = "gaussian"),
    "`formula` must be a formula or `NULL`."
  )

  expect_error(
    sfclust_within_model(formula = y ~ 1),
    "`family` is required."
  )

  expect_error(
    sfclust_within_model(formula = y ~ 1, family = "gaussian", stack_fun = "not a function"),
    "`stack_fun`, if provided, must be a function."
  )

  expect_error(
    sfclust_within_model(formula = y ~ 1, family = "gaussian", fitted_fun = "not a function"),
    "`fitted_fun`, if provided, must be a function."
  )

  expect_output(print(model), "sfclust within-cluster model")
})

test_that("stack_fun list return is limited to stack and formula", {
  stack_fun <- function(data) {
    list(
      stack = structure(list(), class = "inla.data.stack"),
      formula = y ~ 1,
      family = "gaussian"
    )
  }
  model <- sfclust_within_model(formula = y ~ 1, family = "gaussian", stack_fun = stack_fun)
  data <- data.frame(y = 1:2, id = 1:2, ids = 1:2, sid = 1:2, id_time = 1:2)

  expect_error(
    sfclust:::log_mlik_each(1, c(1L, 1L), data, correction = FALSE, detailed = FALSE,
                  inla_args = list(formula = y ~ 1, family = "gaussian"),
                  within_model = model),
    "may only contain `stack` and optional `formula` components"
  )
})

# --- sfclust within_model params ------------------------------------

toy_data <- function() {
  ns <- 2
  nt <- 3
  data <- data.frame(
    id      = seq_len(ns * nt),
    ids     = rep(seq_len(ns), nt),
    id_time = rep(seq_len(nt), each = ns),
    y       = c(-0.3, 0.2, -0.1, 0.4, 0.1, 0.5)
  )
  adjacency <- sparseMatrix(
    i = 1, j = 2, x = 1,
    dims = c(ns, ns), symmetric = TRUE
  )

  list(data = data, adjacency = adjacency)
}

test_that("sfclust: includes within_model", {
  expect_true("within_model" %in% names(formals(sfclust:::sfclust.data.frame)))
  expect_true("within_model" %in% names(formals(sfclust:::sfclust.stars)))
})

test_that("sfclust: original formula and family interface still works", {
  skip_if_not_installed("INLA")
  toy_data_list <- toy_data()

  fit <- suppressWarnings(suppressMessages(
    sfclust(
      toy_data_list$data,
      adjacency = toy_data_list$adjacency,
      fnames = "id_time",
      nclust = 1,
      formula = y ~ 1,
      family = "gaussian",
      correction = FALSE,
      niter = 1,
      move_prob = c(0, 0, 0, 1)
    )
  ))

  expect_s3_class(fit, "sfclust")
  expect_null(attr(fit, "fit_args")$within_model)
  expect_equal(eval(attr(fit, "inla_args")$family), "gaussian")
})

test_that("sfclust: within_model interface works", {
  skip_if_not_installed("INLA")
  toy_data_list <- toy_data()
  model <- sfclust_within_model(
    formula = y ~ 1,
    family = "gaussian"
  )

  fit <- suppressWarnings(suppressMessages(
    sfclust(
      toy_data_list$data,
      adjacency = toy_data_list$adjacency,
      fnames = "id_time",
      nclust = 1,
      within_model = model,
      correction = FALSE,
      niter = 1,
      move_prob = c(0, 0, 0, 1)
    )
  ))

  expect_s3_class(fit, "sfclust")
  expect_s3_class(attr(fit, "fit_args")$within_model, "sfclust_within_model")
  expect_equal(attr(fit, "fit_args")$within_model$family, "gaussian")
})

test_that("sfclust: within_model and original model arguments cannot be mixed", {
  toy_data_list <- toy_data()
  model <- sfclust_within_model(formula = y ~ 1, family = "gaussian")

  expect_error(
    sfclust(
      toy_data_list$data,
      adjacency = toy_data_list$adjacency,
      within_model = model,
      formula = y ~ 1,
      niter = 1
    ),
    "Do not pass `within_model` together with direct model arguments"
  )
})

test_that("sfclust: within_model with stack_fun can fit and extract values", {
  skip_if_not_installed("INLA")
  toy_data_list <- toy_data()

  stack_fun <- function(data) {
    list(
      stack = INLA::inla.stack(
        data = list(y = data$y),
        A = list(1),
        effects = list(data.frame(intercept = 1, id_time = data$id_time)),
        tag = "cluster"
      ),
      formula = y ~ 0 + intercept
    )
  }

  fitted_fun <- function(model, data, cluster, object, aggregate, ...) {
    data.frame(
      cluster = cluster,
      id_time = data$id_time,
      mean = data$y,
      response = "y"
    )
  }

  model <- sfclust_within_model(
    formula = NULL,
    family = "gaussian",
    stack_fun = stack_fun,
    fitted_fun = fitted_fun
  )

  fit <- suppressWarnings(suppressMessages(
    sfclust(
      toy_data_list$data,
      adjacency = toy_data_list$adjacency,
      fnames = "id_time",
      nclust = 1,
      within_model = model,
      correction = FALSE,
      niter = 1,
      move_prob = c(0, 0, 0, 1)
    )
  ))

  expect_s3_class(fit, "sfclust")
  expect_s3_class(attr(fit, "inla_args")$formula, "formula")

  fit_df <- fitted(fit)
  expect_true(is.data.frame(fit_df))
  expect_true(all(c("cluster", "id_time", "mean", "response") %in% names(fit_df)))
  expect_s3_class(plot_clusters_fitted(fit), "ggplot")
})

test_that("plot_clusters_fitted: custom fitted output must be compatible", {
  toy_data_list <- toy_data()
  model <- sfclust_within_model(
    formula = y ~ 1,
    family = "gaussian",
    stack_fun = function(data) data,
    fitted_fun = function(model, data, cluster, object, aggregate, ...) {
      data.frame(cluster = cluster, value = seq_len(nrow(data)))
    }
  )

  fit <- list(
    samples = list(
      membership = matrix(c(1L, 1L), nrow = 1),
      log_mlike = 0,
      move_counts = c(births = 0, deaths = 0, changes = 0, hypers = 0)
    ),
    clust = list(
      id = 1,
      membership = c(1L, 1L),
      models = list(list())
    )
  )
  attr(fit, "fit_args") <- list(
    data = transform(toy_data_list$data, sid = ids),
    within_model = model
  )
  attr(fit, "input_args") <- list(fnames = "id_time")
  attr(fit, "inla_args") <- list(formula = y ~ 1, family = "gaussian")
  class(fit) <- "sfclust"

  expect_error(
    plot_clusters_fitted(fit),
    "requires `fitted_fun` to return columns"
  )
})

test_that("stack_fun validates returned objects", {
  toy_data_list <- toy_data()
  data <- transform(toy_data_list$data, sid = ids)
  membership <- c(1L, 1L)
  inla_args <- list(formula = y ~ 1, family = "gaussian")

  model <- sfclust_within_model(
    formula = y ~ 1,
    family = "gaussian",
    stack_fun = function(data) list(formula = y ~ 1)
  )
  expect_error(
    sfclust:::log_mlik_each(1, membership, data, FALSE, FALSE, inla_args, model),
    "must contain a `stack` component"
  )

  model <- sfclust_within_model(
    formula = y ~ 1,
    family = "gaussian",
    stack_fun = function(data) list(stack = "not a stack")
  )
  expect_error(
    sfclust:::log_mlik_each(1, membership, data, FALSE, FALSE, inla_args, model),
    "must be an INLA stack"
  )

  model <- sfclust_within_model(
    formula = y ~ 1,
    family = "gaussian",
    stack_fun = function(data) {
      list(
        stack = structure(list(), class = "inla.data.stack"),
        formula = "y ~ 1"
      )
    }
  )
  expect_error(
    sfclust:::log_mlik_each(1, membership, data, FALSE, FALSE, inla_args, model),
    "must be a formula"
  )
})

test_that("sfclust: within_model formula works with a direct INLA stack", {
  skip_if_not_installed("INLA")
  toy_data_list <- toy_data()

  stack_fun <- function(data) {
    INLA::inla.stack(
      data = list(y = data$y),
      A = list(1),
      effects = list(data.frame(intercept = rep(1, nrow(data)))),
      tag = "cluster"
    )
  }

  model <- sfclust_within_model(
    formula = y ~ 0 + intercept,
    family = "gaussian",
    stack_fun = stack_fun
  )

  fit <- suppressWarnings(suppressMessages(
    sfclust(
      toy_data_list$data,
      adjacency = toy_data_list$adjacency,
      fnames = "id_time",
      nclust = 1,
      within_model = model,
      correction = FALSE,
      niter = 1,
      move_prob = c(0, 0, 0, 1)
    )
  ))

  expect_s3_class(fit, "sfclust")
  expect_equal(deparse1(eval(attr(fit, "inla_args")$formula)), "y ~ 0 + intercept")
  expect_null(attr(fit, "fit_args")$stack_fun)
  expect_null(attr(fit, "fit_args")$fitted_fun)
})

test_that("stack_fun can provide a cluster-specific formula", {
  skip_if_not_installed("INLA")

  data <- data.frame(
    id = seq_len(9),
    ids = rep(seq_len(3), 3),
    sid = rep(seq_len(3), 3),
    id_time = rep(seq_len(3), each = 3),
    y = c(-0.3, 0.2, -0.1, 0.4, 0.1, 0.5, -0.2, 0.3, 0.6)
  )

  stack_fun <- function(data) {
    n_obs <- nrow(data)
    formula <- as.formula(
      paste0("y ~ 0 + intercept + f(u, model = 'iid', n = ", n_obs, ")")
    )

    list(
      stack = INLA::inla.stack(
        data = list(y = data$y),
        A = list(1),
        effects = list(data.frame(intercept = rep(1, n_obs), u = seq_len(n_obs))),
        tag = "cluster"
      ),
      formula = formula
    )
  }

  model <- sfclust_within_model(
    formula = NULL,
    family = "gaussian",
    stack_fun = stack_fun
  )

  models <- suppressWarnings(suppressMessages(
    sfclust:::log_mlik_all(
      c(1L, 1L, 2L),
      data,
      correction = FALSE,
      detailed = TRUE,
      inla_args = list(family = "gaussian"),
      within_model = model
    )
  ))

  formulas <- vapply(models, function(x) deparse1(x$.args$formula), character(1))
  expect_match(formulas[[1]], "n = 6")
  expect_match(formulas[[2]], "n = 3")
})


test_that("correction works with formula returned by stack_fun", {
  skip_if_not_installed("INLA")

  data <- data.frame(
    id = seq_len(9),
    ids = rep(seq_len(3), 3),
    sid = rep(seq_len(3), 3),
    id_time = rep(seq_len(3), each = 3),
    y = c(-0.3, 0.2, -0.1, 0.4, 0.1, 0.5, -0.2, 0.3, 0.6)
  )

  stack_fun <- function(data) {
    n_obs <- nrow(data)
    formula <- as.formula(
      "y ~ 0 + intercept + f(id_time, model = 'rw1', constr = TRUE)"
    )

    list(
      stack = INLA::inla.stack(
        data = list(y = data$y),
        A = list(1),
        effects = list(data.frame(intercept = rep(1, n_obs), id_time = data$id_time)),
        tag = "cluster"
      ),
      formula = formula
    )
  }

  model <- sfclust_within_model(
    formula = NULL,
    family = "gaussian",
    stack_fun = stack_fun
  )

  log_mlike <- suppressWarnings(suppressMessages(
    sfclust:::log_mlik_all(
      c(1L, 1L, 2L),
      data,
      correction = TRUE,
      detailed = FALSE,
      inla_args = list(family = "gaussian"),
      within_model = model
    )
  ))

  expect_length(log_mlike, 2)
  expect_true(all(is.finite(log_mlike)))
})

test_that("fitted: stack_fun requires fitted_fun", {
  toy_data_list <- toy_data()
  model <- sfclust_within_model(
    formula = y ~ 1,
    family = "gaussian",
    stack_fun = function(data) data
  )

  fit <- list(
    samples = list(
      membership = matrix(c(1L, 1L), nrow = 1),
      log_mlike = 0,
      move_counts = c(births = 0, deaths = 0, changes = 0, hypers = 0)
    ),
    clust = list(
      id = 1,
      membership = c(1L, 1L),
      models = list(list())
    )
  )
  attr(fit, "fit_args") <- list(
    data = transform(toy_data_list$data, sid = ids),
    within_model = model
  )
  attr(fit, "input_args") <- list(fnames = "id_time")
  attr(fit, "inla_args") <- list(formula = y ~ 1, family = "gaussian")
  class(fit) <- "sfclust"

  expect_error(
    fitted(fit),
    "require a user-supplied `fitted_fun`"
  )
})

test_that("update: within_model is preserved", {
  skip_if_not_installed("INLA")
  toy_data_list <- toy_data()
  model <- sfclust_within_model(formula = y ~ 1, family = "gaussian")

  fit <- suppressWarnings(suppressMessages(
    sfclust(
      toy_data_list$data,
      adjacency = toy_data_list$adjacency,
      fnames = "id_time",
      nclust = 1,
      within_model = model,
      correction = FALSE,
      niter = 1,
      move_prob = c(0, 0, 0, 1)
    )
  ))

  fit2 <- suppressWarnings(suppressMessages(update(fit, niter = 1, nmessage = 1)))
  expect_s3_class(attr(fit2, "fit_args")$within_model, "sfclust_within_model")
  expect_equal(attr(fit2, "fit_args")$within_model$family, "gaussian")
  expect_null(attr(fit2, "fit_args")$stack_fun)
  expect_null(attr(fit2, "fit_args")$fitted_fun)
})
