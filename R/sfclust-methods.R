#' Print method for sfclust objects
#'
#' Prints details of an sfclust object, including the (i) within-cluster formula;
#' (ii) hyperparameters used for the MCMC sample such as the number of clusters penalty
#' (q) and the movement probabilities (move_prob); (iii) the number of movement type counts
#' during the MCMC sampling; and (iv) the log marginal likelihood of the model of the last
#' clustering sample.
#'
#' @param x An object of class 'sfclust'.
#' @param ... Additional arguments passed to `print.default`.
#'
#' @return Invisibly returns the input \code{sfclust} object \code{x}. The function also
#' prints a summary of:
#' \itemize{
#'   \item the within-cluster model formula,
#'   \item clustering hyperparameters,
#'   \item movement counts from the MCMC sampler,
#'   \item and the log marginal likelihood of the selected sample.
#' }
#'
#' @export
print.sfclust <- function(x, ...) {
  cat("Within-cluster formula:\n")
  print(eval(attr(x, "inla_args")$formula), showEnv = FALSE, ...)

  cat("\nClustering hyperparameters:\n")
  hypernames <- c("log(1-q)", "birth", "death", "change", "hyper")
  print(setNames(c(attr(x, "args")$logpen, attr(x, "args")$move_prob), hypernames), ...)

  cat("\nClustering movement counts:\n")
  print(x$samples$move_counts, ...)

  cat("\nLog marginal likelihood (sample ", x$clust$id, " out of ",
        length(x$samples$log_mlike), "): ", x$samples$log_mlike[x$clust$id], "\n", sep = "")

  invisible(x)
}


#' Summary method for sfclust objects
#'
#' This function summarizes the cluster assignments from the desired clustering `sample`.
#'
#' @param object An object of class 'sfclust'.
#' @param sample An integer specifying the clustering sample number to be summarized (default
#'        is the last sample).
#' @param sort Logical value indicating if clusters should be relabel based on number of
#'        elements.
#' @param ... Additional arguments passed to `print.default`.
#'
#' @return
#' Invisibly returns a table with the number of regions in each cluster for the selected
#' sample. The function also prints a summary that includes:
#' \itemize{
#'   \item the within-cluster model formula,
#'   \item the total number of MCMC clustering samples,
#'   \item the cluster membership counts for the specified sample (optionally sorted),
#'   \item and the log marginal likelihood of the selected clustering sample.
#' }
#'
#' @export
summary.sfclust <- function(object, sample = object$clust$id, sort = FALSE,...) {

  nsamples <- nrow(object$samples$membership)
  if (sample < 1 || sample > nsamples) {
    stop("`sample` must be between 1 and the total number clustering (membership) samples.")
  }

  cat("Summary for clustering sample", sample, "out of", nsamples, "\n")

  cat("\nWithin-cluster formula:\n")
  print(eval(attr(object, "inla_args")$formula), showEnv = FALSE, ...)

  cat("\nCounts per cluster:")
  membership <- object$samples$membership[sample,]
  if (sort) membership <- sort_membership(object$samples$membership[sample,])
  cluster_summary <- table(membership, deparse.level = 0)
  print(cluster_summary, ...)

  cat("\nLog marginal likelihood: ", object$samples$log_mlike[sample], "\n")
  invisible(cluster_summary)
}

# relabel membership based on the frequency of each cluster
sort_membership <- function(x) {
  clusters_sorted <- order(table(x), decreasing = TRUE)
  clusters_labels <- setNames(seq_along(clusters_sorted), clusters_sorted)
  x <- as.integer(clusters_labels[as.character(x)])
  attr(x, "order") <- clusters_sorted
  x
}

# TRUE when result carries stars metadata (set by sfclust() wrapper)
has_stars_metadata <- function(x) !is.null(attr(x, "stdata"))

#' Update MCMC Clustering Procedure
#'
#' This function continues the MCMC sampling of a `sfclust` object based on previous results or
#' update the model fitting for a specified sample clustering if the argument `sample` is
#' provided.
#'
#' @param object A `sfclust` object.
#' @param niter An integer specifying the number of additional MCMC iterations to perform.
#' @param burnin An integer specifying the number of burn-in iterations to discard.
#' @param thin An integer specifying the thinning interval for recording results.
#' @param nmessage An integer specifying the number of messages to display during the process.
#' @param sample An integer specifying the clustering sample number to be executed.
#'        The default is the last sample (i.e., `nrow(x$samples$membership)`).
#' @param path_save A character string specifying the file path to save the results. If
#'        `NULL`, results are not saved.
#' @param nsave An integer specifying how often to save results. Defaults to `nmessage`.
#' @param ... Additional arguments (currently not used).
#'
#' @details This function takes the last state of the Markov chain from a previous
#'          `sfclust_fit` / `sfclust` execution and uses it as the starting point for
#'          additional MCMC iterations. If `sample` is provided, it simply updates the
#'          within-cluster models for the specified clustering `sample`.
#'
#' @return An updated `sfclust` object with (i) new clustering samples if `sample` is not
#'   specified, or (ii) updated within-cluster model results if `sample` is given.
#'
#' @importFrom stats update
#' @method update sfclust
#' @export
update.sfclust <- function(object, niter = 100, burnin = 0, thin = 1, nmessage = 10, sample = NULL,
                           path_save = NULL, nsave = nmessage, ...) {
  if (!is.null(sample)) {
    update_within(object, sample)
  } else {
    update_sfclust(object, niter, burnin, thin, nmessage, path_save, nsave)
  }
}

update_sfclust <- function(x, niter = 100, burnin = 0, thin = 1, nmessage = 10,
                           path_save = NULL, nsave = nmessage) {
  nsamples <- nrow(x$samples$membership)

  args <- c(attr(x, "args"), attr(x, "inla_args"))
  args$fun_col <- NULL
  args$graphdata$mst        <- attr(x, "mst")[[nsamples]]
  args$graphdata$membership <- x$samples$membership[nsamples, ]

  args$niter     <- niter
  args$burnin    <- burnin
  args$thin      <- thin
  args$nmessage  <- nmessage
  args$path_save <- path_save
  args$nsave     <- nsave

  call   <- as.call(c(list(as.name("sfclust_fit")), args))
  result <- eval(call, envir = parent.frame())

  if (has_stars_metadata(x)) {
    attr(result, "stdata")       <- attr(x, "stdata")
    attr(result, "sp_dims")      <- attr(x, "sp_dims")
    attr(result, "fun_dims")     <- attr(x, "fun_dims")
    attr(result, "valid_ids")    <- attr(x, "valid_ids")
    attr(result, "args")$fun_col <- attr(x, "args")$fun_col
  }

  result
}

update_within <- function(x, sample = nrow(x$samples$membership)) {
  args             <- attr(x, "inla_args")
  args$data        <- attr(x, "args")$data
  args$membership  <- x$samples$membership[sample, ]
  args$correction  <- FALSE
  args$detailed    <- TRUE

  call <- as.call(c(list(as.name("log_mlik_all")), args))
  x$clust$id     <- sample
  x$clust$models <- eval(call, envir = parent.frame())
  x
}

#' Fitted Values for `sfclust` Objects
#'
#' This function calculates the fitted values for a specific clustering sample in an
#' `sfclust` object, based on the estimated models for each cluster. The fitted
#' values are computed using the membership assignments and model parameters
#' associated with the selected clustering sample.
#'
#' @param object An object of class 'sfclust', containing clustering results and models.
#' @param sample An integer specifying the clustering sample number for which
#'        the fitted values should be computed. The default is the `id` of the
#'        current clustering. The value must be between 1 and the total number
#'        of clustering (membership) samples.
#' @param sort Logical value indicating if clusters should be relabel based on number of
#'        elements.
#' @param aggregate Logical value indicating if fitted values are desired at cluster
#'        level. Only supported for vector geometry data (stars wrapper results).
#' @param ... Additional arguments, currently not used.
#' @return For stars wrapper results: a `stars` object with linear predictor fitted values
#'         at region level. For core (`sfclust_fit`) results: a data frame. In case
#'         `aggregate = TRUE` (stars only), the output is a `stars` object at cluster level.
#'
#' @details
#' The function first checks if the provided `sample` value is valid (i.e., it is
#' within the range of available clustering samples). If the specified `sample`
#' does not match the current clustering `id`, the `sfclust` object is updated
#' accordingly. It then retrieves the membership assignments and cluster models
#' for the selected sample, calculates the linear predictions for each cluster,
#' and combines them into a matrix of fitted values.
#'
#' @examples
#'
#' \donttest{
#' library(sfclust)
#'
#' data(stgaus)
#' result <- sfclust(stgaus, formula = y ~ f(idf_time, model = "rw1"), niter = 10,
#'   nmessage = 1)
#'
#' # Estimated values ordering clusters by size
#' df_est <- fitted(result, sort = TRUE)
#'
#' # Estimated values aggregated by cluster
#' df_est <- fitted(result, aggregate = TRUE)
#'
#' # Estimated values using a particular clustering sample
#' df_est <- fitted(result, sample = 3)
#' }
#'
#' @importFrom stats fitted
#' @importFrom sf st_within st_union
#' @export
fitted.sfclust <- function(object, sample = object$clust$id, sort = FALSE, aggregate = FALSE, ...) {
  if (sample < 1 || sample > nrow(object$samples$membership)) {
    stop("`sample` must be between 1 and the total number clustering (membership) samples.")
  }
  if (sample != object$clust$id) object <- update_within(object, sample = sample)

  membership <- object$samples$membership[sample, ]
  if (sort) {
    membership <- sort_membership(object$samples$membership[sample, ])
    object$clust$models <- object$clust$models[attr(membership, "order")]
  }

  data <- attr(object, "args")$data

  pred <- lapply(1:max(membership), linpred_each,
                 membership, object$clust$models, data)
  pred <- do.call(rbind, pred)
  pred <- pred[order(pred$id), setdiff(names(pred), "id")]

  if (!has_stars_metadata(object)) {
    return(pred)
  }

  # reconstruct stars object from data frame
  stdata_ref <- attr(object, "stdata")
  sp_dims    <- attr(object, "sp_dims")
  stars_obj  <- stdata_ref[0]
  valid_ids  <- attr(object, "valid_ids")

  if (!is.null(valid_ids)) {
    n_sp     <- prod(dim(stdata_ref)[sp_dims])
    n_fun    <- prod(dim(stdata_ref)[setdiff(names(dim(stdata_ref)), sp_dims)])
    n_valid  <- length(valid_ids)
    full_idx <- rep(valid_ids, times = n_fun) +
                rep((seq_len(n_fun) - 1L) * n_sp, each = n_valid)
    for (var_name in names(pred)) {
      full_vec <- pred[[var_name]][rep(NA_integer_, n_sp * n_fun)]
      full_vec[full_idx] <- pred[[var_name]]
      stars_obj[[var_name]] <- full_vec
    }
  } else {
    for (var_name in names(pred)) stars_obj[[var_name]] <- pred[[var_name]]
  }

  # aggregate if required (vector geometry only)
  if (aggregate) {
    if (length(sp_dims) > 1) {
      warning("`aggregate = TRUE` is not supported for raster data.", call. = FALSE)
    } else {
      geom_clusters <- lapply(
        split(st_geometry(stdata_ref), membership),
        function(x) st_union(st_geometry(x))
      )
      geom_clusters <- do.call(c, geom_clusters)
      stars_obj <- aggregate(stars_obj[c("mean_cluster", "mean_cluster_inv")], geom_clusters,
        join = st_within, FUN = mean)
    }
  }

  stars_obj
}

linpred_each <- function(k, membership, models, data) {
  # get inverse of linear predictor
  link_name <- tolower(models[[k]]$misc$linkfunctions$names)
  inv_link  <- eval(parse(text = paste0("INLA::inla.link.inv", link_name)))

  # subset data to cluster k; row order matches INLA's linear predictor
  df <- data[data$ids %in% which(membership == k), , drop = FALSE]

  df <- cbind(df, models[[k]]$summary.linear.predictor)
  df$kld      <- NULL
  df$mean_inv <- inv_link(df$mean)

  df$cluster          <- k
  df$mean_cluster     <- linpred_each_corrected(models[[k]])
  df$mean_cluster_inv <- inv_link(df$mean_cluster)
  df
}

# linear predictor only for terms that are defined at cluster level
linpred_each_corrected <- function(x){
  aux <- if ("id" %in% names(x$summary.random)) x$summary.random[["id"]]$mean else 0
  x$summary.linear.predictor$mean - aux
}

#' Plot function for `sfclust` objects
#'
#' This function visualizes the estimated clusters from an `sfclust` object. It can display:
#' (1) a map of regions colored by their assigned cluster,
#' (2) the functional shapes of the linear predictors for each cluster,
#' and (3) a traceplot of the log marginal likelihood.
#' A conditional legend is added if the number of clusters is less than 10.
#'
#' @param x An `sfclust` object containing the clustering results, including the cluster assignments and model parameters.
#' @param sample Integer specifying the clustering sample number to summarize. Defaults to the last sample.
#' @param which Integer vector indicating which plot to display. Options are:
#'        - 1: Map of regions colored by cluster assignment.
#'        - 2: Functional shapes of the linear predictors for each cluster.
#'        - 3: Traceplot of the log marginal likelihood.
#' @param clusters Optional vector specifying which clusters to plot. If `NULL`, all clusters are included.
#' @param sort Logical value indicating whether clusters should be relabeled based on the number of elements. Default is `FALSE`.
#' @param legend Logical value indicating whether a legend should be included in the plot. Default is `FALSE`.
#' @param geom_before An optional `ggplot2` geom layer to add before the cluster map layer (e.g., a background tile). Default is `NULL`.
#' @param ... Additional arguments passed to the underlying plotting functions.
#'
#' @return A composed `patchwork` object displaying the selected subgraphs as specified by `which`.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_continuous scale_y_continuous element_blank
#' @importFrom ggplot2 theme_bw theme_minimal theme labs geom_sf geom_raster
#' @importFrom patchwork wrap_plots
#' @importFrom sf st_sf
#' @export
plot.sfclust <- function(x, sample = x$clust$id, which = 1:3, clusters = NULL, sort = FALSE,
                         legend = FALSE, geom_before = NULL, ...) {

  figs <- list(gg1 = NA, gg2 = NA, gg3 = NA)
  if (1 %in% which) {
    figs$gg1 <- plot_clusters_map(x, sample, clusters, sort, legend, geom_before)
  }
  if (2 %in% which) {
    if (!legend || (1 %in% which)) legend = FALSE
    figs$gg2 <- plot_clusters_fitted(x, sample, clusters, sort, legend)
  }
  if (3 %in% which) {
    figs$gg3 <- plot_log_mlik(x, sample)
  }
  wrap_plots(figs[which])
}

#' Plot a spatial map of cluster assignments
#'
#' Produces a `ggplot2` map of spatial regions colored by their cluster assignment for a
#' given MCMC sample of an `sfclust` object. Requires the result to have been produced
#' by [sfclust()] (stars wrapper); not available for [sfclust_fit()] results.
#'
#' @param x An `sfclust` object.
#' @param sample Integer specifying the clustering sample to display. Defaults to the last sample.
#' @param clusters Optional vector of cluster IDs to include. If `NULL`, all clusters are shown.
#' @param sort Logical; if `TRUE`, clusters are relabeled by decreasing size. Default is `FALSE`.
#' @param legend Logical; if `TRUE`, a fill legend is included. Default is `FALSE`.
#' @param geom_before An optional `ggplot2` geom layer to add before the cluster fill layer. Default is `NULL`.
#' @param ... Additional arguments passed to `geom_sf()`.
#'
#' @return A `ggplot2` object.
#'
#' @importFrom ggplot2 facet_wrap
#' @export
plot_clusters_map <- function(x, sample = x$clust$id, clusters = NULL, sort = FALSE, legend = FALSE, geom_before = NULL, ...) {
  if (!has_stars_metadata(x)) {
    stop("`plot_clusters_map()` requires a stars-based result from `sfclust()`. ",
         "For custom data, build the map directly from `fitted()` output.")
  }
  nsamples <- check_sample_and_get_nsample(x, sample)
  aux      <- get_membership_and_clusters(x, sample, sort, clusters)
  stdata   <- attr(x, "stdata")
  sp_dims  <- attr(x, "sp_dims")

  membership <- aux$membership
  membership[!(membership %in% aux$clusters)] <- NA
  membership <- factor(membership)

  if (length(sp_dims) == 1) {
    geoms <- st_get_dimension_values(stdata, sp_dims)
    gg <- ggplot(st_sf(geometry = geoms, membership = membership))
    if (!is.null(geom_before)) gg <- gg + geom_before
    gg <- gg + geom_sf(aes(fill = membership), shape = 21, ...)
  } else {
    dims    <- expand_dimensions(stdata)
    x_vals  <- dims[[sp_dims[1]]]
    y_vals  <- dims[[sp_dims[2]]]
    sp_grid <- expand.grid(x_vals, y_vals)
    names(sp_grid) <- sp_dims
    valid_ids <- attr(x, "valid_ids")
    if (!is.null(valid_ids)) {
      full_membership <- rep(NA_integer_, nrow(sp_grid))
      full_membership[valid_ids] <- membership
      membership <- factor(full_membership)
    }
    sp_grid$membership <- membership
    gg <- ggplot(sp_grid, aes(!!as.name(sp_dims[1]), !!as.name(sp_dims[2]), fill = membership)) +
      geom_raster() +
      ggplot2::coord_equal()
    if (!is.null(geom_before)) gg <- gg + geom_before
  }
  gg <- gg +
    labs(fill = NULL, subtitle = paste("Clustering:", sample, "/", nsamples)) +
    theme_bw()
  if (!legend) gg <- gg + theme(legend.position = "none")
  gg
}

#' Plot functional shapes of cluster linear predictors
#'
#' Plots the estimated mean functional shape (linear predictor or inverse-link scale) for
#' each cluster in a given MCMC sample of an `sfclust` object.
#'
#' @param x An `sfclust` object.
#' @param sample Integer specifying the clustering sample to display. Defaults to the last sample.
#' @param clusters Optional vector of cluster IDs to include. If `NULL`, all clusters are shown.
#' @param sort Logical; if `TRUE`, clusters are relabeled by decreasing size. Default is `FALSE`.
#' @param legend Logical; if `TRUE`, a color legend is included. Default is `FALSE`.
#' @param inv_link Logical; if `TRUE` (default), values are shown on the inverse-link (mean) scale.
#' @param fun_col Character. Name of the column to use as the x-axis (functional dimension).
#'        If `NULL`, taken from the result's stored args (set automatically by [sfclust()]).
#' @param ... Additional arguments passed to `geom_line()`.
#'
#' @return A `ggplot2` object.
#'
#' @export
plot_clusters_fitted <- function(x, sample = x$clust$id, clusters = NULL, sort = FALSE, legend = FALSE, inv_link = TRUE, fun_col = NULL, ...) {
  nsamples <- check_sample_and_get_nsample(x, sample)
  aux      <- get_membership_and_clusters(x, sample, sort, clusters)
  if (is.null(fun_col)) fun_col <- attr(x, "args")$fun_col
  varname  <- if (!inv_link) "mean_cluster" else "mean_cluster_inv"

  df <- as.data.frame(fitted(x, sample = sample, sort = sort))
  df <- unique(df[df$cluster %in% aux$clusters, c(fun_col, varname, "cluster")])

  gg <- ggplot(df) +
    geom_line(aes(!!as.name(fun_col), !!as.name(varname), color = factor(cluster)), ...) +
    labs(x = NULL, y = "Estimated mean", subtitle = "Cluster functions", color = NULL) +
    theme_bw()
  if (!legend) gg <- gg + theme(legend.position = "none")
  gg
}

#' Plot log marginal likelihood convergence trace
#'
#' Plots the log marginal likelihood across MCMC samples for an `sfclust` object,
#' highlighting the selected sample.
#'
#' @param x An `sfclust` object.
#' @param sample Integer specifying the clustering sample to highlight. Defaults to the last sample.
#' @param ... Additional arguments passed to `geom_line()`.
#'
#' @return A `ggplot2` object.
#'
#' @export
plot_log_mlik <- function(x, sample = x$clust$id, ...) {
  nsamples <- check_sample_and_get_nsample(x, sample)

  gg <- ggplot(mapping = aes(sample, log_mlike)) +
    geom_line(data = data.frame(sample = 1:nsamples, log_mlike = x$samples$log_mlike), ...) +
    geom_point(data = data.frame(sample = sample, log_mlike = x$samples$log_mlike[sample]),
      color = 2) +
    labs(x = "Sample", y = "Log marginal likelihood", subtitle = "Convergence") +
    theme_bw()
  gg
}

#' Plot observed time series by cluster
#'
#' Plots individual region time series faceted by cluster, overlaid with the cluster mean,
#' for a given variable in an `sfclust` object.
#'
#' @param x An `sfclust` object.
#' @param var An unquoted variable name to plot on the y-axis. For stars-based results this
#'        refers to a variable in the original `stars` object; for core results it refers to
#'        a column in the long-format data stored in the result.
#' @param clusters Optional vector of cluster IDs to include. If `NULL`, all clusters are shown.
#' @param sort Logical; if `TRUE`, clusters are relabeled by decreasing size. Default is `FALSE`.
#' @param fun_col Character. Name of the column to use as the x-axis (functional dimension).
#'        If `NULL`, taken from the result's stored args (set automatically by [sfclust()]).
#' @param ... Additional arguments passed to `geom_line()` for individual region series.
#'
#' @return A `ggplot2` object with one facet per cluster.
#'
#' @export
plot_clusters_series <- function(x, var, clusters = NULL, sort = FALSE, fun_col = NULL, ...) {

  if (is.null(fun_col)) fun_col <- attr(x, "args")$fun_col
  fitted_df <- as.data.frame(fitted(x, sort = sort))

  if (has_stars_metadata(x)) {
    stdata  <- attr(x, "stdata")
    sp_dims <- attr(x, "sp_dims")

    stdata$cluster <- fitted(x, sort = sort)$cluster
    if (is.null(clusters)) clusters <- 1:max(stdata$cluster, na.rm = TRUE)

    if (length(sp_dims) == 1) {
      ns <- length(st_get_dimension_values(stdata, sp_dims))
      auxdata <- stdata |>
        st_set_dimensions(sp_dims, values = 1:ns) |>
        as.data.frame()
      sp_id_col <- sp_dims
    } else {
      auxdata <- as.data.frame(stdata)
      nx <- length(unique(auxdata[[sp_dims[1]]]))
      iy_idx <- as.integer(factor(auxdata[[sp_dims[2]]]))
      ix_idx <- as.integer(factor(auxdata[[sp_dims[1]]]))
      auxdata$.sp_id <- (iy_idx - 1L) * nx + ix_idx
      sp_id_col <- ".sp_id"
    }
  } else {
    data    <- attr(x, "args")$data
    auxdata <- data
    auxdata$cluster <- fitted_df$cluster[match(data$id, fitted_df$id)]
    if (is.null(clusters)) clusters <- 1:max(auxdata$cluster, na.rm = TRUE)
    sp_id_col <- "ids"
  }

  fun_sym  <- as.name(fun_col)
  stcluster <- auxdata |>
    dplyr::group_by(!!fun_sym, cluster) |>
    dplyr::summarise(mean_cluster = mean({{ var }}), .groups = "drop")

  dplyr::filter(auxdata, cluster %in% clusters) |>
    ggplot() +
      geom_line(aes(!!fun_sym, {{ var }}, group = !!as.name(sp_id_col)), color = "gray50", linewidth = 0.3, ...) +
      geom_line(aes(!!fun_sym, mean_cluster), dplyr::filter(stcluster, cluster %in% clusters), color = "red", linewidth = 0.4) +
      facet_wrap(~ cluster) +
      theme_bw() +
      labs(x = NULL)
}


check_sample_and_get_nsample <- function(x, sample) {
  nsamples <- nrow(x$samples$membership)
  if (sample < 1 || sample > nsamples) {
    stop("`sample` must be between 1 and the total number clustering (membership) samples.")
  }
  return(nsamples)
}

get_membership_and_clusters <- function(x, sample, sort = FALSE, clusters = NULL) {
  membership <- x$samples$membership[sample,]
  if (sort) membership <- sort_membership(membership)
  if (is.null(clusters)) clusters <- 1:max(membership)
  list(membership = membership, clusters = clusters)
}
