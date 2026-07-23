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
  print(setNames(c(attr(x, "fit_args")$logpen, attr(x, "fit_args")$move_prob), hypernames), ...)

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
  nsamples <- x$clust$id
  fit_args <- attr(x, "fit_args")

  fit_args$graphdata$mst        <- attr(x, "mst")[[nsamples]]
  fit_args$graphdata$membership <- x$samples$membership[nsamples, ]
  fit_args$niter     <- niter
  fit_args$burnin    <- burnin
  fit_args$thin      <- thin
  fit_args$nmessage  <- nmessage
  fit_args$path_save <- path_save
  fit_args$nsave     <- nsave
  fit_args$inla_args  <- attr(x, "inla_args")
  fit_args$save_class <- class(x)
  fit_args$input_args <- attr(x, "input_args")

  do.call(sfclust_fit, fit_args)
}

update_within <- function(x, sample = nrow(x$samples$membership)) {
  x$clust$id     <- sample
  x$clust$models <- log_mlik_all(x$samples$membership[sample, ], attr(x, "fit_args")$data,
                                 FALSE, TRUE, attr(x, "inla_args"))
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
#'        level. Only supported for `sfclust_stars` results.
#' @param ... Additional arguments, currently not used.
#' @return A data frame with fitted values and cluster assignments, keyed by `id`.
#'         For `sfclust_stars` objects, a `stars` object is returned instead.
#'
#' @examples
#'
#' \donttest{
#' library(sfclust)
#'
#' data(stgaus)
#' result <- sfclust(stgaus, formula = y ~ f(id_time, model = "rw1"), niter = 10,
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

  # obtain fitted values
  pred <- lapply(1:max(membership), linpred_each, membership,
                 object$clust$models, attr(object, "fit_args")$data)
  pred <- do.call(rbind, pred)
  pred <- pred[order(pred$id), ]

  if (aggregate) {
    fnames <- attr(object, "input_args")$fnames
    pred <- unique(pred[c("cluster", fnames, "mean_cluster", "mean_cluster_inv")])
  }

  pred
}

#' @importFrom sf st_within st_union
#' @export
fitted.sfclust_stars <- function(object, sample = object$clust$id, sort = FALSE, aggregate = FALSE, ...) {
  pred <- NextMethod(aggregate = FALSE)

  stars_obj <- attr(object, "input_args")$stars
  spnames   <- attr(object, "input_args")$spnames
  id       <- attr(object, "fit_args")$data$id

  # fill variables into a stars object
  n_total  <- prod(dim(stars_obj))
  for (var_name in names(pred)) {
    full_vec <- rep(NA_real_, n_total)
    full_vec[id] <- pred[[var_name]]
    stars_obj[[var_name]] <- full_vec
  }

  if (aggregate) {
    membership <- object$samples$membership[sample, ]
    if (sort) membership <- sort_membership(membership)
    if (length(spnames) > 1) {
      warning("`aggregate = TRUE` is not supported for raster data.", call. = FALSE)
    } else {
    # Project to planar CRS if needed so st_union/st_within use GEOS (no lon/lat warning)
      crs_orig  <- sf::st_crs(stars_obj)
      is_lonlat <- isTRUE(sf::st_is_longlat(stars_obj))
      stars_work <- if (is_lonlat)
        sf::st_transform(stars_obj, sf::st_crs("+proj=laea +datum=WGS84"))
      else
        stars_obj
      geom_clusters <- lapply(
        split(st_geometry(stars_work), membership),
        function(x) st_union(st_geometry(x))
      )
      geom_clusters <- do.call(c, geom_clusters)
      stars_obj <- aggregate(stars_work[c("mean_cluster", "mean_cluster_inv")], geom_clusters,
        join = st_within, FUN = mean)
      if (is_lonlat) stars_obj <- sf::st_transform(stars_obj, crs_orig)
    }
  }

  stars_obj
}

linpred_each <- function(k, membership, models, data) {
  # get inverse of linear predictor
  link_name <- tolower(models[[k]]$misc$linkfunctions$names)
  inv_link  <- eval(parse(text = paste0("INLA::inla.link.inv", link_name)))

  # linear predictor per region
  df <- data[data$sid %in% which(membership == k), , drop = FALSE]
  df <- cbind(df, models[[k]]$summary.linear.predictor)
  df$kld      <- NULL
  df$mean_inv <- inv_link(df$mean)

  # linear predictor per cluster
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
#' Visualizes fitted cluster functions (plot 1) and a log marginal likelihood traceplot
#' (plot 2). For `sfclust_stars` objects, also includes a spatial map of cluster
#' assignments (plot 1 in that method, shifting the others to 2 and 3).
#'
#' @param x An `sfclust` object.
#' @param sample Integer specifying the clustering sample to display. Defaults to the last sample.
#' @param which Integer vector indicating which plots to show. For `sfclust`: 1 = cluster
#'        functions, 2 = log marginal likelihood. For `sfclust_stars`: 1 = map, 2 = cluster
#'        functions, 3 = log marginal likelihood.
#' @param clusters Optional vector of cluster IDs to include. If `NULL`, all clusters are shown.
#' @param sort Logical; if `TRUE`, clusters are relabeled by decreasing size. Default is `FALSE`.
#' @param legend Logical; if `TRUE`, a legend is included. Default is `FALSE`.
#' @param fnames Character. Column name for the x-axis of cluster function plots.
#'        If `NULL`, taken from the result's stored args.
#' @param ... Additional arguments passed to underlying plot functions.
#'
#' @return A composed `patchwork` object.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_continuous scale_y_continuous element_blank
#' @importFrom ggplot2 theme_bw theme_minimal theme labs geom_sf geom_raster
#' @importFrom patchwork wrap_plots
#' @importFrom sf st_sf
#' @export
plot.sfclust <- function(x, sample = x$clust$id, which = 1:2, clusters = NULL, sort = FALSE,
                         legend = FALSE, fnames = NULL, ...) {
  figs <- list(gg1 = NA, gg2 = NA)
  if (1 %in% which) {
    figs$gg1 <- plot_clusters_fitted(x, sample, clusters, sort, legend, fnames = fnames)
  }
  if (2 %in% which) {
    figs$gg2 <- plot_log_mlik(x, sample)
  }
  wrap_plots(figs[which])
}

#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_continuous scale_y_continuous element_blank
#' @importFrom ggplot2 theme_bw theme_minimal theme labs geom_sf geom_raster
#' @importFrom patchwork wrap_plots
#' @importFrom sf st_sf
#' @export
plot.sfclust_stars <- function(x, sample = x$clust$id, which = 1:3, clusters = NULL, sort = FALSE,
                               legend = FALSE, geom_before = NULL, fnames = NULL, ...) {
  figs <- list(gg1 = NA, gg2 = NA, gg3 = NA)
  if (1 %in% which) {
    figs$gg1 <- plot_clusters_map(x, sample, clusters, sort, legend, geom_before)
  }
  if (2 %in% which) {
    if (!legend || (1 %in% which)) legend <- FALSE
    figs$gg2 <- plot_clusters_fitted(x, sample, clusters, sort, legend, fnames = fnames)
  }
  if (3 %in% which) {
    figs$gg3 <- plot_log_mlik(x, sample)
  }
  wrap_plots(figs[which])
}

#' Plot a spatial map of cluster assignments
#'
#' Produces a `ggplot2` map of spatial regions colored by their cluster assignment for a
#' given MCMC sample of an `sfclust_stars` object.
#'
#' @param x An `sfclust_stars` object.
#' @param sample Integer specifying the clustering sample to display. Defaults to the last sample.
#' @param clusters Optional vector of cluster IDs to include. If `NULL`, all clusters are shown.
#' @param sort Logical; if `TRUE`, clusters are relabeled by decreasing size. Default is `FALSE`.
#' @param legend Logical; if `TRUE`, a fill legend is included. Default is `FALSE`.
#' @param geom_before An optional `ggplot2` geom layer to add before the cluster fill layer. Default is `NULL`.
#' @param ... Additional arguments passed to `geom_sf()`.
#'
#' @return A `ggplot2` object.
#'
#' @importFrom stars geom_stars st_as_stars st_dimensions
#' @importFrom ggplot2 facet_wrap
#' @export
plot_clusters_map <- function(x, sample = x$clust$id, clusters = NULL, sort = FALSE, legend = FALSE, geom_before = NULL, ...) {
  if (!inherits(x, "sfclust_stars")) {
    stop("`plot_clusters_map()` requires an `sfclust_stars` object from `sfclust()`.")
  }
  nsamples <- check_sample_and_get_nsample(x, sample)
  stars_obj <- attr(x, "input_args")$stars
  spnames   <- attr(x, "input_args")$spnames

  # filter memebership if required
  aux      <- get_membership_and_clusters(x, sample, sort, clusters)
  membership <- aux$membership
  membership[!(membership %in% aux$clusters)] <- NA

  # create spatial only stars object
  sp_dims   <- st_dimensions(stars_obj)[spnames]
  mem_array <- array(NA_integer_, dim(sp_dims))
  mem_array[attr(x, "fit_args")$graphdata$valid_ids] <- membership
  sp_stars <- st_as_stars(membership = mem_array, dimensions = sp_dims)
  sp_stars$membership <- factor(sp_stars$membership)

  # visualize
  geom_vals <- st_get_dimension_values(sp_stars, spnames[[1]])
  is_point  <- inherits(geom_vals, "sfc_POINT") || inherits(geom_vals, "sfc_MULTIPOINT")

  gg <- ggplot()
  if (!is.null(geom_before)) gg <- gg + geom_before
  if (is_point) {
    gg <- gg + geom_stars(data = sp_stars, aes(fill = membership), shape = 21, color = "black", ...)
  } else {
    gg <- gg + geom_stars(data = sp_stars, aes(fill = membership), ...)
  }
  gg <- gg +
      labs(color = NULL, fill = NULL, subtitle = paste("Clustering:", sample, "/", nsamples)) +
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
#' @param fnames Character. Name of the column to use as the x-axis (functional dimension).
#'        If `NULL`, taken from the result's stored args (set automatically by [sfclust()]).
#' @param ... Additional arguments passed to `geom_line()`.
#'
#' @return A `ggplot2` object.
#'
#' @export
plot_clusters_fitted <- function(x, sample = x$clust$id, clusters = NULL, sort = FALSE, legend = FALSE, inv_link = TRUE, fnames = NULL, ...) {
  nsamples <- check_sample_and_get_nsample(x, sample)
  aux      <- get_membership_and_clusters(x, sample, sort, clusters)
  if (is.null(fnames)) fnames <- resolve_fnames(x)
  varname  <- if (!inv_link) "mean_cluster" else "mean_cluster_inv"

  df <- as.data.frame(fitted(x, sample = sample, sort = sort))
  df <- unique(df[df$cluster %in% aux$clusters, c(fnames, varname, "cluster")])

  gg <- ggplot(df) +
    geom_line(aes(!!as.name(fnames), !!as.name(varname), color = factor(cluster)), ...) +
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
#' @param var An unquoted variable name to plot on the y-axis.
#' @param clusters Optional vector of cluster IDs to include. If `NULL`, all clusters are shown.
#' @param sort Logical; if `TRUE`, clusters are relabeled by decreasing size. Default is `FALSE`.
#' @param fnames Character. Name of the column to use as the x-axis (functional dimension).
#'        If `NULL`, taken from the result's stored args (set automatically by [sfclust()]).
#' @param ... Additional arguments passed to `geom_line()` for individual region series.
#'
#' @return A `ggplot2` object with one facet per cluster.
#'
#' @export
plot_clusters_series <- function(x, var, clusters = NULL, sort = FALSE, fnames = NULL, ...) {
  UseMethod("plot_clusters_series")
}

#' @export
plot_clusters_series.sfclust <- function(x, var, clusters = NULL, sort = FALSE, fnames = NULL, ...) {
  if (is.null(fnames)) fnames <- resolve_fnames(x)

  # data with clusters
  fitted_df <- fitted.sfclust(x, sort = sort)
  auxdata <- attr(x, "fit_args")$data
  auxdata$cluster <- fitted_df$cluster[match(auxdata$id, fitted_df$id)]
  if (is.null(clusters)) clusters <- 1:max(auxdata$cluster, na.rm = TRUE)

  # data summary per group
  fun_sym   <- as.name(fnames)
  stcluster <- auxdata |>
    dplyr::group_by(!!fun_sym, cluster) |>
    dplyr::summarise(mean_cluster = mean({{ var }}), .groups = "drop")

  # visualize
  dplyr::filter(auxdata, cluster %in% clusters) |>
    ggplot() +
      geom_line(aes(!!fun_sym, {{ var }}, group = !!as.name("ids")), color = "gray50", linewidth = 0.3, ...) +
      geom_line(aes(!!fun_sym, mean_cluster), dplyr::filter(stcluster, cluster %in% clusters), color = "red", linewidth = 0.4) +
      facet_wrap(~ cluster) +
      theme_bw() +
      labs(x = NULL)
}

resolve_fnames <- function(x) {
  fnames <- attr(x, "input_args")$fnames
  if (length(fnames) != 1L)
    stop("Requires exactly one functional dimension (", length(fnames), " found). Set `fnames`.")
  fnames
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
