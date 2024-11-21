#' Print method for sfclust objects
#'
#' Prints details of an sfclust object, including the (i) within-cluster formula;
#' (ii) hyperparameters used for the MCMC sample such as the number of clusters penalty
#' (q) and the movement probabilities (move_prob); (iii) the number of movement type dones
#' during the MCMC sampling; and (iv) the log marginal likelihood of the model of the last
#' clustering sample.
#'
#' @param x An object of class 'sfclust'.
#' @param ... Additional arguments passed to `print.default`.
#' @export
print.sfclust <- function(x, ...) {
  cat("Within-cluster formula:\n")
  print(eval(attr(x, "inla_args")$formula), showEnv = FALSE, ...)

  cat("\nClustering hyperparameters:\n")
  hypernames <- c("q", "birth", "death", "change", "hyper")
  print(setNames(c(attr(x, "args")$q, attr(x, "args")$move_prob), hypernames), ...)

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
#' @param x A `sfclust` object.
#' @param niter An integer specifying the number of additional MCMC iterations to perform.
#' @param burnin An integer specifying the number of burn-in iterations to discard.
#' @param thin An integer specifying the thinning interval for recording results.
#' @param nmessage An integer specifying the number of messages to display during the process.
#' @param sample An integer specifying the clustering sample number to be executed.
#'        The default is the last sample (i.e., `nrow(x$samples$membership)`).
#' @param path_save A character string specifying the file path to save the results. If
#'        `NULL`, results are not saved.
#' @param nsave An integer specifying how often to save results. Defaults to `nmessage`.
#'
#' @details This function takes the last state of the Markov chain from a previous
#'          `sfclust` execution and uses it as the starting point for additional MCMC
#'          iterations. If `sample` is provided, it simply udpates the within-cluster
#'          models for the specified clustering `sample`.
#'
#' @method update sfclust
#' @export
update.sfclust <- function(x, niter = 100, burnin = 0, thin = 1, nmessage = 10, sample = NULL,
                           path_save = NULL, nsave = nmessage) {
  if (!is.null(sample)) {
    update_within(x, nrow(x$samples$membership))
  } else {
    update_sfclust(x, niter, burnin, thin, nmessage, path_save, nsave)
  }
}

update_sfclust <- function(x, niter = 100, burnin = 0, thin = 1, nmessage = 10,
                           path_save = NULL, nsave = nmessage) {
  nsamples <- nrow(x$samples$membership)

  args <- c(attr(x, "args"), attr(x, "inla_args"))
  args$graphdata$mst <- attr(x, "mst")[[nsamples]]
  args$graphdata$membership <- x$samples$membership[nsamples,]

  args$niter <- niter
  args$burnin <- burnin
  args$thin <- thin
  args$nmessage <- nmessage
  args$path_save <- path_save
  args$nsave <- nsave

  call <- as.call(c(list(as.name("sfclust")), args))
  eval(call, envir = parent.frame())
}

update_within <- function(x, sample = nrow(x$samples$membership)) {
  args <- attr(x, "inla_args")
  args$membership <- x$samples$membership[sample,]
  args$stdata <- attr(x, "args")$stdata
  args$stnames <- attr(x, "args")$stnames
  args$correction <- FALSE
  args$detailed <- TRUE

  call <- as.call(c(list(as.name("log_mlik_all")), args))
  x$clust$id <- sample
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
#'        level.
#' @param ... Additional arguments, currently not used.
#' @return A `stars` object with linear predictor fitted values at regions levels. In case
#'         `aggregate = TRUE`, the `output` produces an `stars` objecto at cluster levels.
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
#' \dontrun{
#'
#' # Assuming 'mod' is a pre-existing sfclust object
#'
#' # Estimated values ordering clusters by size
#' df_est <- fitted(mod, order = TRUE)
#'
#' # Estimated values aggregated by cluster
#' df_est <- fitted(mod, aggregate = TRUE)
#'
#' # Estimated values using a particular clustering sample
#' df_est <- fitted(mod, sample = 3)
#'
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

  membership <- object$samples$membership[sample,]
  if (sort) {
    membership <- sort_membership(object$samples$membership[sample,])
    object$clust$models <- object$clust$models[attr(membership, "order")]
  }

  # obtain fitted values
  clusters <- 1:max(membership)
  pred <- lapply(
    1:max(membership), linpred_each, membership, object$clust$models,
    attr(object, "args")$stdata, attr(object, "args")$stnames
  )
  pred <- do.call(rbind, pred)
  pred <- pred[order(pred$id), setdiff(names(pred), "id")]

  # save as stars object
  stars_obj <- attr(object, "args")$stdata[0]
  for (var_name in names(pred)) stars_obj[[var_name]] <- pred[[var_name]]

  # aggregate if required
  if (aggregate) {
    geom_clusters <- lapply(
      split(st_geometry(attr(object, "args")$stdata), membership),
      function(x) st_union(st_geometry(x))
    )
    geom_clusters <- do.call(c, geom_clusters)
    stars_obj <- aggregate(stars_obj[c("mean", "mean_cluster")], geom_clusters,
      join = st_within, FUN = mean)
  }

  stars_obj
}

linpred_each <- function(k, membership, models, stdata, stnames){
  df <- data_each(k, membership, stdata, stnames)[c("id")]
  df <- cbind(df, models[[k]]$summary.linear.predictor)
  df$lkd <- NULL
  df$mean_cluster <- linpred_each_corrected(models[[k]])
  df$cluster <- k
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
#' @param ... Additional arguments passed to the underlying plotting functions.
#'
#' @return A plot displaying the selected subgraphs as specified by `which`.
#'
#' @importFrom graphics par points
#' @export
plot.sfclust <- function(x, sample = x$clust$id, which = 1:3, clusters = NULL, sort = FALSE, legend = FALSE, ...) {

  nsamples <- nrow(x$samples$membership)
  if (sample < 1 || sample > nsamples) {
    stop("`sample` must be between 1 and the total number clustering (membership) samples.")
  }

  # get geometries, selected membership and clusters to plot
  geoms <- st_get_dimension_values(attr(x, "args")$stdata, attr(x, "args")$stnames[1])
  membership <- x$samples$membership[sample,]
  if (sort) membership <- sort_membership(x$samples$membership[sample,])
  if (is.null(clusters)) clusters <- 1:max(membership)

  # visualize
  which <- which[which %in% 1:3]
  if (length(which) > 1) opar <- par(mfrow = c(1, length(which)))
  if (1 %in% which) { # spatial clustering membership
    membership[!(membership %in% clusters)] <- NA
    membership <- factor(membership)
    plot(geoms, col = membership, border = "gray50",
      main = paste("Clustering sample", sample, "out of", nsamples),
      xlab = "Longitude", ylab = "Latitude", ...
    )
    if (legend) {
      legend("topright", legend = levels(membership), fill = 1:length(levels(membership)),
        border = "black", title = "Cluster", cex = 0.8, bty = "n")
    }
  }
  if (3 %in% which) { # functional shapes
    df <- fitted(x, sample = sample, sort = sort)
    df <- filter(df, !!as.name(attr(x, "args")$stnames[1]) %in% which(membership %in% clusters))
    df <- st_set_dimensions(df[c("cluster", "mean_cluster")], attr(x, "args")$stnames[1],
      values = seq_len(length(st_get_dimension_values(df, attr(x, "args")$stnames[1])))) |>
        as.data.frame()
    plot(df[[attr(x, "args")$stnames[2]]], df$mean_cluster, col = df$cluster,
      main = "Cluster mean functions", xlab = "Time", ylab = "Cluster linear predictor",
      pch = 19)
  }
  if (2 %in% which) { # marginal likelihood convergence
    plot(x$samples$log_mlike, type = "l", main = "Log marginal likelihood convergence",
      xlab = "Sample", ylab = "Log marginal likelihood", col = "blue", lwd = 2)
    points(sample, x$samples$log_mlike[sample], col = "red", pch = 19)
  }
  if (length(which) > 1) par(opar)

  invisible(NULL)
}
