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

  cat("\nLog marginal likelihood (sample ", x$clustering$id, " out of ",
        length(x$samples$log_mlike), "): ", x$samples$log_mlike[x$clustering$id], "\n", sep = "")

  invisible(x)
}


#' Summary method for sfclust objects
#'
#' This function summarizes the cluster assignments from the desired clustering `sample`.
#'
#' @param x An object of class 'sfclust'.
#' @param sample An integer specifying the clustering sample number to be summarized (default
#' is the last sample).
#' @param ... Additional arguments passed to `print.default`.
#' @export
summary.sfclust <- function(x, sample = x$clustering$id, sort = FALSE,...) {

  nsamples <- nrow(x$samples$membership)
  if (sample < 1 || sample > nsamples) {
    stop("`sample` must be between 1 and the total number clustering (membership) samples.")
  }

  cat("Summary for clustering sample", sample, "out of", nsamples, "\n")

  cat("\nWithin-cluster formula:\n")
  print(eval(attr(x, "inla_args")$formula), showEnv = FALSE, ...)

  cat("\nCounts per cluster:")
  membership <- x$samples$membership[sample,]
  if (sort) membership <- sort_membership(x$samples$membership[sample,])
  cluster_summary <- table(membership, deparse.level = 0)
  print(cluster_summary, ...)

  cat("\nLog marginal likelihood: ", x$samples$log_mlike[sample], "\n")
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

#' Update method for `sfclust` Objects
#'
#' This function updates the `sfclust` object with the specified clustering sample
#' and prepares the necessary arguments for further analysis.
#'
#' @param x An object of class 'sfclust'.
#' @param sample An integer specifying the clustering sample number to be summarized.
#'        The default is the last sample (i.e., `nrow(x$samples$membership)`).
#' @param ... Additional arguments passed to `print.default`. These are not used directly
#'        in this method but are retained for compatibility.
#' @return The updated `sfclust` object with the modified clustering sample and associated
#'         arguments.
#' @export
update.sfclust <- function(x, sample = nrow(x$samples$membership)) {
  args <- attr(x, "inla_args")
  args$membership <- x$samples$membership[sample,]
  args$stdata <- attr(x, "args")$stdata
  args$stnames <- attr(x, "args")$stnames
  args$correction <- FALSE
  args$detailed <- TRUE

  call <- as.call(c(list(as.name("log_mlik_all")), args))
  x$clustering$id <- sample
  x$clustering$models <- eval(call, envir = parent.frame())
  x
}

#' Fitted Values for `sfclust` Objects
#'
#' This function calculates the fitted values for a specific clustering sample in an
#' `sfclust` object, based on the estimated models for each cluster. The fitted
#' values are computed using the membership assignments and model parameters
#' associated with the selected clustering sample.
#'
#' @param x An object of class 'sfclust', containing clustering results and models.
#' @param sample An integer specifying the clustering sample number for which
#'        the fitted values should be computed. The default is the `id` of the
#'        current clustering. The value must be between 1 and the total number
#'        of clustering (membership) samples.
#' @return A matrix of fitted values, where each row corresponds to a data point
#'         and each column corresponds to a cluster.
#' @details
#' The function first checks if the provided `sample` value is valid (i.e., it is
#' within the range of available clustering samples). If the specified `sample`
#' does not match the current clustering `id`, the `sfclust` object is updated
#' accordingly. It then retrieves the membership assignments and cluster models
#' for the selected sample, calculates the linear predictions for each cluster,
#' and combines them into a matrix of fitted values.
#'
#' @examples
#' # Assuming 'sfclust_obj' is a pre-existing sfclust object
#' fitted_values <- fitted.sfclust(sfclust_obj, sample = 3)
#'
#' @export
fitted.sfclust <- function(x, sample = x$clustering$id, sort = FALSE) {
  if (sample < 1 || sample > nrow(x$samples$membership)) {
    stop("`sample` must be between 1 and the total number clustering (membership) samples.")
  }
  if (sample != x$clustering$id) x <- update(x, sample = sample)

  membership <- x$samples$membership[sample,]
  if (sort) {
    membership <- sort_membership(x$samples$membership[sample,])
    x$clustering$models <- x$clustering$models[attr(membership, "order")]
  }

  clusters <- 1:max(membership)
  predlist <- lapply(
    1:max(membership), linpred_each, membership, x$clustering$models,
    attr(x, "args")$stdata, attr(x, "args")$stnames
  )
  do.call(rbind, predlist)
}

linpred_each <- function(k, membership, models, stdata, stnames){
  df <- data_each(k, membership, stdata, stnames)
  df$linpred <- linpred_each_corrected(models[[k]])
  df$cluster <- k
  df
}

# linear predictor only for terms that are defined at cluster level
linpred_each_corrected <- function(x){
  aux <- if ("id" %in% names(x$summary.random)) x$summary.random[["id"]]$mean else 0
  x$summary.linear.predictor$mean - aux
}

#' Plot function for sfclust objects with a conditional legend and log marginal likelihood
#'
#' This function plots the map with estimated clusters from an sfclust object and adds a legend if the number of clusters is less than 10.
#' It also plots the log marginal likelihood.
#'
#' @param map An sf object provided by the user.
#' @param output An sfclust object containing the clustering results.
#' @param sample The row of the cluster matrix to use for plotting (default is the last row).
#' @param title A title for the plot (default is "Estimated Clusters").
#' @export
plot.sfclust <- function(x, sample = x$clustering$id, which = 1:3, clusters = NULL, sort = FALSE, legend = FALSE, ...) {

  nsamples <- nrow(x$samples$membership)
  if (sample < 1 || sample > nsamples) {
    stop("`sample` must be between 1 and the total number clustering (membership) samples.")
  }

  # get geometries, selected membership and clusters to plot
  geoms <- st_get_dimension_values(attr(x, "args")$stdata, attr(out, "args")$stnames[1])
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
    df <- fitted(out, sample = sample, sort = sort)
    df <- subset(df, cluster %in% clusters)
    plot(df$time, df$linpred, col = df$cluster, main = "Cluster mean functions",
      xlab = "Time", ylab = "Cluster linear predictor", pch = 19)
  }
  if (2 %in% which) { # marginal likelihood convergence
    plot(x$samples$log_mlike, type = "l", main = "Log marginal likelihood convergence",
      xlab = "Sample", ylab = "Log marginal likelihood", col = "blue", lwd = 2)
    points(sample, x$samples$log_mlike[sample], col = "red", pch = 19)
  }
  if (length(which) > 1) par(opar)

  invisible(NULL)
}
