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
  as.integer(clusters_labels[as.character(x)])
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
fitted.sfclust <- function(x, sample = x$clustering$id) {
  if (sample < 1 || sample > nrow(x$samples$membership)) {
    stop("`sample` must be between 1 and the total number clustering (membership) samples.")
  }
  if (sample != x$clustering$id) x <- update(x, sample = sample)

  membership <- x$samples$membership[sample,]
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
    plot(geoms, col = membership,
      main = paste("Clustering sample", sample, "out of", nsamples),
      xlab = "Longitude", ylab = "Latitude", border = "black", ...
    )
    if (legend) {
      legend("topright", legend = levels(membership), fill = 1:length(levels(membership)),
        border = "black", title = "Cluster", cex = 0.8, bty = "n")
    }
  }
  if (3 %in% which) { # functional shapes
    # if (sample != x$clustering$id) x <- update(x, sample = sample)
    df <- fitted(out, sample = sample)
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

#' Plot function for sfclustm objects
#'
#' This function plots the fitted values (mean) from the a-th model in an sfclustm object.
#'
#' @param x An sfclustm object containing the models for each cluster.
#' @param a The index of the model to plot (default is 1).
#' @param title A title for the plot (default is "Fitted Values").
#' @export
plot.sfclustm <- function(x, a = 1, title = "Fitted Values") {

  # Ensure that a is valid (within the range of available models)
  if (a < 1 || a > length(x)) {
    stop("The value of a must be between 1 and the number of models in the sfclustm object.")
  }

  # Extract the fitted values (mean) from the a-th model
  fitted_values <- x[[a]]$summary.fitted.values$mean

  # Ensure that fitted_values exists
  if (is.null(fitted_values)) {
    stop("Fitted values not found in the selected model.")
  }

  # Plot the fitted values
  plot(fitted_values, type = "l", col = "blue", lwd = 2,
       xlab = "Index", ylab = "Fitted Mean",
       main = paste(title, "- Model", a))

  invisible(NULL)  # Return nothing, just plot
}


#' Fit function for sfclust objects with cluster sorting
#'
#' This function applies the log_mlik_each function to fit models for each cluster in sfclust.
#' The clusters are sorted in descending order based on the number of members in each cluster.
#'
#' @param output An sfclust object containing the clustering results.
#' @param k The row of the cluster matrix to use (default is the last row).
#' @param Y The response variable.
#' @param X Covariates or predictor variables.
#' @param N Number of observations or other input (optional).
#' @param formula A model formula.
#' @param family A family object specifying the distribution (default is "normal").
#' @param correction Logical, indicating whether to apply corrections (default is FALSE).
#' @param detailed Logical, indicating whether detailed results should be returned (default is TRUE).
#' @param ... Additional arguments passed to log_mlik_each.
#' @return A list of models, each of class 'sfclustm', fit for each sorted cluster.
#' @export
fit.sfclust <- function(output, k = nrow(output$cluster), Y, X = NULL, N = NULL, formula, family = "normal", correction = FALSE, detailed = TRUE, ...) {

  # Ensure that output is an sfclust object
  if (!inherits(output, "sfclust")) {
    stop("The output must be an object of class 'sfclust'.")
  }

  # Ensure that k is valid (within the range of available rows in the cluster matrix)
  if (k < 1 || k > nrow(output$cluster)) {
    stop("The value of k must be between 1 and the number of rows in the cluster matrix.")
  }

  # Extract the cluster assignments from the k-th row
  cluster <- output$cluster[k, ]

  # Determine the number of clusters (p) based on the max of the cluster vector
  p <- max(cluster)

  # Count the number of members in each cluster
  cluster_sizes <- table(cluster)

  # Sort the clusters by size (descending order)
  sorted_clusters <- as.numeric(names(sort(cluster_sizes, decreasing = TRUE)))

  # Apply log_mlik_each to each cluster (based on the sorted clusters)
  models <- lapply(sorted_clusters, function(c) {
    log_mlik_each(Y, cluster == c, X, N, formula, family, correction = correction, detailed = detailed, ...)
  })

  # Assign the class 'sfclustm' to the result
  class(models) <- "sfclustm"

  return(models)
}

#' Print method for sfclustm objects
#'
#' This function prints the summary of each model within an sfclustm object.
#'
#' @param x An object of class 'sfclustm'.
#' @param ... Additional arguments (currently unused).
#' @export
print.sfclustm <- function(x, ...) {
  cat("sfclustm object containing models for each cluster:\n")
  for (i in seq_along(x)) {
    cat("\nModel for cluster", i, ":\n")
    print(summary(x[[i]]))
  }
}

#' Summary method for sfclustm objects
#'
#' This function provides a summary of each model within an sfclustm object.
#'
#' @param x An object of class 'sfclustm'.
#' @param ... Additional arguments (currently unused).
#' @return A list of summaries for each model.
#' @export
summary.sfclustm <- function(x, ...) {
  # Ensure that x is an sfclustm object
  # if (!inherits(x, "sfclustm")) {
  #   stop("The object must be of class 'sfclustm'.")
  # }

  # Initialize a list to store the summaries
  summaries <- lapply(seq_along(x), function(i) {
    cat("\nSummary of Model for Cluster", i, ":\n")
    summary(x[[i]])  # Call the summary function for each model in the sfclustm object
  })

  invisible(summaries)  # Return the summaries invisibly
}
