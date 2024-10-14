#' Detect Spatial Functional Clusters Based on Bayesian Spanning Tree
#'
#' This function implements a Reversible Jump Markov Chain Monte Carlo (RJMCMC) algorithm
#' for detecting spatial functional clusters based on a Bayesian analysis of spanning trees.
#' It handles various types of moves including birth, death, change, and hyperparameter updates
#' to explore the space of possible cluster configurations.
#'
#' @param data A stars object contains response, covariates, and other needed data
#' @param formula A formula object representing the model to be fitted.
#' @param graph Initial spanning tree used for the Bayesian model.
#' @param init_val List of initial values for parameters 'trees', 'beta', and 'cluster'.
#' @param q q is the penalty of the cluster number.
#' @param niter MCMC Integer, number of MCMC iterations to perform.
#' @param burnin Integer, number of burn-in iterations to discard.
#' @param thin Integer, thinning interval for recording the results.
#' @param path_save Character, the path where results should be saved.
#' @param nsave the number of the gap of saved results in the chain.
#' @param move_prob a vector of probabilities of the birth, death, change steps.
#'
#' @return NULL The function primarily outputs results to a specified path and does not return anything.
#'
#' @export
sfclust <- function(data, formula, graphdata = list(graph = NULL, mst = NULL, cluster = NULL),
                    family = "normal", q = 0.5, correction = FALSE, niter = 100, burnin = 0, thin = 1,
                    path_save = NULL, nsave = 10, move_prob = c(0.425, 0.425, 0.1), ...) {
  ## Setup
  # Dimensions
  ns <- length(unique(data$geometry))

  # Initial values
  graph <- graphdata[["graph"]]
  mstgraph <- graphdata[["mst"]]
  cluster <- graphdata[["cluster"]]
  k <- max(cluster)

  # Hyperparameters
  c <- q

  # Movement counts
  hyper_cnt <- 0
  birth_cnt <- 0
  death_cnt <- 0
  change_cnt <- 0

  ## Initialize
  # Initialize log likelihood vector
  log_mlike_vec <- log_mlik_all(data, cluster, formula, family, correction, ...)
  log_mlike <- sum(log_mlike_vec)

  # Determine if an edge in graph is within a cluster or between two clusters
  edge_status <- getEdgeStatus(cluster, mstgraph)

  ## Prepare output
  cluster_out <- array(0, dim = c((niter - burnin) / thin, ns))
  mst_out <- list()
  log_mlike_out <- numeric((niter - burnin) / thin)

  ## MCMC sampling
  for (iter in 1:niter) {
    rhy <- 0.05
    if (k == 1) {
      rb <- 0.95
      rd <- 0
      rc <- 0
    } else if (k == ns) {
      rb <- 0
      rd <- 0.85
      rc <- 0.1
    } else {
      rb <- move_prob[1]
      rd <- move_prob[2]
      rc <- move_prob[3]
    }

    move_choice <- sample(4, 1, prob = c(rb, rd, rc, rhy))

    if (move_choice == 1) { ## Birth move
      split_res <- splitCluster(mstgraph, k, cluster)
      split_res$k <- k + 1
      membership_new <- split_res$cluster

      if (k == ns - 1) {
        rd_new <- 0.85
      } else {
        rd_new <- move_prob[2]
      }
      log_P <- log(rd_new) - log(rb)
      log_A <- log(1 - c)
      log_L_new <- log_mlik_ratio("split", log_mlike_vec, split_res, data, formula, family, correction, FALSE, ...)
      log_L <- log_L_new$ratio
      acc_prob <- min(0, log_A + log_P + log_L)
      acc_prob <- exp(acc_prob)
      if (runif(1) < acc_prob) {
        cluster <- membership_new
        k <- k + 1
        log_mlike_vec <- log_L_new$log_mlike_vec
        log_mlike <- sum(log_mlike_vec)
        edge_status <- getEdgeStatus(cluster, mstgraph)
        birth_cnt <- birth_cnt + 1
      }
    }

    if (move_choice == 2) { ## Death move
      merge_res <- mergeCluster(mstgraph, edge_status, cluster)
      membership_new <- merge_res$cluster
      cid_rm <- merge_res$cluster_rm

      if (k == 2) {
        rb_new <- 0.85
      } else {
        rb_new <- 0.425
      }
      log_P <- log(rb_new) - log(rd)
      log_A <- -log(1 - c)
      log_L_new <- log_mlik_ratio("merge", log_mlike_vec, merge_res, data, formula, family, correction, FALSE, ...)
      log_L <- log_L_new$ratio
      acc_prob <- min(0, log_A + log_P + log_L)
      acc_prob <- exp(acc_prob)
      if (runif(1) < acc_prob) {
        cluster <- membership_new
        k <- k - 1
        log_mlike_vec <- log_L_new$log_mlike_vec
        log_mlike <- sum(log_mlike_vec)
        edge_status <- getEdgeStatus(cluster, mstgraph)
        death_cnt <- death_cnt + 1
      }
    }

    if (move_choice == 3) { ## Change move
      merge_res <- mergeCluster(mstgraph, edge_status, cluster)
      membership_new <- merge_res$cluster
      cid_rm <- merge_res$cluster_rm
      k <- k - 1

      log_L_new_merge <- log_mlik_ratio("merge", log_mlike_vec, merge_res, data, formula, family, correction, FALSE, ...)

      split_res <- splitCluster(mstgraph, k, merge_res$cluster)
      split_res$k <- k + 1
      membership_new <- split_res$cluster
      k <- k + 1

      log_L_new <- log_mlik_ratio(
        "split", log_L_new_merge$log_mlike_vec, split_res, data,
        formula, family, correction, FALSE, ...
      )
      log_L <- log_L_new$ratio + log_L_new_merge$ratio

      acc_prob <- min(0, log_L)
      acc_prob <- exp(acc_prob)
      if (runif(1) < acc_prob) {
        cluster <- membership_new
        log_mlike_vec <- log_L_new$log_mlike_vec
        log_mlike <- sum(log_mlike_vec)
        edge_status <- getEdgeStatus(cluster, mstgraph)
        change_cnt <- change_cnt + 1
      }
    }

    if (move_choice == 4) { ## Hyper move
      edge_status_G <- getEdgeStatus(cluster, graph)
      mstgraph <- proposeMST(graph, edge_status_G)
      V(mstgraph)$vid <- 1:ns
      edge_status <- getEdgeStatus(cluster, mstgraph)
      hyper_cnt <- hyper_cnt + 1
    }

    ## Store estimates

    if (iter %% 10 == 0) {
      # cat("Iteration ", iter, ": clusters = ", k, ", births = ", birth_cnt, ", deaths = ",
      #   death_cnt, ", changes = ", change_cnt, ", hypers = ", hyper_cnt, ", log_mlike = ", log_mlike, "\n",
      #   sep = ""
      # )
      message("Iteration ", iter, ": clusters = ", k, ", births = ", birth_cnt, ", deaths = ",
        death_cnt, ", changes = ", change_cnt, ", hypers = ", hyper_cnt, ", log_mlike = ", log_mlike, "\n",
        sep = ""
      )
    }

    if (iter > burnin & (iter - burnin) %% thin == 0) {
      mst_out[[(iter - burnin) / thin]] <- mstgraph
      cluster_out[(iter - burnin) / thin, ] <- cluster
      log_mlike_out[(iter - burnin) / thin] <- log_mlike
    }

    if (iter %% nsave == 0) {
      if (!is.null(path_save)) {
        saveRDS(
          list(
            cluster = cluster_out, log_mlike = log_mlike_out, mst = mst_out,
            counts = c(births = birth_cnt, deaths = death_cnt, changes = change_cnt, hypers = hyper_cnt)
          ),
          file = path_save
        )
      }
    }
  }

  p <- max(cluster)
  sorted_cluster <- as.numeric(names(sort(table(cluster), decreasing = TRUE)))
  final_model <- lapply(1:p, log_mlik_each, data, sorted_cluster, formula, family, correction = FALSE, detailed = TRUE, time_var = time_var, N_var = N_var, ...)

  class(final_model) <- "sfclustm"
  # Final result
  output <- list(
    cluster = cluster_out, log_mlike = log_mlike_out, mst = mst_out,
    counts = c(births = birth_cnt, deaths = death_cnt, changes = change_cnt, hypers = hyper_cnt),
    model = final_model,
    formula = formula,
    q = q
  )
  if (!is.null(path_save)) saveRDS(output, file = path_save)
  class(output) <- "sfclust"
  return(output)
}

#' Continue MCMC Clustering Procedure
#'
#' This function continues a Bayesian spatial fusion clustering (sfclust) process based on a previous result.
#' It is specifically designed to extend the MCMC iterations from a stopping point using the prior state
#' of cluster assignments and graphical model parameters.
#'
#' @param result List containing results from a previous sfclust function, expected to include
#'        a minimum spanning tree (`mst`) and cluster assignments.
#' @param data A stars object containing spatial data, including covariates and response variables.
#' @param membership Integer vector, indicating the cluster membership of each observation.
#' @param formula An object of class \code{\link[stats]{formula}}, specifying the model used in the analysis.
#' @param family Character string specifying the family of distributions to use for the model.
#'        Defaults to "normal".
#' @param q q is the penalty of the number of cluster.
#' @param correction Logical indicating whether correction for overdispersion or other factors
#'        should be applied. Defaults to `FALSE`.
#' @param niter Integer specifying the number of additional MCMC iterations to perform.
#' @param burnin Integer specifying the number of burn-in iterations to discard in this continuation.
#' @param thin Integer specifying the thinning interval for recording the results.
#' @param path_save Character string specifying the file path to save the continuation results.
#'        If `NULL`, results may not be saved to file.
#' @param time_var The name of the time variable in the stars object.
#' @param N_var The name of the variable used as exposure or trials (for Poisson, binomial, etc.).
#' @param move_prob Numeric vector specifying the probabilities for different move types in the MCMC process.
#' @param ... Additional arguments passed to the underlying `sfclust` function.
#'
#' @details The function takes the last state of the Markov chain from a previous `sfclust` execution and
#'          uses it as the starting point for additional MCMC iterations. This is useful for extending
#'          the analysis without restarting the process, thereby saving computational resources and time.
#'
#' @return The function does not return a value within R but may output results to files if `path_save`
#'         is specified. The results include updated cluster assignments and model parameters after
#'         the additional MCMC iterations.
#'
#' @export
continue_sfclust <- function(result, data, membership,
                          formula, family = "normal", q = 0.5,
                          correction = FALSE, niter = 100, burnin = 0, thin = 1, path_save = NULL, nsave = 10,
                          time_var, N_var, move_prob = c(0.425, 0.425, 0.1), ...) {
  n <- length(result[["mst"]])
  cluster <- result[['cluster']][n,]
  mst <- result[["mst"]][[n]]
  c = q
  sfclust(data = data,
          graphdata = list(graph = graph, mst = mst, cluster = cluster),
          formula = formula, family = family, q = c,
          correction = correction, niter = niter, burnin = burnin, thin = thin,
          path_save = path_save, nsave = nsave, time_var = time_var, N_var = N_var, move_prob = move_prob, ...
  )
}
