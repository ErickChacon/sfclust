#' Bayesian spatial functional clustering
#'
#' Bayesian detection of neighboring spatial regions with similar functional shapes using
#' spanning trees and latent Gaussian models. It ensures spatial contiguity in the
#' clusters, handles a large family of latent Gaussian models supported by `inla`, and
#' allows to work with non-Gaussian likelihoods.
#'
#' @param stdata A stars object containing response variables, covariates, and other necessary data.
#' @param graphdata A list containing the initial graph used for the Bayesian model.
#'        It should include components like `graph`, `mst`, and `membership` (default is `NULL`).
#' @param stnames A character vector specifying the spatio-temporal dimension names of
#'        `stdata` that represent spatial geometry and time, respectively (default is
#'        `c("geometry", "time")`).
#' @param move_prob A numeric vector of probabilities for different types of moves in the MCMC process:
#'        birth, death, change, and hyperparameter moves (default is `c(0.425, 0.425, 0.1, 0.05)`).
#' @param q A numeric value representing the penalty for the number of clusters (default is `0.5`).
#' @param correction A logical indicating whether correction to compute the marginal
#'        likelihoods should be applied (default is `TRUE`). This depend of the type of
#'        effect inclused in the `INLA` model.
#' @param niter An integer specifying the number of MCMC iterations to perform (default is `100`).
#' @param burnin An integer specifying the number of burn-in iterations to discard (default is `0`).
#' @param thin An integer specifying the thinning interval for recording the results (default is `1`).
#' @param nmessage An integer specifying how often progress messages should be printed (default is `10`).
#' @param path_save A character string specifying the file path to save the results (default is `NULL`).
#' @param nsave An integer specifying the number of iterations between saved results in the chain
#'        (default is `nmessage`).
#' @param ... Additional arguments such as `formula`, `family`, and others that are passed
#'        to the `inla` function.
#'
#' @details
#' This implementation draws inspiration from the methods described in the paper:
#' *"Bayesian Clustering of Spatial Functional Data with Application to a Human Mobility
#' Study During COVID-19"* by Bohai Zhang, Huiyan Sang, Zhao Tang Luo, and Hui Huang,
#' published in *The Annals of Applied Statistics*, 2023. For further details on the
#' methodology, please refer to:
#' - The paper: <https://doi.org/10.1214/22-AOAS1643>
#' - Supplementary material: <https://doi.org/10.1214/22-AOAS1643SUPPB>
#'
#' The MCMC algorithm in this implementation is largely based on the supplementary
#' material provided in the paper. However, we have generalized the computation of the
#' marginal likelihood ratio by leveraging INLA (Integrated Nested Laplace Approximation).
#' This generalization enables integration over all parameters and hyperparameters,
#' allowing for inference within a broader family of distribution functions and model
#' terms, thereby extending the scope and flexibility of the original approach.
#' Further details of our approach can be found in our paper *"Bayesian spatial functional
#' data clustering: applications in disease surveillance"* by Ruiman Zhong, Erick A.
#' Chacón-Montalván, Paula Moraga:
#' - The paper: <https://arxiv.org/abs/2407.12633>
#'
#' @return
#' An `sfclust` object containing two main lists: `samples` and `clust`.
#' - The `samples` list includes details from the sampling process, such as:
#'   - `membership`: The cluster membership assignments for each sample.
#'   - `log_marginal_likelihood`: The log marginal likelihood for each sample.
#'   - `move_counts`: The counts of each type of move during the MCMC process.
#' - The `clust` list contains information about the selected clustering, including:
#'   - `id`: The identifier of the selected sample (default is the last sample).
#'   - `membership`: The cluster assignments for the selected sample.
#'   - `models`: The fitted models for each cluster in the selected sample.
#'
#' @author
#' Ruiman Zhong \email{ruiman.zhong@kaust.edu.sa},
#' Erick A. Chacón-Montalván \email{erick.chaconmontalvan@kaust.edu.sa},
#' Paula Moraga \email{paula.moraga@kaust.edu.sa}
#'
#' @examples
#'
#' \dontrun{
#'
#' library(stars)
#'
#' # Gaussian model
#' out <- sfclust(stdata, formula = cases ~ temperature, niter = 100)
#' print(out)
#' summary(out)
#' plot(out)
#'
#'
#' # Poisson model
#' out <- sfclust(stdata, formula = cases ~ temperature + f(id, model = 'iid'),
#'     family = "poisson", E = expected, niter = 100)
#' print(out)
#' summary(out)
#' plot(out)
#'
#' }
#'
#' @export
sfclust <- function(stdata, graphdata = NULL, stnames = c("geometry", "time"),
                    move_prob = c(0.425, 0.425, 0.1, 0.05), q = 0.5, correction = TRUE,
                    niter = 100, burnin = 0, thin = 1, nmessage = 10, path_save = NULL, nsave = nmessage, ...) {

  inla_args <- match.call(expand.dots = FALSE)$...

  # number of regions
  geoms <- st_get_dimension_values(stdata, stnames[1])
  ns <- length(geoms)

  # check if correction is required
  if (correction) {
    if (length(correction_required(eval(inla_args$formula))) == 0) {
      correction <- FALSE
      warning("Log marginal-likelihood correction not required.", immediate. = TRUE)
    }
  }

  # initial clustering
  if (is.null(graphdata)) graphdata <- genclust(geoms)
  graph <- graphdata[["graph"]]
  mstgraph <- graphdata[["mst"]]
  membership <- graphdata[["membership"]]

  args <- list(stdata = stdata, graphdata = graphdata, stnames = stnames,
    move_prob = move_prob, q = q, correction = correction)

  nclust <- max(membership)
  edge_status <- getEdgeStatus(membership, mstgraph) # edge is within or between clusters
  log_mlike_vec <- log_mlik_all(membership, stdata, stnames, correction, detailed = FALSE, ...)
  log_mlike <- sum(log_mlike_vec)

  # output objects
  niter <- floor((niter - burnin - 1) / thin) * thin + 1 + burnin
  nsamples <- (niter - burnin - 1) / thin + 1

  membership_out <- array(0, dim = c(nsamples, ns))
  log_mlike_out <- numeric(nsamples)
  mst_out <- list()

  birth_cnt <- death_cnt <- change_cnt <- hyper_cnt <- 0

  # MCMC sampling
  for (iter in 1:niter) {
    rhy <- move_prob[4]
    if (nclust == 1) {
      rb <- 1 - rhy
      rd <- 0
      rc <- 0
    } else if (nclust == ns) {
      rb <- 0
      rd <- 0.9 - rhy
      rc <- 0.1
    } else {
      rb <- move_prob[1]
      rd <- move_prob[2]
      rc <- move_prob[3]
    }

    move_choice <- sample(4, 1, prob = c(rb, rd, rc, rhy))

    if (move_choice == 1) { # birth move
      split_res <- splitCluster(mstgraph, nclust, membership)

      if (nclust == ns - 1) {
        rd_new <-  0.9 - rhy
      } else {
        rd_new <- move_prob[2]
      }
      log_P <- log(rd_new) - log(rb)
      log_A <- log(1 - q)
      log_L_new <- log_mlik_ratio("split", split_res, log_mlike_vec, stdata, stnames, correction, ...)
      acc_prob <- exp(min(0, log_A + log_P + log_L_new$ratio))

      if (runif(1) < acc_prob) {
        membership <- split_res$membership
        edge_status <- getEdgeStatus(membership, mstgraph)
        nclust <- nclust + 1
        log_mlike_vec <- log_L_new$log_mlike_vec
        log_mlike <- sum(log_mlike_vec)
        birth_cnt <- birth_cnt + 1
      }
    }

    if (move_choice == 2) { # death move
      merge_res <- mergeCluster(mstgraph, edge_status, membership)

      if (nclust == 2) {
        rb_new <- 1 - rhy
      } else {
        rb_new <- move_prob[1]
      }
      log_P <- log(rb_new) - log(rd)
      log_A <- -log(1 - q)
      log_L_new <- log_mlik_ratio("merge", merge_res, log_mlike_vec, stdata, stnames, correction, ...)
      acc_prob <- exp(min(0, log_A + log_P + log_L_new$ratio))

      if (runif(1) < acc_prob) {
        membership <- merge_res$membership
        edge_status <- getEdgeStatus(membership, mstgraph)
        nclust <- nclust - 1
        log_mlike_vec <- log_L_new$log_mlike_vec
        log_mlike <- sum(log_mlike_vec)
        death_cnt <- death_cnt + 1
      }
    }

    if (move_choice == 3) { # change move
      merge_res <- mergeCluster(mstgraph, edge_status, membership)
      split_res <- splitCluster(mstgraph, nclust - 1, merge_res$membership)

      log_L_new_merge <- log_mlik_ratio("merge", merge_res, log_mlike_vec, stdata, stnames, correction, ...)
      log_L_new <- log_mlik_ratio("split", split_res, log_L_new_merge$log_mlike_vec, stdata, stnames, correction, ...)
      acc_prob <- exp(min(0, log_L_new_merge$ratio + log_L_new$ratio))

      if (runif(1) < acc_prob) {
        membership <- split_res$membership
        edge_status <- getEdgeStatus(membership, mstgraph)
        log_mlike_vec <- log_L_new$log_mlike_vec
        log_mlike <- sum(log_mlike_vec)
        change_cnt <- change_cnt + 1
      }
    }

    if (move_choice == 4) { ## hyper move
      mstgraph <- proposeMST(graph, getEdgeStatus(membership, graph))
      V(mstgraph)$vid <- 1:ns
      edge_status <- getEdgeStatus(membership, mstgraph)
      hyper_cnt <- hyper_cnt + 1
    }

    # report status
    if (iter %% nmessage == 0) {
      message("Iteration ", iter, ": clusters = ", nclust, ", births = ", birth_cnt, ", deaths = ",
        death_cnt, ", changes = ", change_cnt, ", hypers = ", hyper_cnt, ", log_mlike = ", log_mlike, "\n",
        sep = ""
      )
    }

    # store estimates
    if (iter > burnin & (iter - burnin - 1) %% thin == 0) {
      isample <- (iter - burnin - 1) / thin + 1
      membership_out[isample, ] <- membership
      log_mlike_out[isample] <- log_mlike
      mst_out[[isample]] <- mstgraph
    }

    # save to file
    if (iter > burnin & (iter - burnin - 1) %% nsave == 0) {
      if (!is.null(path_save)) {
        output <- list(
           samples = list(
             membership = membership_out,
             log_mlike = log_mlike_out,
             move_counts = c(births = birth_cnt, deaths = death_cnt, changes = change_cnt, hypers = hyper_cnt)
           ),
           clust = list(
             id = isample,
             membership = membership,
             models = NULL
           )
         )
         attr(output, "mst") <- mst_out
         attr(output, "args") <- args
         attr(output, "inla_args") <- inla_args
         class(output) <- "sfclust"

        saveRDS(output, file = path_save)
      }
    }
  }

  # final outcome
  membership <- membership_out[nrow(membership_out),]
  output <- list(
    samples = list(
      membership = membership_out,
      log_mlike = log_mlike_out,
      move_counts = c(births = birth_cnt, deaths = death_cnt, changes = change_cnt, hypers = hyper_cnt)
    ),
    clust = list(
      id = nrow(membership_out),
      membership = membership,
      models = log_mlik_all(membership, stdata, stnames, correction = FALSE, detailed = TRUE, ...)
    )
  )
  attr(output, "mst") <- mst_out
  attr(output, "args") <- args
  attr(output, "inla_args") <- inla_args
  class(output) <- "sfclust"

  if (!is.null(path_save)) saveRDS(output, file = path_save)
  return(output)
}

log_mlik_ratio <- function(move_type, move, log_mlike_vec, stdata, stnames = c("geometry", "time"),
                           correction = TRUE, ...) {
  # update local marginal likelihoods for split move
  if (move_type == "split") {
    log_like_vec_new <- log_mlike_vec
    M1 <- log_mlik_each(move$cluster_old, move$membership, stdata, stnames, correction, detailed = FALSE, ...)
    M2 <- log_mlik_each(move$cluster_new, move$membership, stdata, stnames, correction, detailed = FALSE, ...)
    log_like_vec_new[move$cluster_old] <- M1
    log_like_vec_new[move$cluster_new] <- M2
    llratio <- M1 + M2 - log_mlike_vec[move$cluster_old]
  }

  # update local marginal likelihoods for merge move
  if (move_type == "merge") {
    log_like_vec_new <- log_mlike_vec[- move$cluster_rm]
    M <- log_mlik_each(move$cluster_new, move$membership, stdata, stnames, correction, detailed = FALSE, ...)
    log_like_vec_new[move$cluster_new] <- M
    llratio <- M - sum(log_mlike_vec[c(move$cluster_rm, move$cluster_new)])
  }

  return(list(ratio = llratio, log_mlike_vec = log_like_vec_new))
}
