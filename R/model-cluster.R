#' Detect Spatial Functional Clusters Based on Bayesian Spanning Tree
#'
#' This function implements a Reversible Jump Markov Chain Monte Carlo (RJMCMC) algorithm
#' for detecting spatial functional clusters based on a Bayesian analysis of spanning trees.
#' It handles various types of moves including birth, death, change, and hyperparameter updates
#' to explore the space of possible cluster configurations.
#'
#' @param data A stars object contains response, covariates, and other needed data
#' @param graph Initial spanning tree used for the Bayesian model, list(graph = NULL, mst = NULL, cluster = NULL).
#' @param init_val List of initial values for parameters 'trees', 'beta', and 'cluster'.
#' @param q q is the penalty of the cluster number.
#' @param niter MCMC Integer, number of MCMC iterations to perform.
#' @param burnin Integer, number of burn-in iterations to discard.
#' @param thin Integer, thinning interval for recording the results.
#' @param path_save Character, the path where results should be saved.
#' @param nsave the number of the gap of saved results in the chain.
#' @param time_var the variable name of the time dimension in the data.
#' @param N_var the variable name of the N dimension in the data when the it is necessary.
#' @param move_prob a vector of probabilities for birth, death, change and hyper moves
#' @param formula A formula object representing the model to be fitted.
#' respectively.
#'
#' @return NULL The function primarily outputs results to a specified path and does not return anything.
#'
#' @export
sfclust <- function(stdata, graphdata = NULL, stnames = c("geometry", "time"),
                    move_prob = c(0.425, 0.425, 0.1, 0.05), q = 0.5, correction = TRUE,
                    niter = 100, burnin = 0, thin = 1, nmessage = 10, path_save = NULL, nsave = nmessage, ...) {

  # args <- list(stdata = stdata, stnames = stnames, move_prob = move_prob, q = q, correction = correction)
  # inla_args <- rlang::enquos(...)
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
