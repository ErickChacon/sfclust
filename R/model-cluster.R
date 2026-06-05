#' Bayesian spatial functional clustering from long-format data
#'
#' Core MCMC algorithm for Bayesian detection of neighboring spatial units with similar
#' functional shapes. Works directly on a pre-computed long-format data frame and an
#' adjacency matrix, with no dependency on spatial data structures.
#'
#' For working with `stars` spatio-temporal objects, use the [sfclust()] wrapper, which
#' converts the stars object to long format and calls this function internally.
#'
#' @param data A long-format data frame with at least columns `id` (unique row index)
#'        and `ids` (integer spatial unit index, 1 to `ns`), plus any response and
#'        covariate columns referenced in `formula`.
#' @param adjacency A square pre-weighted adjacency matrix (ns × ns) encoding spatial
#'        contiguity and edge weights. Can be a dense `matrix` or a sparse `Matrix`.
#'        Typically obtained via [igraph::as_adjacency_matrix()] on the graph returned
#'        by [genclust()] or [genclust_adj()].
#' @param graphdata A list with components `graph`, `mst`, and `membership` as returned
#'        by [genclust()] or [genclust_adj()]. If `NULL`, [genclust_adj()] is called on
#'        `adjacency` using `nclust` as the initial number of clusters.
#' @param nclust Integer. Initial number of clusters when `graphdata = NULL`. Ignored if
#'        `graphdata` is provided (default is `10`).
#' @param move_prob A numeric vector of probabilities for different types of moves in the MCMC process:
#'        birth, death, change, and hyperparameter moves (default is `c(0.425, 0.425, 0.1, 0.05)`).
#' @param logpen A negative numeric value representing the log-scale penalty for
#'        increasing the number of clusters by one. The number of clusters is assumed to
#'        follow a geometric prior with probability `q`, making this penalty equal to
#'        `log(1 - q)`. For example, if `logpen = -50`, then a proposal that increases the
#'        number of clusters will only be favored if it improves the log marginal
#'        likelihood by more than 50.
#' @param correction A logical indicating whether correction to compute the marginal
#'        likelihoods should be applied (default is `TRUE`). This depends on the type of
#'        effects included in the `INLA` model.
#' @param niter An integer specifying the number of MCMC iterations after burn-in (default is `100`).
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
#' - The paper: \doi{doi:10.1214/22-AOAS1643}
#' - Supplementary material: \doi{doi:10.1214/22-AOAS1643SUPPB}
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
#' - The paper: <https://doi.org/10.1002/sim.70597>
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
#' \donttest{
#' library(sfclust)
#' library(stars)
#'
#' data(stgaus)
#' df      <- data_all(stgaus)
#' geodata <- genclust(stgaus)
#' adj     <- igraph::as_adjacency_matrix(geodata$graph)
#'
#' result <- sfclust_fit(df, adj, graphdata = geodata,
#'   fun_col = "time", formula = y ~ f(idf_time, model = "rw1"),
#'   niter = 10, nmessage = 1)
#' print(result)
#' summary(result)
#' }
#'
#' @importFrom igraph graph_from_adjacency_matrix mst V E vcount ecount delete_edges components
sfclust_fit <- function(data, adjacency, graphdata = NULL,
                        move_prob = c(0.425, 0.425, 0.1, 0.05), logpen = log(1 - 0.5),
                        nclust = 10,
                        correction = TRUE, niter = 100, burnin = 0, thin = 1,
                        nmessage = 10, path_save = NULL, nsave = nmessage, ...) {

  inla_args <- match.call(expand.dots = FALSE)$...
  # number of regions
  ns <- nrow(adjacency)

  # check if correction is required
  if (correction) {
    if (length(correction_required(eval(inla_args$formula))) == 0) {
      correction <- FALSE
      warning("Log marginal-likelihood correction not required.", immediate. = TRUE)
    }
  }

  # initial clustering
  if (is.null(graphdata)) graphdata <- genclust_adj(adjacency, nclust = nclust)
  graph      <- graphdata[["graph"]]
  mstgraph   <- graphdata[["mst"]]
  membership <- graphdata[["membership"]]

  args <- list(data = data, adjacency = adjacency,
               graphdata = graphdata, move_prob = move_prob, logpen = logpen,
               correction = correction)

  nclust     <- max(membership)
  edge_status <- getEdgeStatus(membership, mstgraph)
  log_mlike_vec <- log_mlik_all(membership, data, correction, detailed = FALSE, ...)
  log_mlike  <- sum(log_mlike_vec)

  # output objects
  niter_total <- niter + burnin
  nsamples    <- floor((niter - 1) / thin) + 1

  membership_out <- array(0, dim = c(nsamples, ns))
  log_mlike_out  <- numeric(nsamples)
  mst_out        <- list()

  birth_cnt <- death_cnt <- change_cnt <- hyper_cnt <- 0

  # MCMC sampling
  for (iter in 1:niter_total) {
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
        rd_new <- 0.9 - rhy
      } else {
        rd_new <- move_prob[2]
      }
      log_P   <- log(rd_new) - log(rb)
      log_A   <- logpen
      log_L_new <- log_mlik_ratio("split", split_res, log_mlike_vec, data, correction, ...)
      acc_prob  <- exp(min(0, log_A + log_P + log_L_new$ratio))

      if (runif(1) < acc_prob) {
        membership  <- split_res$membership
        edge_status <- getEdgeStatus(membership, mstgraph)
        nclust      <- nclust + 1
        log_mlike_vec <- log_L_new$log_mlike_vec
        log_mlike   <- sum(log_mlike_vec)
        birth_cnt   <- birth_cnt + 1
      }
    }

    if (move_choice == 2) { # death move
      merge_res <- mergeCluster(mstgraph, edge_status, membership)

      if (nclust == 2) {
        rb_new <- 1 - rhy
      } else {
        rb_new <- move_prob[1]
      }
      log_P   <- log(rb_new) - log(rd)
      log_A   <- -logpen
      log_L_new <- log_mlik_ratio("merge", merge_res, log_mlike_vec, data, correction, ...)
      acc_prob  <- exp(min(0, log_A + log_P + log_L_new$ratio))

      if (runif(1) < acc_prob) {
        membership  <- merge_res$membership
        edge_status <- getEdgeStatus(membership, mstgraph)
        nclust      <- nclust - 1
        log_mlike_vec <- log_L_new$log_mlike_vec
        log_mlike   <- sum(log_mlike_vec)
        death_cnt   <- death_cnt + 1
      }
    }

    if (move_choice == 3) { # change move
      merge_res <- mergeCluster(mstgraph, edge_status, membership)
      split_res <- splitCluster(mstgraph, nclust - 1, merge_res$membership)

      log_L_new_merge <- log_mlik_ratio("merge", merge_res, log_mlike_vec, data, correction, ...)
      log_L_new       <- log_mlik_ratio("split", split_res, log_L_new_merge$log_mlike_vec, data, correction, ...)
      acc_prob        <- exp(min(0, log_L_new_merge$ratio + log_L_new$ratio))

      if (runif(1) < acc_prob) {
        membership  <- split_res$membership
        edge_status <- getEdgeStatus(membership, mstgraph)
        log_mlike_vec <- log_L_new$log_mlike_vec
        log_mlike   <- sum(log_mlike_vec)
        change_cnt  <- change_cnt + 1
      }
    }

    if (move_choice == 4) { # hyper move
      mstgraph    <- proposeMST(graph, getEdgeStatus(membership, graph))
      V(mstgraph)$vid <- 1:ns
      edge_status <- getEdgeStatus(membership, mstgraph)
      hyper_cnt   <- hyper_cnt + 1
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
      log_mlike_out[isample]    <- log_mlike
      mst_out[[isample]]        <- mstgraph
    }

    # save to file
    if (iter > burnin & (iter - burnin - 1) %% nsave == 0) {
      if (!is.null(path_save)) {
        output <- list(
          samples = list(
            membership  = membership_out,
            log_mlike   = log_mlike_out,
            move_counts = c(births = birth_cnt, deaths = death_cnt, changes = change_cnt, hypers = hyper_cnt)
          ),
          clust = list(
            id       = isample,
            membership = membership,
            models   = NULL
          )
        )
        attr(output, "mst")       <- mst_out
        attr(output, "args")      <- args
        attr(output, "inla_args") <- inla_args
        class(output) <- "sfclust"
        saveRDS(output, file = path_save)
      }
    }
  }

  # final outcome
  membership <- membership_out[nrow(membership_out), ]
  output <- list(
    samples = list(
      membership  = membership_out,
      log_mlike   = log_mlike_out,
      move_counts = c(births = birth_cnt, deaths = death_cnt, changes = change_cnt, hypers = hyper_cnt)
    ),
    clust = list(
      id       = nrow(membership_out),
      membership = membership,
      models   = log_mlik_all(membership, data, correction = FALSE, detailed = TRUE, ...)
    )
  )
  attr(output, "mst")       <- mst_out
  attr(output, "args")      <- args
  attr(output, "inla_args") <- inla_args
  class(output) <- "sfclust"

  if (!is.null(path_save)) saveRDS(output, file = path_save)
  return(output)
}


#' Bayesian spatial functional clustering from stars objects
#'
#' Stars-facing wrapper around [sfclust_fit()]. Converts a `stars` spatio-temporal
#' object to long format, builds the spatial graph, runs the core MCMC algorithm, and
#' attaches stars metadata to the result so that spatial plot methods work.
#'
#' @param stdata A stars object containing response variables, covariates, and other necessary data.
#' @param graphdata A list containing the initial graph used for the Bayesian model.
#'        It should include components like `graph`, `mst`, and `membership` (default is `NULL`).
#' @param sp_dims Character vector with the names of the spatial dimensions of `stdata`.
#'        Use a single name (e.g. `"geometry"`) for vector geometry data, or two names
#'        (e.g. `c("x", "y")`) for raster data. Default is `"geometry"`.
#' @param fun_dims Character vector with the names of the functional dimensions of `stdata`
#'        (e.g. `"time"`). Default is `"time"`.
#' @param move_prob A numeric vector of probabilities for different types of moves in the MCMC process:
#'        birth, death, change, and hyperparameter moves (default is `c(0.425, 0.425, 0.1, 0.05)`).
#' @param logpen A negative numeric value representing the log-scale penalty for
#'        increasing the number of clusters by one. The number of clusters is assumed to
#'        follow a geometric prior with probability `q`, making this penalty equal to
#'        `log(1 - q)`. For example, if `logp = -50`, then a proposal that increases the
#'        number of clusters will only be favored if it improves the log marginal
#'        likelihood by more than 50.
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
#' @return
#' An `sfclust` object (see [sfclust_fit()]) with additional attributes:
#' - `stdata`: the original stars object (used by spatial plot methods).
#' - `sp_dims`: spatial dimension name(s).
#' - `fun_dims`: functional dimension name(s).
#' - `valid_ids`: integer vector of non-NA spatial indices (raster only; `NULL` for vector data).
#'
#' @author
#' Ruiman Zhong \email{ruiman.zhong@kaust.edu.sa},
#' Erick A. Chacón-Montalván \email{erick.chaconmontalvan@kaust.edu.sa},
#' Paula Moraga \email{paula.moraga@kaust.edu.sa}
#'
#' @examples
#'
#' \donttest{
#' library(sfclust)
#'
#' # Clustering with Gaussian data
#' data(stgaus)
#' result <- sfclust(stgaus, formula = y ~ f(idf_time, model = "rw1"),
#'   niter = 10, nmessage = 1)
#' print(result)
#' summary(result)
#' plot(result)
#'
#' # Clustering with binomial data
#' data(stbinom)
#' result <- sfclust(stbinom, formula = cases ~ poly(idf_time, 2) + f(id),
#'   family = "binomial", Ntrials = population, niter = 10, nmessage = 1)
#' print(result)
#' summary(result)
#' plot(result)
#' }
#'
#' @importFrom igraph as_adjacency_matrix
#' @importFrom stars st_get_dimension_values
#' @export
sfclust <- function(stdata, graphdata = NULL, sp_dims = "geometry", fun_dims = "time",
                    move_prob = c(0.425, 0.425, 0.1, 0.05), logpen = log(1 - 0.5),
                    correction = TRUE, niter = 100, burnin = 0, thin = 1,
                    nmessage = 10, path_save = NULL, nsave = nmessage, ...) {

  # capture inla_args here so match.call() sees named arguments directly,
  # not ..1/..2/..3 references that arise when ... is forwarded through sfclust_fit
  inla_args <- match.call(expand.dots = FALSE)$...

  # initial clustering
  if (is.null(graphdata)) {
    if (length(sp_dims) == 1) {
      graphdata <- genclust(st_get_dimension_values(stdata, sp_dims))
    } else {
      graphdata <- genclust(stdata, sp_dims = sp_dims)
    }
  }

  # attach valid_ids to stdata so data_all() filters correctly
  valid_ids <- graphdata[["valid_ids"]]
  if (!is.null(valid_ids)) {
    attr(stdata, "valid_ids") <- valid_ids
  }

  # convert stars to long-format and extract adjacency
  data       <- data_all(stdata, sp_dims, fun_dims)
  adjacency  <- igraph::as_adjacency_matrix(graphdata$graph)

  result <- sfclust_fit(data, adjacency, graphdata,
                        move_prob = move_prob, logpen = logpen,
                        correction = correction, niter = niter,
                        burnin = burnin, thin = thin, nmessage = nmessage,
                        path_save = path_save, nsave = nsave, ...)

  # override inla_args with the correctly captured version from this wrapper call
  attr(result, "inla_args") <- inla_args
  attr(result, "args")$fun_col <- fun_dims[1]

  # attach stars metadata for spatial plot methods and fitted() reconstruction
  attr(result, "stdata")    <- stdata
  attr(result, "sp_dims")   <- sp_dims
  attr(result, "fun_dims")  <- fun_dims
  attr(result, "valid_ids") <- valid_ids

  result
}


log_mlik_ratio <- function(move_type, move, log_mlike_vec, data, correction = TRUE, ...) {
  # update local marginal likelihoods for split move
  if (move_type == "split") {
    log_like_vec_new <- log_mlike_vec
    M1 <- log_mlik_each(move$cluster_old, move$membership, data, correction, detailed = FALSE, ...)
    M2 <- log_mlik_each(move$cluster_new, move$membership, data, correction, detailed = FALSE, ...)
    log_like_vec_new[move$cluster_old] <- M1
    log_like_vec_new[move$cluster_new] <- M2
    llratio <- M1 + M2 - log_mlike_vec[move$cluster_old]
  }

  # update local marginal likelihoods for merge move
  if (move_type == "merge") {
    log_like_vec_new <- log_mlike_vec[- move$cluster_rm]
    M <- log_mlik_each(move$cluster_new, move$membership, data, correction, detailed = FALSE, ...)
    log_like_vec_new[move$cluster_new] <- M
    llratio <- M - sum(log_mlike_vec[c(move$cluster_rm, move$cluster_new)])
  }

  return(list(ratio = llratio, log_mlike_vec = log_like_vec_new))
}
