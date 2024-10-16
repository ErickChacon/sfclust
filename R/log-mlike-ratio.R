#' Evaluate Log-Likelihood Ratio for a Proposed Move
#'
#' This function calculates the log-likelihood ratio for a proposed "split" or "merge"
#' move within a clustering or partitioning algorithm using the Integrated Nested
#' Laplace Approximation (INLA) approach. It updates local likelihoods based on
#' whether a split or merge operation has been proposed and supports multiple family
#' distributions.
#'
#' @param move A character string specifying the type of move; either "split" or "merge".
#' @param log_like_vec A numeric vector of existing log likelihoods for all clusters.
#' @param control.move A list containing control parameters and information about the move,
#'        such as `clust_old`, `k`, `cluster_newid`, and `cluster_rm`.
#' @param stars_obj A stars spatial temporal objects including responses, covariates, and other data.
#' @param family The family of distributions to use, e.g., "normal", "poisson", etc.
#' @param formula An object of class \code{\link[stats]{formula}} specifying the model to be used in INLA.
#' @param correction Logical, specifying whether correction calculations should be executed.
#' @param detailed Logical, indicating whether detailed output from INLA is required.
#' @param time_var The name of the time variable in the stars object.
#' @param N_var The name of the variable used as exposure or trials (for Poisson, binomial, etc.).
#'
#' @return A list containing the updated log likelihood vector and the log likelihood ratio (`llratio`).
#' @export
log_mlik_ratio <- function(
    move, log_mlike_vec, proposed_move, stars_obj,
    formula, family = "normal", correction = FALSE, time_var, N_var, detailed = F, ...) {
  # update local marginal likelihoods for split move
  if (move == "split") {
    log_like_vec_new <- log_mlike_vec
    M1 <- log_mlik_each(proposed_move$clust_old, stars_obj, proposed_move$cluster, formula, family, correction, time_var, N_var, detailed, ...)
    M2 <- log_mlik_each(proposed_move$k, stars_obj, proposed_move$cluster, formula, family, correction, time_var, N_var, detailed, ...)
    log_like_vec_new[proposed_move$clust_old] <- M1
    log_like_vec_new[proposed_move$k] <- M2
    llratio <- M1 + M2 - log_mlike_vec[proposed_move$clust_old]
  }

  # update local margina likelihoods for merge move
  if (move == "merge") {
    log_like_vec_new <- log_mlike_vec[-proposed_move$cluster_rm]
    M1 <- log_mlik_each(proposed_move$cluster_newid, stars_obj, proposed_move$cluster, formula, family, correction, time_var, N_var, detailed, ...)
    log_like_vec_new[proposed_move$cluster_newid] <- M1
    llratio <- M1 - sum(log_mlike_vec[c(proposed_move$cluster_rm, proposed_move$cluster_newid)])
  }

  return(list(ratio = llratio, log_mlike_vec = log_like_vec_new))
}
