library(sfclust)
library(stars)

set.seed(7)
system.time(
  sfclust(stbinom, formula = cases ~ poly(time, 2) + f(id),
    family = "binomial", Ntrials = population,
    niter = 2000, path_save = file.path("dev", "vigdata", "full-binomial-mcmc.rds"))
)
#      user    system   elapsed
# 10781.142  1368.382  2794.899

# Reduce size of object
result <- readRDS(file.path("dev", "vigdata", "full-binomial-mcmc.rds"))
pseudo_inla <- function(x) {
  list(
    summary.random = list(id = x$summary.random$id["mean"]),
    summary.linear.predictor = x$summary.linear.predictor["mean"],
    misc = list(linkfunctions = list(names = "logit"))
  )
}
result$clust$models <- lapply(result$clust$models, pseudo_inla)
saveRDS(result, file.path("inst", "vigdata", "binomial-mcmc.rds"))
cat("binomial-mcmc.rds saved.\n")
