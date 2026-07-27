library(sfclust)
library(stars)

data(chapa)

set.seed(7)
print(system.time(
result <- sfclust(chapa, nclust = 335, formula = ndwi2 ~ poly(time, 2),
    niter = 4000, thin = 10, nmessage = 100, nsave = 1000,
    path_save = file.path("dev", "vigdata", "full-chapa-mcmc.rds")
)
))
# Iteration 4000: clusters = 46, births = 40, deaths = 329, changes = 43, hypers = 212, log_mlike = 46972.2812633967
#
#     user   system  elapsed
# 2302.479 1555.784 3219.032

# Reduce size of object
result <- readRDS(file.path("dev", "vigdata", "full-chapa-mcmc.rds"))
pseudo_inla <- function(x) {
  list(
    summary.linear.predictor = x$summary.linear.predictor["mean"],
    misc = list(linkfunctions = list(names = "identity"))
  )
}
result$clust$models <- lapply(result$clust$models, pseudo_inla)
saveRDS(result, file.path("tools", "data", "chapa-mcmc.rds"))
cat("chapa-mcmc.rds saved.\n")
