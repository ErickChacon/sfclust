library(sfclust)
library(stars)

set.seed(123)
formula <- y ~ f(id_time, model = "rw1",
  hyper = list(prec = list(prior = "normal", param = c(-2, 1))))

system.time(
  result <- sfclust(stgaus, 20, formula = formula, logpen = -50,
    niter = 50, burnin = 10, thin = 2, nmessage = 10,
    path_save = file.path("dev", "vigdata", "full-gaussian-mcmc1.rds"))
)
#    user  system elapsed
#  58.105  37.554  73.656

system.time(
  result2 <- update(
    result, niter = 1000, nsave = 500,
    path_save = file.path("dev", "vigdata", "full-gaussian-mcmc2.rds"))
)
#     user   system  elapsed
# 712.609  463.280  903.903

# Reduce size of objects
result1 <- readRDS(file.path("dev", "vigdata", "full-gaussian-mcmc1.rds"))
result2 <- readRDS(file.path("dev", "vigdata", "full-gaussian-mcmc2.rds"))
pseudo_inla <- function(x) {
  list(
    summary.linear.predictor = x$summary.linear.predictor["mean"],
    misc = list(linkfunctions = list(names = "identity"))
  )
}
result1$clust$models <- NULL
result2$clust$models <- lapply(result2$clust$models, pseudo_inla)
saveRDS(result1, file.path("inst", "vigdata", "gaussian-mcmc1.rds"))
saveRDS(result2, file.path("inst", "vigdata", "gaussian-mcmc2.rds"))
cat("gaussian-mcmc1.rds and gaussian-mcmc2.rds saved.\n")
