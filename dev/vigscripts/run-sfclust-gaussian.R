library(sfclust)
library(stars)

set.seed(123)
initial_cluster <- genclust(stgaus, nclust = 20)

system.time(
  result <- sfclust(stgaus, graphdata = initial_cluster, logpen = -50,,
    formula = y ~ f(id_time, model = "rw1",
                    hyper = list(prec = list(prior = "normal", param = c(-2, 1)))),
    niter = 50, burnin = 10, thin = 2, nmessage = 10,
    path_save = file.path("dev", "vigdata", "full-gaussian-mcmc1.rds"))
)
#    user  system elapsed
# 331.552  60.037  58.307

system.time(
  result2 <- update(
    result, niter = 1000, nsave = 500,
    path_save = file.path("dev", "vigdata", "full-gaussian-mcmc2.rds")
  )
)
#     user   system  elapsed
# 4685.726  859.919  830.404

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
