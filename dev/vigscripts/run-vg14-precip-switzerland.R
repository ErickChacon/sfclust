library(sfclust)
library(stars)
library(pbs)
library(here)

# load worldclim data (prepared in scripts/30-process/44-switzerland-precip.Rmd)
cheprec <- readRDS(here::here(file.path("tools", "data", "cheprec.rds")))
cat("n valid cells:", sum(!is.na(cheprec[["prec"]][, , 1])), "\n")

# model
set.seed(7)
formula <- prec ~ pbs(month, df = 5, Boundary.knots = c(0.5, 12.5))
cat("formula:", deparse(formula), "\n")

print(system.time(
result <- sfclust(cheprec, nclust = 823, spnames = c("x", "y"), formula = formula,
    niter = 4000, thin = 10, nmessage = 100, nsave = 1000,
    path_save = here::here("dev", "vigdata", "full-cheprec-mcmc.rds")
)
))

# reduce size of object
result <- readRDS(here::here("dev", "vigdata", "full-cheprec-mcmc.rds"))
pseudo_inla <- function(x) {
  list(
    summary.linear.predictor = x$summary.linear.predictor["mean"],
    misc = list(linkfunctions = list(names = "identity"))
  )
}
result$clust$models <- lapply(result$clust$models, pseudo_inla)
saveRDS(result, here::here("tools", "data", "cheprec-mcmc.rds"))
cat("cheprec-mcmc.rds saved.\n")
