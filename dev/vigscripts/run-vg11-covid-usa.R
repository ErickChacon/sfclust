library(sfclust)
library(stars)
library(here)

usacovid <- readRDS(here::here(file.path("tools", "data", "usacovid.rds")))


set.seed(123)
formula <- cases ~ 1 + poly(time, 3) + f(id_time, model = "ar1") + f(id)
result <- sfclust(usacovid, nclust = 49, spnames = "space",
    formula = formula, family = "poisson", E = expected,
    niter = 4000, burnin = 0, thin = 10, nmessage = 10, nsave = 1000,
    control.inla = list(control.vb = list(enable = FALSE), int.strategy = "eb"),
    path_save = here::here(file.path("tools", "data", "usacovid-mcmc.rds")))
cat("usacovid-mcmc.rds saved.\n")
