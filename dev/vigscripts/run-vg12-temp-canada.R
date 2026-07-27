library(sfclust)
library(stars)
library(dplyr)
library(rnaturalearth)
library(here)

data(CanadianWeather, package = "fda")

# prepare data
stations <- as.data.frame(CanadianWeather$coordinates) |>
  mutate(longitud = - `W.longitude`) |>
  rename(latitud = "N.latitude") |>
  select(longitud, latitud) |>
  st_as_sf(coords = c("longitud", "latitud"), crs = st_crs(4326))

time <- seq(as.Date("1977-01-01"), as.Date("1977-12-31"), by = "1 day")
canweather <- st_as_stars(
    temp = t(CanadianWeather$dailyAv[, , 1]),
    ztemp = t(scale(CanadianWeather$dailyAv[, , 1])),
    dimensions = st_dimensions(geometry = st_geometry(stations), time = time, point = TRUE)
)

# create regions domain
stations2 <- st_transform(stations, st_crs(3857))
boundary <- st_convex_hull(st_union(stations2)) |>
  st_buffer(units::set_units(1000, "km"))
domain <- st_cast(st_voronoi(st_union(stations2), boundary)) |>
  st_intersection(boundary)
domain <- domain[as.numeric(st_within(stations2, domain))] |>
  st_transform(st_crs(stations))

# model
set.seed(7)
geodata <- genclust(as(st_touches(domain), "sparseMatrix"), nclust = 35)

set.seed(123)
print(system.time(
result <- sfclust(canweather, graphdata = geodata, logpen = -300,
    formula = ztemp ~ poly(time, 2) + f(id_time, model = "rw1"),
    niter = 4000, burnin = 0, thin = 10, nmessage = 10, nsave = 1000,
    control.inla = list(control.vb = list(enable = FALSE)),
    path_save = here::here(file.path("tools", "data", "canweather-mcmc.rds")))
))
cat("canweather-mcmc.rds saved.\n")
#     user   system  elapsed
# 3295.286 1592.867 3487.686

