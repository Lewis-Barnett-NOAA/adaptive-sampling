library(sdmTMB)
library(ggplot2)
library(INLA)
library(mgcv)
library(sdmTMB)
library(dplyr)
library(gratia)
library(raster)
library(sp)

#read data
m <- readRDS("Data/ice_coldpool_lm.RDS")
m

set.seed(3)

#sim from cold pool data 10 sims
sims <- simulate(m, nsim = 10)
sims

#observed sea ice measurements
m$model$march_sea_ice

#plot compare observed sea ice to sim 1 cold pool
plot(m$model$march_sea_ice, sims$sim_1)

#bind observed sea ice and cold pool sims
newdata <- cbind(m$model$march_sea_ice, sims)


#Simulations with added norm dist noise, 10 sims
for(i in 1:10) {
  tmp <- predict(m) + rnorm(length(predict(m)), 0, 0.5*summary(m)$sigma)
  res <- if (i==1) tmp else cbind(res, tmp)
}
head(res)

#Simulate from 0 to 1 sea ice values
simulateX <- function(object, nsim = 1, seed = NULL, X, ...) {
  object$fitted.values <- predict(object, X)
  simulate(object = object, nsim = nsim, seed = seed, ...)
}
x <- data.frame(march_sea_ice = seq(0, 1, by = 0.05))
x
## Try it out
order_sims <- simulateX(m,  nsim = 4, X =x)
ordersimsx <- cbind(order_sims, x)
ordersimsx



predictor_dat <- expand.grid(
  X = seq(0, 1, length.out = 100),
  Y = seq(0, 1, length.out = 100)
)

gradient_raster <- poly_gradient(predictor_dat, 180)

#Operating model 1 constant state
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = poisson(link = "log"),
  range = 1.4,
  sigma_O = 0.2,
  seed = 1,
  B = 1  # B0 = intercept  # Add strata to the operating model
)

ggplot(sim_dat, aes(X, Y)) +
  geom_raster(aes(fill = exp(eta))) + # mean without observation error
  geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)

################

Polygon1 <- Polygon(rbind(c(4498482, 2668272), c(4498482, 2669343),
                          c(4499991, 2669343), c(4499991, 2668272)))
Polygon1 <- Polygons(list(Polygon1), 1)
Polygon1 <- SpatialPolygons(list(Polygon1))

## Funcs ##############
inv_atan2 <- function(angle) {
  ## deg to rad
  radian = angle * pi / 180 
  ## invert atan2
  x1 = cos(radian)
  y1 = sin(radian)
  return(cbind(x1, y1))
}
poly_gradient <- function(shape, angle) {
  a = inv_atan2(angle)[1]
  b = inv_atan2(angle)[2]
  x = seq(-1, 1, len = 100)
  y = seq(-1, 1, len = 100)
  ras = raster(outer(x, y, function(x,y) {a * x + b * y}))
  extent(ras) <- extent(shape)
  projection(ras) <- proj4string(shape)
  invisible(ras)
}

## Tests ##############
r12 <- poly_gradient(Polygon1, angle = 180)
plot(Polygon1)
plot(r12, add = T)

