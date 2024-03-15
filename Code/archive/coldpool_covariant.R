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


