library(sdmTMB)
library(ggplot2)

set.seed(123)

#Sampling units
n <- 1000

#make data frame including fake "depth" gradient
predictor_dat <- data.frame(
  expand.grid(X = 1:(n/10), Y = 1:(n/10)),
  depth = rep(rev(seq(1, n/10)), each = n/10)
)

ggplot(predictor_dat, aes(X, Y)) +
  geom_tile(aes(fill = depth)) +
  scale_color_gradient2()

#get triangulated mesh to simulate from
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 200)
plot(mesh)

#define range to use for simulation (higher = smoother)
ranges <- seq(10, 100, 10) 
i <- 5 # choosing a specific range value for example (but can do all in a loop)

# ## Depth coefficient slope = 0.9 (strong gradient)
# sim_dat <- sdmTMB_simulate(
#   formula = ~ 1 + depth,
#   data = predictor_dat,
#   mesh = mesh,
#   family = gaussian(link = "identity"),
#   B = c(0.1, 0.9), # B0 = intercept, B1 = depth coefficient slope
#   range = ranges[i],
#   sigma_O = 0.2,
#   phi = 5, # SD of observation error in this (Gaussian) case
#   seed = 42
# )
# 
# # Shift values so that minimum abundance is 0 (or fit with a different family)
# sim_dat$eta <- (sim_dat$eta - min(sim_dat$eta))
# sim_dat$observed <- (sim_dat$observed - min(sim_dat$observed))
# 
# # Plot simulated "true" data
# ggplot(sim_dat, aes(X, Y)) +
#   geom_tile(aes(fill = eta)) +
#   scale_color_gradient2()
# 
# 
# # Plot simulated "observed" data
# ggplot(sim_dat, aes(X, Y)) +
#   geom_tile(aes(fill = observed)) +
#   scale_color_gradient2()
# 
# 
# 
# #Depth coefficient slope = 0.6 (mid gradient)
# sim_dat2 <- sdmTMB_simulate(
#   formula = ~ 1 + depth,
#   data = predictor_dat,
#   mesh = mesh,
#   family = gaussian(link = "identity"),
#   B = c(0.1, 0.6), # B0 = intercept, B1 = depth coefficient slope
#   range = ranges[i],
#   sigma_O = 0.2,
#   phi = 5, # SD of observation error in this (Gaussian) case
#   seed = 42
# )
# 
# sim_dat2$eta <- (sim_dat2$eta - min(sim_dat2$eta))
# sim_dat2$observed <- (sim_dat2$observed - min(sim_dat2$observed))
# 
# # Plot simulated "observed" data
# ggplot(sim_dat2, aes(X, Y)) +
#   geom_tile(aes(fill = observed)) +
#   scale_color_gradient2()
# 
# 
# 
# 
# ##depth coefficient slope = 0.2 (weak gradient)
# sim_dat3 <- sdmTMB_simulate(
#   formula = ~ 1 + depth,
#   data = predictor_dat,
#   mesh = mesh,
#   family = gaussian(link = "identity"),
#   B = c(0.1, 0.2), # B0 = intercept, B1 = depth coefficient slope
#   range = ranges[i],
#   sigma_O = 0.2,
#   phi = 5, # SD of observation error in this (Gaussian) case
#   seed = 42
# )
# 
# sim_dat3$eta <- (sim_dat3$eta - min(sim_dat3$eta))
# sim_dat3$observed <- (sim_dat3$observed - min(sim_dat3$observed))
# 
# 
# # Plot simulated "observed" data
# ggplot(sim_dat3, aes(X, Y)) +
#   geom_tile(aes(fill = observed)) +
#   scale_color_gradient2()



#COLD POOL SIM
m <- readRDS("Data/ice_coldpool_lm.RDS")

simulateX <- function(object, nsim = 1, seed = NULL, X, ...) {
  object$fitted.values <- predict(object, X)
  simulate(object = object, nsim = nsim, seed = seed, ...)
}
x <- data.frame(march_sea_ice = seq(0, 1, by = 0.05))
x
## Try it out
order_sims <- simulateX(m,  nsim = 4, X =x)
ordersimsx <- cbind(order_sims, x)
ordersimsx[ordersimsx < 0] <- 0 # set negative cold pool extents to 0
ordersimsx

# Draw a random march sea ice proportion
ice_value <- sample(ordersimsx$march_sea_ice,1)

# Draw a cold pool extent simulated from a random march sea ice proportion
cold_pool_value <- as.numeric(ordersimsx[ordersimsx$march_sea_ice == ice_value,][sample(1:(ncol(ordersimsx)-1), 1, replace = TRUE)])
cold_pool_value

# Define the ranges for high, mid, and low sea ice values to guide sample allocation
ice_cat <- quantile(m$model$march_sea_ice, probs = c(0, .5, 1))
low_ice <- c(0, ice_cat[2])
mid_ice <- c(ice_cat[2], ice_cat[3])
high_ice <- c(ice_cat[3], 1)

# Define the ranges for high, mid, and low cold pool extent to determine the operating model
cp_cat <- quantile(m$model$area_lte2_km2, probs = c(0, .333, .666, 1))
low_cp <- c(cp_cat[1], cp_cat[2])
mid_cp <- c(cp_cat[2], cp_cat[3])
high_cp <- c(cp_cat[3], Inf)

# Function to determine the operating model based on sea ice value
get_operating_model <- function(cold_pool_value) {
  if (cold_pool_value >= high_cp[1] && cold_pool_value <= high_cp[2]) {
    sim_dat <- sdmTMB_simulate(
      formula = ~ 1 + depth,
      data = predictor_dat,
      mesh = mesh,
      family = gaussian(link = "identity"),
      B = c(0.1, 0.9), # B0 = intercept, B1 = depth coefficient slope
      range = ranges[i],
      sigma_O = 0.2,
      phi = 5, # SD of observation error in this (Gaussian) case
      seed = #42
    )
    
    sim_dat$eta <- (sim_dat$eta - min(sim_dat$eta))
    sim_dat$observed <- (sim_dat$observed - min(sim_dat$observed))
   
    return(sim_dat)
     
  } else if (cold_pool_value >= mid_cp[1] && cold_pool_value <= mid_cp[2]) {
    sim_dat <- sdmTMB_simulate(
      formula = ~ 1 + depth,
      data = predictor_dat,
      mesh = mesh,
      family = gaussian(link = "identity"),
      B = c(0.1, 0.6), # B0 = intercept, B1 = depth coefficient slope
      range = ranges[i],
      sigma_O = 0.2,
      phi = 5, # SD of observation error in this (Gaussian) case
      seed = #42
    )
    
    sim_dat$eta <- (sim_dat$eta - min(sim_dat$eta))
    sim_dat$observed <- (sim_dat$observed - min(sim_dat$observed))
    
    return(sim_dat)
    
  } else if (cold_pool_value >= low_cp[1] && cold_pool_value <= low_cp[2]) {
    sim_dat <- sdmTMB_simulate(
      formula = ~ 1 + depth,
      data = predictor_dat,
      mesh = mesh,
      family = gaussian(link = "identity"),
      B = c(0.1, 0.2), # B0 = intercept, B1 = depth coefficient slope
      range = ranges[i],
      sigma_O = 0.2,
      phi = 5, # SD of observation error in this (Gaussian) case
      seed = #42
    )
    
    sim_dat$eta <- (sim_dat$eta - min(sim_dat$eta))
    sim_dat$observed <- (sim_dat$observed - min(sim_dat$observed))
    
    return(sim_dat)
    
  } else {
    return("No Matching Operating Model")
  }
}

d <- get_operating_model(cold_pool_value)
d

# Simulate sampling
abundance <- function(d, ice_value) {
  
  d$strata <- ifelse(d$Y >= 50, 1, 2) # Define strata (not evenly split somehow?)
  
  if (ice_value >= high_ice[1] && ice_value <= high_ice[2]) {
    samples_north <- sample(d[d$strata == 1, "observed"], 20)
    samples_south <- sample(d[d$strata == 2, "observed"], 100)

  } else if (ice_value >= mid_ice[1] && ice_value <= mid_ice[2]) {
    samples_north <- sample(d[d$strata == 1, "observed"], 40)
    samples_south <- sample(d[d$strata == 2, "observed"], 80)
    
  } else if (ice_value >= low_ice[1] && ice_value <= low_ice[2]) {
    samples_north <- sample(d[d$strata == 1, "observed"], 60)
    samples_south <- sample(d[d$strata == 2, "observed"], 60)
  }
  
  # estimate abundance
  est <- sum(mean(samples_north) * nrow(d[d$strata == 1,]), 
             mean(samples_south) * nrow(d[d$strata == 2,])
             )
  
  # true abundance
  truth <- sum(d$eta)
  
  return(as.data.frame(cbind(ice_value, est, truth)))
}

result <- abundance(d, ice_value)
result


# SORRY, I REMOVED THIS PLOT AND HAVEN'T REINTEGRATED IT YET
# plot <- ggplot(sim_dat, aes(X, Y)) +
#   geom_raster(aes(fill = eta)) +
#   geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
#   #geom_hline(aes(yintercept = 50)) +
#   scale_fill_viridis_c() +
#   scale_size_area() +
#   coord_cartesian(expand = FALSE)