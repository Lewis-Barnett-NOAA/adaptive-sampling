library(sdmTMB)
library(ggplot2)

set.seed(123)

#Sampling units
n <- 1000

#make data frame including fake "depth" gradient
predictor_dat <- data.frame(
  expand.grid(X = 1:(n/10), Y = 1:(n/10)),
  depth = rep(rev(seq(1, n/10)), each = n/10)
  #depth = rep(c(1:(n/20), rev(1:(n/20))), n/10) * rev(seq(1, n*10))/max(seq(1, n*10))
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

## Depth coefficient slope = 0.9 (strong gradient)
sim_dat <- sdmTMB_simulate(
  formula = ~ 1 + depth,
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(link = "identity"),
  B = c(0.1, 0.9), # B0 = intercept, B1 = depth coefficient slope
  range = ranges[i],
  sigma_O = 0.2,
  phi = 5, # SD of observation error in this (Gaussian) case
  seed = 42
)

# Shift values so that minimum abundance is 0 (or fit with a different family)
sim_dat$eta <- (sim_dat$eta - min(sim_dat$eta))
sim_dat$observed <- (sim_dat$observed - min(sim_dat$observed))

# Plot simulated "true" data
ggplot(sim_dat, aes(X, Y)) +
  geom_tile(aes(fill = eta)) +
  scale_color_gradient2()


# Plot simulated "observed" data
ggplot(sim_dat, aes(X, Y)) +
  geom_tile(aes(fill = observed)) +
  scale_color_gradient2()



#Depth coefficient slope = 0.6 (mid gradient)
sim_dat2 <- sdmTMB_simulate(
  formula = ~ 1 + depth,
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(link = "identity"),
  B = c(0.1, 0.6), # B0 = intercept, B1 = depth coefficient slope
  range = ranges[i],
  sigma_O = 0.2,
  phi = 5, # SD of observation error in this (Gaussian) case
  seed = 42
)

sim_dat2$eta <- (sim_dat2$eta - min(sim_dat2$eta))
sim_dat2$observed <- (sim_dat2$observed - min(sim_dat2$observed))

# Plot simulated "observed" data
ggplot(sim_dat2, aes(X, Y)) +
  geom_tile(aes(fill = observed)) +
  scale_color_gradient2()




##depth coefficient slope = 0.2 (weak gradient)
sim_dat3 <- sdmTMB_simulate(
  formula = ~ 1 + depth,
  data = predictor_dat,
  mesh = mesh,
  family = gaussian(link = "identity"),
  B = c(0.1, 0.2), # B0 = intercept, B1 = depth coefficient slope
  range = ranges[i],
  sigma_O = 0.2,
  phi = 5, # SD of observation error in this (Gaussian) case
  seed = 42
)

sim_dat3$eta <- (sim_dat3$eta - min(sim_dat3$eta))
sim_dat3$observed <- (sim_dat3$observed - min(sim_dat3$observed))


# Plot simulated "observed" data
ggplot(sim_dat3, aes(X, Y)) +
  geom_tile(aes(fill = observed)) +
  scale_color_gradient2()



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

# Draw a cold pool extent simulated from a random march sea ice proportion
cold_pool_value <- as.numeric(ordersimsx[ordersimsx$march_sea_ice == sample(ordersimsx$march_sea_ice,1),][sample(1:(ncol(ordersimsx)-1), 1, replace = TRUE)])
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

strata <- ifelse(predictor_dat$Y >= 50, 1, 2) # Define strata

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
   
    samples_north <- sim_dat[strata == 1, ][sample(seq_len(sum(strata == 1)), 20), ]
    samples_south <- sim_dat[strata == 2, ][sample(seq_len(sum(strata == 2)), 100), ]
    
    
    # Combine samples
    sim_dat_obs <- rbind(samples_north,samples_south)
   
     plot1 <- ggplot(sim_dat, aes(X, Y)) +
      geom_raster(aes(fill = eta)) +
      geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
      scale_fill_viridis_c() +
      scale_size_area() +
      coord_cartesian(expand = FALSE)
    return(list("Strong Gradient",plot1, sim_dat_obs))
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
    
    samples_north <- sim_dat[strata == 1, ][sample(seq_len(sum(strata == 1)), 30), ]
    samples_south <- sim_dat[strata == 2, ][sample(seq_len(sum(strata == 2)), 90), ]
    
    
    # Combine samples
    sim_dat_obs <- rbind(samples_north,samples_south)
    
    plot <- ggplot(sim_dat, aes(X, Y)) +
      geom_raster(aes(fill = eta)) +
      geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
      scale_fill_viridis_c() +
      scale_size_area() +
      coord_cartesian(expand = FALSE)
    
    return(list("Mid Gradient",plot, sim_dat_obs))
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
    
    samples_north <- sim_dat[strata == 1, ][sample(seq_len(sum(strata == 1)), 60), ]
    samples_south <- sim_dat[strata == 2, ][sample(seq_len(sum(strata == 2)), 60), ]
    
    
    # Combine samples
    sim_dat_obs <- rbind(samples_north,samples_south)
    
    plot <- ggplot(sim_dat, aes(X, Y)) +
      geom_raster(aes(fill = eta)) +
      geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
      scale_fill_viridis_c() +
      scale_size_area() +
      coord_cartesian(expand = FALSE)
    return(list("Low Operating Model",plot, sim_dat_obs))
  } else {
    return("No Matching Operating Model")
  }
}

obs <- get_operating_model(sample(ordersimsx$march_sea_ice,1))





get_estimate <- function(sea_ice_prop, obs) {
 
  return(list(sea_ice_prop, cold_pool_value, est))
}

get_estimate(sample(ordersimsx$march_sea_ice,1))
