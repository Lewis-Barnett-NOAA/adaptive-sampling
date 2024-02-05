set.seed(123)

#Sampling units
n <- 1000

#data frame
predictor_dat <- data.frame(
  expand.grid(X = 1:(n/10), Y = 1:(n/10)),
  depth = rep(rev(seq(1, n/10)), each = n/10)
  #depth = rep(c(1:(n/20), rev(1:(n/20))), n/10) * rev(seq(1, n*10))/max(seq(1, n*10))
)

predictor_dat

strata <- ifelse(predictor_dat$Y >= 50, 1, 2)

ggplot(predictor_dat, aes(X, Y)) +
  geom_tile(aes(fill = depth)) +
  scale_color_gradient2()

mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 300)



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
ordersimsx





cold_pool_value <- sample(ordersimsx$march_sea_ice,1)
cold_pool_value

# Define the ranges for high, mid, and low cold pool values
high_range <- c(0.8, 1)
mid_range <- c(0.4, 0.8)
low_range <- c(0, 0.4)

# Function to determine the operating model based on cold pool value
get_operating_model <- function(cold_pool_value) {
  if (cold_pool_value >= high_range[1] && cold_pool_value <= high_range[2]) {
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
  } else if (cold_pool_value >= mid_range[1] && cold_pool_value <= mid_range[2]) {
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
  } else if (cold_pool_value >= low_range[1] && cold_pool_value <= low_range[2]) {
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

get_operating_model(sample(ordersimsx$march_sea_ice,1))


get_fit <- function(cold_pool_value) {
  # Call get_operating_model function to get the simulation data and plot
  result <- get_operating_model(cold_pool_value)
  

  sim_dat_obs <- result[[3]]
  
  mesh <- make_mesh(sim_dat_obs, xy_cols = c("X", "Y"), cutoff = 0.05)
  
  fit <- sdmTMB(
    observed ~ 1,
    data = sim_dat_obs,
    mesh = mesh,
    family = poisson()
  )

  
  return(list(plot(mesh), result[[1]], sanity(fit)))
}

get_fit(sample(ordersimsx$march_sea_ice,1))
