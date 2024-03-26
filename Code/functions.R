# simulate observations of sea ice and cold pool
simulateX <- function(object, nsim = 1, seed = NULL, X, ...) {
  object$fitted.values <- predict(object, X)
  simulate(object = object, nsim = nsim, seed = seed, ...)
}

# Function to determine the operating model based on sea ice value
get_operating_model <- function(cold_pool_value, params) {
  if (cold_pool_value >= high_cp[1] && cold_pool_value <= high_cp[2]) {
    sim_dat <- sdmTMB_simulate(
      formula = ~ 1 + depth,
      data = predictor_dat,
      mesh = mesh,
      family = gaussian(link = "identity"),
      B = c(0.1, params$B1_high), # B0 = intercept, B1 = depth coefficient slope
      range = params$range,
      sigma_O = 0.2,
      phi = params$phi # SD of observation error in this (Gaussian) case
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
      B = c(0.1, params$B1_mid), # B0 = intercept, B1 = depth coefficient slope
      range = params$range,
      sigma_O = 0.2,
      phi = params$phi # SD of observation error in this (Gaussian) case
    )
    
    sim_dat$eta <- (sim_dat$eta - min(sim_dat$eta))
    sim_dat$observed <- (sim_dat$observed - min(sim_dat$observed))
    
    return(sim_dat)
    
  } else if (cold_pool_value >= 0 && cold_pool_value <= low_cp[2]) {
    sim_dat <- sdmTMB_simulate(
      formula = ~ 1 + depth,
      data = predictor_dat,
      mesh = mesh,
      family = gaussian(link = "identity"),
      B = c(0.1, params$B1_low), # B0 = intercept, B1 = depth coefficient slope
      range = params$range,
      sigma_O = 0.2,
      phi = params$phi # SD of observation error in this (Gaussian) case
    )
    
    sim_dat$eta <- (sim_dat$eta - min(sim_dat$eta))
    sim_dat$observed <- (sim_dat$observed - min(sim_dat$observed))
    
    return(sim_dat)
    
  } else {
    return("No Matching Operating Model")
  }
}

# Simulate sampling
abundance <- function(d, ice_value, n, design = c(
  "adaptive stratified","proportional stratified",
  "simple random", "south stratum only"
)) {
  
  d$strata <- ifelse(d$Y > 50, 1, 2)
  
  if(design == "adaptive stratified") {
       if (ice_value >= high_ice[1] && ice_value <= high_ice[2]) {
      samples_north <- sample(d[d$strata == 1, "observed"], n*0.1)
      samples_south <- sample(d[d$strata == 2, "observed"], n*0.9)
      
    } else if (ice_value >= mid_ice[1] && ice_value <= mid_ice[2]) {
      samples_north <- sample(d[d$strata == 1, "observed"], n*0.25)
      samples_south <- sample(d[d$strata == 2, "observed"], n*0.75)
      
    } else if (ice_value >= low_ice[1] && ice_value <= low_ice[2]) {
      samples_north <- sample(d[d$strata == 1, "observed"], n*0.5)
      samples_south <- sample(d[d$strata == 2, "observed"], n*0.5)
    }
    # estimate abundance
    est <- sum(mean(samples_north) * nrow(d[d$strata == 1,]), 
               mean(samples_south) * nrow(d[d$strata == 2,])
    )
  }
 
  if(design == "proportional stratified") {
    samples_north <- sample(d[d$strata == 1, "observed"], n*0.5)
    samples_south <- sample(d[d$strata == 2, "observed"], n*0.5)
    est <- sum(mean(samples_north) * nrow(d[d$strata == 1,]), 
               mean(samples_south) * nrow(d[d$strata == 2,])
    )
  }
  
  if(design == "simple random") {
    samples <- sample(d[, "observed"], n)
    est <- mean(samples) * nrow(d)
  }
  
  if(design == "south stratum only") {
    samples_north <- sample(d[d$strata == 1, "observed"], 0)
    samples_south <- sample(d[d$strata == 2, "observed"], n)
    est <- mean(samples_south) * nrow(d)
  }
  
  # true abundance
  truth <- sum(d$eta)
  
  return(as.data.frame(cbind(n, ice_value, est, truth, design)))
}

# simulate observations for plotting
sampling <- function(d, ice_value, n) {
  # Define strata based on Y values
  d$strata <- ifelse(d$Y > 50, 1, 2)
  
  # Define sampling proportions based on ice_value
  if (ice_value >= high_ice[1] && ice_value <= high_ice[2]) {
    prop_north <- 0.1
    prop_south <- 0.9
  } else if (ice_value >= mid_ice[1] && ice_value <= mid_ice[2]) {
    prop_north <- 0.25
    prop_south <- 0.75
  } else if (ice_value >= low_ice[1] && ice_value <= low_ice[2]) {
    prop_north <- 0.5
    prop_south <- 0.5
  }
  
  # Sample row indices for north and south strata
  indices_north <- sample(which(d$strata == 1), floor(n * prop_north))
  indices_south <- sample(which(d$strata == 2), floor(n * prop_south))
  
  # Extract sampled rows from the data frame
  samples_north <- d[indices_north, ]
  samples_south <- d[indices_south, ]
  
  # Add a region column to each data frame
  samples_north$region <- "north"
  samples_south$region <- "south"
  
  # Combine the data frames into a single data frame
  combined_samples <- dplyr::bind_rows(samples_north, samples_south)
  
  return(combined_samples)
}