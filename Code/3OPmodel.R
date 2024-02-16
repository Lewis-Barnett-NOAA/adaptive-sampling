library(sdmTMB)
library(ggplot2)
library(dplyr)
library(tidyr)

# COLD POOL SIM FUNCTIONS ----
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
abundance <- function(d, ice_value, n) {
  
  d$strata <- ifelse(d$Y >= 50, 1, 2) # TODO: redefine strata (not evenly split somehow?)
  
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
  
  # true abundance
  truth <- sum(d$eta)
  
  return(as.data.frame(cbind(n, ice_value, est, truth)))
}


# Sea Ice - Cold pool simulations ----
set.seed(123)

# Sample size
n <- 100

#Sampling units in domain
N <- 1000

#make data frame including fake "depth" covariate gradient
predictor_dat <- data.frame(
  expand.grid(X = 1:(N/10), Y = 1:(N/10)),
  depth = rep(rev(seq(1, N/10)), each = N/10)
)

ggplot(predictor_dat, aes(X, Y)) +
  geom_tile(aes(fill = depth)) +
  scale_color_gradient2()

#get triangulated mesh to simulate from
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 200)
plot(mesh)

#define parameters to loop over in operating models
ranges <- seq(20, 100, 20) # spatial range (higher = smoother, lower = patchier)
phis <- seq(0.1, 0.5, 0.1) # observation error SD
B1_lows <- seq(0, 0.3, 0.1) # slope of depth-density relationship in low cold pool scenario
params <- as.data.frame(expand.grid(range=ranges, phi=phis, B1_low=B1_lows))
params$B1_mid <- params$B1_low + 0.3
params$B1_high <- params$B1_low + 0.6

# replicate parameter df once per simulation replicate
n_rep <- 100
params <- replicate_df(params, time_name = "sim_id", time_values = 1:n_rep)

# define empty object to house results dataframe
results <- data.frame(matrix(NA, nrow(params), ncol(params) + 5))
colnames(results) <- c(colnames(params), "coldpool", "n", "ice_value", "est", "truth")

# simulate cold pool extent from random march sea ice values
m <- readRDS("Data/ice_coldpool_lm.RDS") # load linear fit cold pool:sea ice
x <- data.frame(march_sea_ice = seq(0, max(m$model$march_sea_ice), by = 0.05))
order_sims <- simulateX(m,  nsim = n_rep, X = x)
ordersimsx <- cbind(order_sims, x)
ordersimsx[ordersimsx < 0] <- 0 # set negative cold pool extents to 0
ordersimsx

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

# Visualize the simulated data relative to the thresholds for these categories
vars <- ordersimsx %>% 
          pivot_longer(
            cols = starts_with("sim_"),
            names_to = "sim_id",
            names_prefix = "sim_",
            values_to = "coldpool"
          )
        
p_vars <- ggplot(vars, aes(as.factor(march_sea_ice), coldpool)) + 
  geom_boxplot() +
  labs(x = "March Sea Ice Proportion", y = "Cold Pool Area") +
  theme_bw()
p_vars
ggsave("Figures/sea_ice_coldpool_sims.pdf")

# Loop over parameter scenarios for the operating model ----
for(i in 1:nrow(params)){
  # Draw a random march sea ice proportion
  ice_value <- sample(ordersimsx$march_sea_ice,1)
  # Draw a cold pool extent simulated from a random march sea ice proportion
  cold_pool_value <- as.numeric(ordersimsx[ordersimsx$march_sea_ice == ice_value,][sample(1:(ncol(ordersimsx)-1), 1)])
  
  d <- get_operating_model(cold_pool_value, params[i, ])
  
  # Append results to params
  results[i, ] <- cbind(params[i, ], 
                        cold_pool_value, 
                        abundance(d, ice_value, n)
                        )
}

results
#saveRDS(results, "results.RDS")

# plot bias and RRMSE of abundance estimates, grouped by scenario ----
#results <- readRDS("results.RDS") #load results

p_range_bias <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  ggplot(aes(as.factor(range), bias)) + 
  geom_boxplot() +
  labs(x = "Spatial Range", y = "Bias") +
  theme_bw()
p_range_bias
ggsave("Figures/range_bias.pdf")
  
p_range_rrmse <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  group_by(range) %>%
  summarise(rrmse = sqrt(mean((truth - est) ^ 2) / sum(est ^ 2))) %>%
  ungroup() %>%
  ggplot(aes(range, rrmse)) + 
  geom_point() +
  geom_line() +
  labs(x = "Spatial Range", y = "RRMSE") +
  theme_bw()
p_range_rrmse
ggsave("Figures/range_rrmse.pdf")

p_obserr_bias <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  ggplot(aes(as.factor(phi), bias)) + 
  geom_boxplot() +
  labs(x = "Observation Error SD", y = "Bias") +
  theme_bw()
p_obserr_bias
ggsave("Figures/obserr_bias.pdf")

p_obserr_rrmse <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  group_by(phi) %>%
  summarise(rrmse = sqrt(mean((truth - est) ^ 2) / sum(est ^ 2))) %>%
  ungroup() %>%
  ggplot(aes(phi, rrmse)) + 
  geom_point() +
  geom_line() +
  labs(x = "Observation Error SD", y = "RRMSE") +
  theme_bw()
p_range_rrmse
ggsave("Figures/obserr_rrmse.pdf")           

p_gradient_bias <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  ggplot(aes(as.factor(B1_low), bias)) + 
  geom_boxplot() +
  labs(x = "True Population Density Gradient", y = "Bias") +
  theme_bw()
p_gradient_bias
ggsave("Figures/gradient_bias.pdf")

p_gradient_rrmse <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  group_by(B1_low) %>%
  summarise(rrmse = sqrt(mean((truth - est) ^ 2) / sum(est ^ 2))) %>%
  ungroup() %>%
  ggplot(aes(B1_low, rrmse)) + 
  geom_point() +
  geom_line() +
  labs(x = "True Population Density Gradient", y = "RRMSE") +
  theme_bw()
p_gradient_rrmse
ggsave("Figures/gradient_rrmse.pdf")

# SORRY, I REMOVED THIS PLOT AND HAVEN'T REINTEGRATED IT YET ----
# plot <- ggplot(sim_dat, aes(X, Y)) +
#   geom_raster(aes(fill = eta)) +
#   geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
#   #geom_hline(aes(yintercept = 50)) +
#   scale_fill_viridis_c() +
#   scale_size_area() +
#   coord_cartesian(expand = FALSE)


# OLD EXPERIMENTAL CODE
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