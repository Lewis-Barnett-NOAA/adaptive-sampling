library(sdmTMB)
library(ggplot2)

set.seed(123)

# Simulating from scratch ----

# scale number of primary sampling units in sampling universe
n <- 1000 # produces 10,000 per year
# set number of years
nyears <- 6
  
# make fake predictor (depth) and sampling locations:
predictor_dat <- data.frame(
  expand.grid(X = 1:(n/10), Y = 1:(n/10)),
  depth = rep(c(1:(n/20), rev(1:(n/20))), n/10) * rep(c(1,1.1,0.9,1.2,0.95,0.8), each = n*10), # * rlnorm(n*10, sdlog = 0.1) 
  year = rep(1:nyears, each = n*10)            
)
saveRDS(predictor_dat, "grid_depth.RDS")

# plot domain for single year
ggplot(dplyr::filter(predictor_dat, year == 1), aes(X, Y)) +
  geom_tile(aes(fill = depth)) +
  scale_color_gradient2()
#ggsave("depth.png")

# ggplot(filter(predictor_dat, year == 1), aes(X, Y, colour = depth)) +
#   geom_point() +
#   scale_color_gradient(low = "lightblue", high = "darkblue")

#---- fit model
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 300)

ranges <- seq(5, 100, 5) # specify values of range to simulate from

for(i in 1:length(ranges)){
  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + depth,
    data = predictor_dat,
    mesh = mesh,
    family = gaussian(link = "identity"),
    time = "year",
    B = c(0.1, 0.7), # B0 = intercept, B1 = depth coefficient slope
    range = ranges[i],
    rho = NULL, # or value between -1 and 1 to add spatiotemporal correlation
    sigma_O = 0.2,
    sigma_E = 10,
    phi = 5, # SD of observation error in this (Gaussian) case
    seed = 42
  )

  sim_dat$eta <- (sim_dat$eta - min(sim_dat$eta))
  sim_dat$observed <- (sim_dat$observed - min(sim_dat$observed))
  
  # Plot simulated "true" data
  ggplot(sim_dat, aes(X, Y)) +
    geom_tile(aes(fill = eta)) +
    facet_wrap(~year) +
    scale_color_gradient2()
  ggsave(paste0("true_density_range_",i,".png"))
  
  # Plot simulated "observed" data
  ggplot(sim_dat, aes(X, Y)) +
    geom_tile(aes(fill = observed)) +
    facet_wrap(~year) +
    scale_color_gradient2()
  ggsave(paste0("observed_density_range_",i,".png"))
  
  saveRDS(sim_dat, paste0("sim_dat_",i,".RDS"))
}