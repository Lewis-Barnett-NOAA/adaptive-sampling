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

saveRDS(predictor_dat, "grid_depth.RDS")

ggplot(predictor_dat, aes(X, Y)) +
  geom_tile(aes(fill = depth)) +
  scale_color_gradient2()

#Thats what the hell I'm talking about


mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 300)

plot(mesh)

ranges <- seq(5, 100, 5)

for(i in 1:length(ranges)){
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
  ggsave(paste0("true_density_range_",i,".png"))
  
  # Plot simulated "observed" data
  ggplot(sim_dat, aes(X, Y)) +
    geom_tile(aes(fill = observed)) +
    scale_color_gradient2()
  ggsave(paste0("observed_density_range_",i,".png"))
  
  saveRDS(sim_dat, paste0("sim_dat_",i,".RDS"))
}


sim_dat <- readRDS("sim_dat/sim_dat_1.RDS")

set.seed(123)
sim_dat_obs <- sim_dat[sample(seq_len(nrow(sim_dat)), 100), ]

# plot samples over surface of mean without observation error (taken as "true")
ggplot(sim_dat, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)


mesh <- make_mesh(sim_dat_obs, xy_cols = c("X", "Y"), cutoff = 0.1)



###this doesnt work, need to figure out how to remove time before fitting
fit <- sdmTMB(
  formula = observed ~ 0, # or year + depth if including covariate
  data = sim_dat_obs,
  mesh = mesh,
  family = gaussian(link = "identity"), # default, could be left out of class example to simplify or leave ambiguity
  spatial = "on", # c("on", "off")
)

fit