#-----  fit index model to the simulated data
sim_dat <- readRDS("sim_dat.RDS")

# sample 100 locations using simple random sampling
set.seed(123)
sim_dat_obs <- sim_dat[sample(seq_len(nrow(sim_dat)), 100), ]

# plot samples over surface of mean without observation error (taken as "true")
ggplot(sim_dat, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
  facet_wrap(~year) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)

# fit new model
mesh <- make_mesh(sim_dat_obs, xy_cols = c("X", "Y"), cutoff = 0.1)

fit <- sdmTMB(
  formula = observed ~ 0 + as.factor(year), # or year + depth if including covariate
  data = sim_dat_obs,
  mesh = mesh,
  time = "year",
  family = gaussian(link = "identity"), # default, could be left out of class example to simplify or leave ambiguity
  spatial = "on", # c("on", "off")
  spatiotemporal = "iid", # c("iid", "ar1", "rw", "off")
)

fit

#---- inspect model and fit, get abundance index (following vignette...)
#https://pbs-assess.github.io/sdmTMB/articles/index-standardization.html