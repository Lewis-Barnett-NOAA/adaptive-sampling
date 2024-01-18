library(sdmTMB)
library(ggplot2)
library(INLA)
library(mgcv)
library(sdmTMB)
library(dplyr)
library(gratia)


cold_pool <- readRDS("Data/ice_coldpool_lm.RDS")
cold_pool


#Make Data grid
predictor_dat <- expand.grid(
  X = seq(0, 1, length.out = 100),
  Y = seq(0, 1, length.out = 100)
)

plot(predictor_dat)
# Make mesh
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.5)
plot(mesh)
# Define strata variable based on the Y coordinate (1: North, 2: South)
strata <- ifelse(predictor_dat$Y >= 0.5, 1, 2)

# Operating model
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = poisson(link = "log"),
  range = 0.3,
  sigma_O = 0.4,
  seed = 1,
  B = 1  # B0 = intercept  # Add strata to the operating model
)



# Allocate samples based on strata
samples_north <- sim_dat[strata == 1, ][sample(seq_len(sum(strata == 1)), 20), ]


samples_south <- sim_dat[strata == 2, ][sample(seq_len(sum(strata == 2)), 100), ]


# Combine samples
sim_dat_obs <- rbind(samples_north,samples_south)

# Explore the results

ggplot(sim_dat, aes(X, Y)) +
  geom_raster(aes(fill = exp(eta))) + # mean without observation error
  geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)#+
  #geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") #+
 # geom_text(aes(y = 0.5, x = 0.5, label = "North/South Boundary"), color = "red", vjust = -1)

#Make mesh for sim_dat_obs 
mesh2 <- make_mesh(sim_dat_obs, xy_cols = c("X", "Y"), cutoff = 0.05)
plot(mesh2)
#Thats sick asf

#fit model to sim data samoles
fit <- sdmTMB(
  observed ~ 1,
  data = sim_dat_obs,
  mesh = mesh2,
  family = poisson()
)

fit
sanity(fit)

#check for residuals
sim_dat_obs$resid <- residuals(fit)

ggplot(sim_dat_obs, aes(X, Y, colour = resid)) +
  geom_point(size = 0.5) +
  coord_fixed() +
  scale_colour_gradient2()

qqnorm(sim_dat_obs$resid)
qqline(sim_dat_obs$resid)


#Index
new_data <- sim_dat_obs
p <- predict(fit, newdata = new_data, return_tmb_object = TRUE)
p


#INDex for one year because data has no temporal component
index <- get_index( p, bias_correct = FALSE)


