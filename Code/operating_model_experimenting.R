set.seed(123)

N = 1000

predictor_dat <- data.frame(
  expand.grid(X = 1:(N/10), Y = 1:(N/10)),
  temperature = rep(rev(seq(1, N/10)), each = N/10)
)
predictor_dat$temperature_scaled <- scale(predictor_dat$temperature)

#get triangulated mesh to simulate from
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 200)
#plot(mesh)

# #define parameters to loop over in operating models
# ranges <- c(10, 50, 100) # spatial range (higher = smoother, lower = patchier)
# phis <- c(0.01, 0.05, 0.1) # observation error or dispersion depending on distribution
# B1_lows <- c(0, 0.2, 0.4) # slope of temperature-density relationship in low cold pool scenario
# params <- as.data.frame(expand.grid(range=ranges, phi=phis, B1_low=B1_lows))
# params$B1_mid <- params$B1_low + 0.3
# params$B1_high <- params$B1_low + 0.6

sim_dat <- sdmTMB_simulate(
  formula = ~ 1 + temperature_scaled,
  data = predictor_dat,
  mesh = mesh,
  family = tweedie(),
  range = 200,
  phi = 0.1, # dispersion
  sigma_O = 0.2,
  tweedie_p = 1.9,
  B = c(0.2, 0.2) # B0 = intercept, B1 = a1 slope
)

hist(exp(sim_dat$eta))
hist(sim_dat$observed)

ggplot(sim_dat, aes(X, Y)) +
 geom_raster(aes(fill = exp(eta))) +
                        geom_hline(aes(yintercept = 50)) +
                        scale_fill_viridis_c() +
                        scale_size_area() +
                        coord_cartesian(expand = FALSE)