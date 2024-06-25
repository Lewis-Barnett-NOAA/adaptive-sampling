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

# plot operating model map with observations ----


#strongest gradient with medium range and medium observation error
params
ice_value <- sample(ordersimsx$march_sea_ice,1)
cold_pool_value <- as.numeric(ordersimsx[ordersimsx$march_sea_ice == ice_value,][sample(1:(ncol(ordersimsx)-1), 1)])

d <- get_operating_model_tw(cold_pool_value, params[88, ])

sim_dat <- sampling(d,ice_value,100)
#op <- get_operating_model_tw(150000, params[88,])

ggplot(d, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21) +
  geom_hline(aes(yintercept = 50)) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)
#weakest gradient with medium range and medium observation error
op <- get_operating_model_tw(150000, params[8,])

ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_hline(aes(yintercept = 50)) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)

#medium gradient with lowest range and medium observation error
op <- get_operating_model_tw(150000, params[46,])

ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_hline(aes(yintercept = 50)) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)
#medium gradient with highest range and medium observation error
op <- get_operating_model_tw(150000, params[50,])

ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_hline(aes(yintercept = 50)) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)
#lowest observation error, medium gradient and range
op <- get_operating_model_tw(150000, params[43,])

ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_hline(aes(yintercept = 50)) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)
#highest observation error, medium gradient and range
op <- get_operating_model_tw(150000, params[158,])

ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_hline(aes(yintercept = 50)) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)