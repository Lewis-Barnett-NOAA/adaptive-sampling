library(sdmTMB)
library(ggplot2)

set.seed(123)

N = 1000

predictor_dat <- data.frame(
  expand.grid(X = 1:(N/10), Y = 1:(N/10)),
  temperature = rep(rev(seq(1, N/10)), each = N/10)
)
predictor_dat$temperature_scaled <- scale(predictor_dat$temperature)

#get triangulated mesh to simulate from
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 350)
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
  range = 10,
  phi = 0.5, # dispersion
  sigma_O = 0.2,
  tweedie_p = 1.9,
  B = c(0.2, 0.2) # B0 = intercept, B1 = a1 slope
)

ggplot(sim_dat, aes(X, Y)) +
 geom_raster(aes(fill = exp(eta))) +
                        geom_hline(aes(yintercept = 50)) +
                        scale_fill_viridis_c() +
                        scale_size_area() +
                        coord_cartesian(expand = FALSE)

hist(exp(sim_dat$eta))
hist(sim_dat$observed)

# plot operating model map with observations ----


#strongest gradient with medium range and medium observation error
params <- 
ice_value <- sample(ordersimsx$march_sea_ice,1)
cold_pool_value <- as.numeric(ordersimsx[ordersimsx$march_sea_ice == ice_value,][sample(1:(ncol(ordersimsx)-1), 1)])

d <- get_operating_model_tw(cold_pool_value, params[88, ])

sim_dat <- sampling(d,ice_value,35)
#op <- get_operating_model_tw(150000, params[88,])

strng_grad <- ggplot(d, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21, color = "white") +
  geom_hline(aes(yintercept = 50), color = "red") +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_fixed(ratio = 1,expand = FALSE)+
  labs(title = "Strong Gradient") +  # Add title for Plot A
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(hjust = 0.5)
  )

#intermediate
op <- get_operating_model_tw(cold_pool_value, params[53,])
sim_dat <- sampling(op,ice_value,35)

mid_grad <- ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21, color = "white") +
  geom_hline(aes(yintercept = 50), color = "red") +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_fixed(ratio = 1, expand = FALSE)+
  labs(title = "Mid Gradient")+
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(hjust = 0.5)
  )

#weakest gradient with medium range and medium observation error
op <- get_operating_model_tw(cold_pool_value, params[8,])
sim_dat <- sampling(op,ice_value,35)

wk_grad <- ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21, color = "white") +
  geom_hline(aes(yintercept = 50), color = "red") +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_fixed(ratio = 1, expand = FALSE)+
  labs(title = "Weak Gradient")+
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(hjust = 0.5)
  )

gradient <- plot_grid(strng_grad, mid_grad, wk_grad, labels = "AUTO", ncol = 3)
ggsave("Figures/SimFigs/gradient_vis.pdf")

#medium gradient with lowest range and medium observation error
op <- get_operating_model_tw(cold_pool_value, params[46,])
sim_dat <- sampling(op,ice_value,35)

low_range <- ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21, color = "white") +
  geom_hline(aes(yintercept = 50), color = "red") +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_fixed(ratio = 1, expand = FALSE)+
  labs(title = "Low Range") +  # Add title for Plot A
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(hjust = 0.5)
  )



#mid range
op <- get_operating_model_tw(cold_pool_value, params[48,])
sim_dat <- sampling(op,ice_value,35)

mid_range <- ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21, color = "white") +
  geom_hline(aes(yintercept = 50), color = "red") +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_fixed(ratio = 1, expand = FALSE)+
  labs(title = "Mid Range") +  # Add title for Plot A
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(hjust = 0.5)
  )
#high range
op <- get_operating_model_tw(cold_pool_value, params[50,])
sim_dat <- sampling(op,ice_value,35)

high_range <- ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21, color = "white") +
  geom_hline(aes(yintercept = 50), color = "red") +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_fixed(ratio = 1, expand = FALSE)+
  labs(title = "High Range") +  # Add title for Plot A
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(hjust = 0.5)
  )

range <- plot_grid(low_range, mid_range, high_range, labels = "AUTO", ncol = 3)
ggsave("Figures/SimFigs/range_vis.pdf")
#lowest observation error, medium gradient and range

op <- get_operating_model_tw(cold_pool_value, params[43,])
sim_dat <- sampling(op,ice_value,35)

low_obv_err <- ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21, color = "white") +
  geom_hline(aes(yintercept = 50), color = "red") +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_fixed(ratio = 1, expand = FALSE)+
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(hjust = 0.5)
  )

#highest observation error, medium gradient and range
op <- get_operating_model_tw(cold_pool_value, params[158,])
sim_dat <- sampling(op,ice_value,35)

high_obv_err <- ggplot(op, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21, color = "white") +
  geom_hline(aes(yintercept = 50), color = "red") +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_fixed(ratio = 1, expand = FALSE)+
  labs(title = "High Observation Error") +  # Add title for Plot A
  theme(
    legend.position = "none",  # Remove legend
    plot.title = element_text(hjust = 0.5)
  )

obserr <- plot_grid(low_obv_err, high_obv_err, labels = "AUTO", ncol = 1)
ggsave("Figures/SimFigs/obv_error_vis.pdf")

# Extract the legend from one plot
legend <- get_legend(
  ggplot(d, aes(X, Y)) +
    geom_raster(aes(fill = eta)) +
    geom_point(aes(size = observed), data = sim_dat, pch = 21, color = "white") +
    geom_hline(aes(yintercept = 50), color = "red") +
    scale_fill_viridis_c() +
    scale_size_area() +
    coord_fixed(ratio = 1, expand = FALSE) +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.5, 'lines'),  # Shrink legend size
      legend.text = element_text(size = 8)   # Shrink legend text
    )
)


combined_plot <- plot_grid(
  plot_grid(strng_grad, mid_grad, wk_grad, high_range, mid_range, low_range, ncol = 3)
)

combined_plot

combined_plot_leg <- plot_grid(
  combined_plot,
  legend,
  ncol = 2,
  rel_widths = c(3, 0.4)  # Adjust the relative widths to minimize empty space
)

combined_plot_leg

ggsave("Figures/SimFigs/group_plot.pdf")

#######
#observation error histograms

low_phi <- filter(results_example, results_example[,2]==0.01)
low_phi

high_phi <- filter(results_example, results_example[,2]==0.3)
high_phi

pdf("Figures/SimFigs/phi_abund.pdf")
par(mfrow= c(1,2))

hist(low_phi$est)
hist(high_phi$est)

dev.off()

