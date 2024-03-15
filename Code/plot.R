# Plots for cold pool sea ice simulations

library(dplyr)
library(reshape2)

results <- readRDS("results.RDS") #load results

# plot bias and RRMSE of abundance estimates, grouped by scenario ----
results_bias <- mutate(results, bias = est-truth)
range(results_bias$bias)

p_range_bias <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  ggplot(aes(as.factor(range), bias)) + 
  geom_boxplot() +
  labs(x = "Spatial Range", y = "Bias") +
  theme_bw()
p_range_bias
ggsave("Figures/SimFigs/range_bias.pdf")

p_range_rrmse <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  group_by(range) %>%
  summarise(rrmse = (sqrt(mean(bias ^ 2)) / mean(est)) * 100) %>%
  ungroup() %>%
  ggplot(aes(range, rrmse)) + 
  geom_point() +
  geom_line() +
  labs(x = "Spatial Range", y = "RRMSE %") +
  theme_bw()
p_range_rrmse
ggsave("Figures/SimFigs/range_rrmse.pdf")

p_obserr_bias <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  ggplot(aes(as.factor(phi), bias)) + 
  geom_boxplot() +
  labs(x = "Observation Error SD", y = "Bias") +
  theme_bw()
p_obserr_bias
ggsave("Figures/SimFigs/obserr_bias.pdf")

p_obserr_rrmse <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  group_by(phi) %>%
  summarise(rrmse = (sqrt(mean(bias ^ 2)) / mean(est)) * 100) %>%
  ungroup() %>%
  ggplot(aes(phi, rrmse)) + 
  geom_point() +
  geom_line() +
  labs(x = "Observation Error SD", y = "RRMSE %") +
  theme_bw()
p_obserr_rrmse
ggsave("Figures/SimFigs/obserr_rrmse.pdf")           

p_gradient_bias <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  ggplot(aes(as.factor(B1_low), bias)) + 
  geom_boxplot() +
  labs(x = "True Population Density Gradient", y = "Bias") +
  theme_bw()
p_gradient_bias
ggsave("Figures/SimFigs/gradient_bias.pdf")

p_gradient_rrmse <- drop_na(results) %>% 
  mutate(bias = est - truth) %>%
  group_by(B1_low) %>%
  summarise(rrmse = (sqrt(mean(bias ^ 2)) / mean(est)) * 100) %>%
  ungroup() %>%
  ggplot(aes(B1_low, rrmse)) + 
  geom_point() +
  geom_line() +
  labs(x = "True Population Density Gradient", y = "RRMSE %") +
  theme_bw()
p_gradient_rrmse
ggsave("Figures/SimFigs/gradient_rrmse.pdf")


# plot operating model map with observations ----

# exploratory
#range and estimate
ggplot(data = results, aes(x = range, y = est, group = range)) +
  geom_boxplot()
#range and truth
ggplot(data = results, aes(x = range, y = truth, group = range)) +
  geom_boxplot()

#range
results_range <- select(results, range, est, truth)
long_results_r <- melt(results_range, id.vars = "range")
long_results

ggplot(long_results_r,aes(range,value,fill=variable))+
  geom_bar(stat="identity",position="dodge")
ggsave("Figures/SimFigs/rangecomp.pdf")

#phi
results_phi <- select(results, phi, est, truth)
long_results_p <- melt(results_phi, id.vars = "phi")
long_results

ggplot(long_results_p,aes(phi,value,fill=variable))+
  geom_bar(stat="identity",position="dodge")

#b1_low
results_bl <- select(results, B1_low, est, truth)
long_results_bl <- melt(results_bl, id.vars = "B1_low")
long_results

ggplot(long_results_bl,aes(B1_low,value,fill=variable))+
  geom_bar(stat="identity",position="dodge")

#b1_mid
results_bm <- select(results, B1_mid, est, truth)
long_results_bm <- melt(results_bm, id.vars = "B1_mid")
long_results

ggplot(long_results_bm,aes(B1_mid,value,fill=variable))+
  geom_bar(stat="identity",position="dodge")

#b1_high
results_bh <- select(results, B1_high, est, truth)
long_results_bh <- melt(results_bh, id.vars = "B1_high")
long_results

ggplot(long_results_bh,aes(B1_high,value,fill=variable))+
  geom_bar(stat="identity",position="dodge")

#ice_value #this is not an explanatory variable in this example I dont think
results_ice <- select(results, ice_value, est, truth)
long_results_ice <- melt(results_ice, id.vars = "ice_value")
long_results

ggplot(long_results_ice,aes(ice_value,value,fill=variable))+
  geom_bar(stat="identity",position="dodge")

####
sampling <- function(d, ice_value, n) {
  # Define strata based on Y values
  d$strata <- ifelse(d$Y >= 50, 1, 2)
  
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

#representative range values
params[16,]
op_low <- get_operating_model(100000, params[16,])
sim_dat <- sampling(op_low,0.5,100)

range20 <- ggplot(op_low, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21) +
  geom_hline(aes(yintercept = 50)) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)
ggsave("Figures/SimFigs/range20_plot.png", range20)

params[20,]
op_low <- get_operating_model(100000, params[20,])
sim_dat <- sampling(op_low,0.5,100)

range100 <- ggplot(op_low, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21) +
  geom_hline(aes(yintercept = 50)) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)
ggsave("Figures/SimFigs/range100_plot.png", range100)

#representative density values
params[13,]
op_low <- get_operating_model(100000, params[13,])
sim_dat <- sampling(op_low,0.5,100)

wk_grad <- ggplot(op_low, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21) +
  geom_hline(aes(yintercept = 50)) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)
ggsave("Figures/SimFigs/wk_grad_plot.png", wk_grad)

params[88,]
op_low <- get_operating_model(100000, params[88,])
sim_dat <- sampling(op_low,0.5,100)

strng_grad <- ggplot(op_low, aes(X, Y)) +
  geom_raster(aes(fill = eta)) +
  geom_point(aes(size = observed), data = sim_dat, pch = 21) +
  geom_hline(aes(yintercept = 50)) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)
ggsave("Figures/SimFigs/strng_grad_plot.png", strng_grad)