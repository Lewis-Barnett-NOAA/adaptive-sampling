source("Code/functions.R")
library(sdmTMB)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(cowplot)

# General simulation settings ----
set.seed(77)

# Sample size
n <- 30

#Sampling units in domain
N <- 1000

#make data frame including fake "temperature" covariate gradient
predictor_dat <- data.frame(
  expand.grid(X = 1:(N/10), Y = 1:(N/10)),
  temperature = rep(rev(seq(1, N/10)), each = N/10)
)
predictor_dat$temperature_scaled <- scale(predictor_dat$temperature)

#get triangulated mesh to simulate from
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), type = "cutoff_search", n_knots = 350)
plot(mesh)

#define parameters to loop over in operating models
ranges <- c(5, 50, 100, 200) # spatial range (higher = smoother, lower = patchier)
phis <- c(0.01, 0.1, 0.3, 0.5) # observation error or dispersion depending on distribution
B1_lows <- c(-0.25, -0.1, 0, 0.1, 0.25) # slope of temperature-density relationship in low cold pool scenario
params <- as.data.frame(expand.grid(range=ranges, phi=phis, B1_low=B1_lows))
params$B1_mid <- params$B1_low + 0.3
params$B1_high <- params$B1_low + 0.6

# replicate parameter df once per simulation replicate
n_rep <- 150
params <- replicate_df(params, time_name = "sim_id", time_values = 1:n_rep)
saveRDS(params, "parameters.RDS")

# define empty object to house results dataframe
results_adapt <- data.frame(matrix(NA, nrow(params), ncol(params) + 6))
colnames(results_adapt) <- c(colnames(params), "coldpool", "n", "ice_value", "est", "truth", "design")
results_noadapt <- results_sonly <- results_srs <- results_sextrap <- results_adapt_perf <- results_adapt 


# simulate cold pool extent from random march sea ice values
m <- readRDS("Data/ice_coldpool_lm.RDS") # load linear fit cold pool:sea ice
x <- data.frame(march_sea_ice = seq(0, max(m$model$march_sea_ice), by = 0.05))
order_sims <- simulateX(m,  nsim = n_rep, X = x)
ordersimsx <- cbind(order_sims, x)
ordersimsx[ordersimsx < 0] <- 0 # set negative cold pool extents to 0
ordersimsx

# Define the ranges for high, mid, and low sea ice values to guide sample allocation
ice_cat <- quantile(m$model$march_sea_ice, probs = c(0, .5, 1))
ice_cat
low_ice <- c(0, ice_cat[2])
mid_ice <- c(ice_cat[2], ice_cat[3])
high_ice <- c(ice_cat[3], 1)

# Define the ranges for high, mid, and low cold pool extent to determine the operating model
cp_cat <- quantile(m$model$area_lte2_km2, probs = c(0, .333, .666, 1))
cp_cat
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

# Loop over parameter scenarios for the operating model ----
for(i in 1:nrow(params)){
  # Draw a random march sea ice proportion
  ice_value <- sample(ordersimsx$march_sea_ice,1)
  # Draw a cold pool extent simulated from a random march sea ice proportion
  cold_pool_value <- as.numeric(ordersimsx[ordersimsx$march_sea_ice == ice_value,][sample(1:(ncol(ordersimsx)-1), 1)])
  
  d <- get_operating_model_tw(cold_pool_value, params[i, ])
  
  # Append results to params
  results_adapt[i, ] <- bind_cols(params[i, ], 
                              cold_pool_value, 
                              abundance(d, ice_value, cold_pool_value, n, design = "adaptive stratified")
                              )
  results_adapt_perf[i, ] <- bind_cols(params[i, ], 
                                  cold_pool_value, 
                                  abundance(d, ice_value, cold_pool_value, n, design = "adaptive stratified perfect")
  )
  results_noadapt[i, ] <- bind_cols(params[i, ], 
                              cold_pool_value, 
                              abundance(d, ice_value, cold_pool_value, n, design = "proportional stratified")
                             )
  results_srs[i, ] <- bind_cols(params[i, ], 
                              cold_pool_value, 
                              abundance(d, ice_value, cold_pool_value, n, design = "simple random")
                            )
  results_sonly[i, ] <- bind_cols(params[i, ], 
                             cold_pool_value, 
                             abundance(d, ice_value, cold_pool_value, n, design = "south stratum only")
  )
  results_sextrap[i, ] <- bind_cols(params[i, ], 
                                  cold_pool_value, 
                                  abundance(d, ice_value, cold_pool_value, n, design = "south stratum extrapolated")
  )
}

results <- bind_rows(results_adapt, results_sonly, results_noadapt, results_srs, results_sextrap, results_adapt_perf)
saveRDS(results, "results_tw_p9_omega2_nrep150_n30_newparams_scaled.RDS")


# Plots for cold pool sea ice simulations -----

#results <- readRDS("results_tw_p9_omega2_nrep100_n50_newparams_scaled.RDS") #load results

p_gradient <-
  ggplot(predictor_dat, aes(X, Y)) +
    geom_tile(aes(fill = temperature_scaled)) +
    scale_color_gradient2() +
    theme_bw()
p_gradient
ggsave("Figures/SimFigs/gradient_fake_covar.pdf")

p_ice <- ggplot(vars, aes(as.factor(march_sea_ice), coldpool)) + 
  geom_boxplot() +
  labs(x = "March Sea Ice Proportion", y = "Cold Pool Area") +
  theme_bw()
p_ice
ggsave("Figures/SimFigs/sea_ice_coldpool_sims.pdf")

# plot bias and RRMSE of abundance estimates, grouped by scenario ----
# separate panels for the south only scenarios
res <- drop_na(results) %>% 
  filter(!design %in% c("south stratum only", "south stratum extrapolated", 
                        "adaptive stratified perfect")) %>% 
  mutate(bias = est - truth)
res_s <- drop_na(results) %>% 
  filter(design %in% c("south stratum only", "south stratum extrapolated")) %>% 
  mutate(bias = est - truth)

p_range_bias <- res %>% 
  ggplot(aes(as.factor(range), bias, fill = as.factor(design))) + 
  geom_boxplot() + #outlier.shape=NA
  #coord_cartesian(ylim = quantile(res$bias, c(0.225, 0.799))) +
  labs(x = "Spatial Range", y = "Bias") +
  theme_bw()
p_range_bias_s <- res_s %>% 
  ggplot(aes(as.factor(range), bias, fill = as.factor(design))) + 
  geom_boxplot() + 
  #coord_cartesian(ylim = quantile(res_s$bias, c(0.115, 0.879))) +
  labs(x = "Spatial Range", y = "Bias") +
  theme_bw()
plot_grid(p_range_bias, p_range_bias_s, labels = "AUTO", ncol = 1)
ggsave("Figures/SimFigs/range_bias.pdf")

p_range_rrmse <- res %>% 
  group_by(range, design) %>%
  summarise(rrmse = (sqrt(mean((est - truth) ^ 2)) / mean(est)) * 100) %>%
  ungroup() %>%
  ggplot(aes(range, rrmse, group = design, color = design)) + 
  geom_point() +
  geom_line() +
  labs(x = "Spatial Range", y = "RRMSE %") +
  theme_bw()
p_range_rrmse_s <- res_s %>% 
  group_by(range, design) %>%
  summarise(rrmse = (sqrt(mean((est - truth) ^ 2)) / mean(est)) * 100) %>%
  ungroup() %>%
  ggplot(aes(range, rrmse, group = design, color = design)) + 
  geom_point() +
  geom_line() +
  labs(x = "Spatial Range", y = "RRMSE %") +
  theme_bw()
plot_grid(p_range_rrmse, p_range_rrmse_s, labels = "AUTO", ncol = 1)
ggsave("Figures/SimFigs/range_rrmse.pdf")

p_obserr_bias <- res %>% 
  ggplot(aes(as.factor(phi), bias, fill = as.factor(design))) + 
  geom_boxplot() + #outlier.shape=NA
  #coord_cartesian(ylim = quantile(res$bias, c(0.225, 0.799))) +
  labs(x = "Dispersion", y = "Bias") +
  theme_bw()
p_obserr_bias_s <- res_s %>% 
  ggplot(aes(as.factor(phi), bias, fill = as.factor(design))) + 
  geom_boxplot() + 
  #coord_cartesian(ylim = quantile(res_s$bias, c(0.05, 0.877))) +
  labs(x = "Dispersion", y = "Bias") +
  theme_bw()
plot_grid(p_obserr_bias, p_obserr_bias_s, labels = "AUTO", ncol = 1)
ggsave("Figures/SimFigs/obserr_bias.pdf")

p_obserr_rrmse <- res %>% 
  group_by(phi, design) %>%
  summarise(rrmse = (sqrt(mean((est - truth) ^ 2)) / mean(est)) * 100) %>%
  ungroup() %>%
  ggplot(aes(phi, rrmse, group = design, color = design)) + 
  geom_point() +
  geom_line() +
  labs(x = "Dispersion", y = "RRMSE %") +
  theme_bw()
p_obserr_rrmse_s <- res_s %>% 
  group_by(phi, design) %>%
  summarise(rrmse = (sqrt(mean((est - truth) ^ 2)) / mean(est)) * 100) %>%
  ungroup() %>%
  ggplot(aes(phi, rrmse, group = design, color = design)) + 
  geom_point() +
  geom_line() +
  labs(x = "Dispersion", y = "RRMSE %") +
  theme_bw()
plot_grid(p_obserr_rrmse, p_obserr_rrmse_s, labels = "AUTO", ncol = 1)
ggsave("Figures/SimFigs/obserr_rrmse.pdf")           

p_gradient_bias <- res %>% 
  ggplot(aes(as.factor(B1_low), bias, fill = as.factor(design))) + 
  geom_boxplot() + 
  #coord_cartesian(ylim = quantile(res$bias, c(0.12, 0.88))) +
  labs(x = "Temperature Coefficient at Low Cold Pool Values", y = "Bias") +
  theme_bw()
p_gradient_bias_s <- res_s %>% 
  ggplot(aes(as.factor(B1_low), bias, fill = as.factor(design))) + 
  geom_boxplot() + 
  #coord_cartesian(ylim = quantile(res_s$bias, c(0.1, 0.87))) +
  labs(x = "Temperature Coefficient at Low Cold Pool Values", y = "Bias") +
  theme_bw()
plot_grid(p_gradient_bias, p_gradient_bias_s, labels = "AUTO", ncol = 1)
ggsave("Figures/SimFigs/gradient_bias.pdf")

p_gradient_rrmse <- res %>% 
  group_by(B1_low, design) %>%
  summarise(rrmse = (sqrt(mean((est - truth) ^ 2)) / mean(est)) * 100) %>%
  ungroup() %>%
  ggplot(aes(B1_low, rrmse, group = design, color = design)) + 
  geom_point() +
  geom_line() +
  labs(x = "Temperature Coefficient at Low Cold Pool Values", y = "RRMSE %") +
  theme_bw()
p_gradient_rrmse_s <- res_s %>% 
  group_by(B1_low, design) %>%
  summarise(rrmse = (sqrt(mean((est - truth) ^ 2)) / mean(est)) * 100) %>%
  ungroup() %>%
  ggplot(aes(B1_low, rrmse, group = design, color = design)) + 
  geom_point() +
  geom_line() +
  labs(x = "Temperature Coefficient at Low Cold Pool Values", y = "RRMSE %") +
  theme_bw()
plot_grid(p_gradient_rrmse, p_gradient_rrmse_s, labels = "AUTO", ncol = 1)
ggsave("Figures/SimFigs/gradient_rrmse.pdf")


# plot estimate and truth across parameters ----

# # exploratory
# #range and estimate
# ggplot(data = results, aes(x = range, y = est, group = range)) +
#   geom_boxplot()
# #range and truth
# ggplot(data = results, aes(x = range, y = truth, group = range)) +
#   geom_boxplot()

# #range
# results_range <- select(results, range, est, truth)
# long_results_r <- melt(results_range, id.vars = "range")
# 
# ggplot(long_results_r,aes(range,value,fill=variable))+
#   geom_bar(stat="identity",position="dodge")
# #ggsave("Figures/SimFigs/rangecomp.pdf")
# 
# #phi
# results_phi <- select(results, phi, est, truth)
# long_results_p <- melt(results_phi, id.vars = "phi")
# 
# ggplot(long_results_p,aes(phi,value,fill=variable))+
#   geom_bar(stat="identity",position="dodge")
# 
# #b1_low
# results_b1 <- select(results, B1_low, est, truth)
# long_results_b1 <- melt(results_b1, id.vars = "B1_low")
# 
# ggplot(long_results_b1,aes(B1_low,value,fill=variable))+
#   geom_bar(stat="identity",position="dodge")
# 
# #b1_mid
# results_bm <- select(results, B1_mid, est, truth)
# long_results_bm <- melt(results_bm, id.vars = "B1_mid")
# 
# ggplot(long_results_bm,aes(B1_mid,value,fill=variable))+
#   geom_bar(stat="identity",position="dodge")
# 
# #b1_high
# results_bh <- select(results, B1_high, est, truth)
# long_results_bh <- melt(results_bh, id.vars = "B1_high")
# 
# ggplot(long_results_bh,aes(B1_high,value,fill=variable))+
#   geom_bar(stat="identity",position="dodge")
# 
# #ice_value
# results_ice <- select(results, ice_value, est, truth)
# long_results_ice <- melt(results_ice, id.vars = "ice_value")
# 
# ggplot(long_results_ice,aes(ice_value,value,fill=variable))+
#   geom_bar(stat="identity",position="dodge")
# 
# 
# # plot operating model map with observations ----
# 
# #representative range values
# op_low <- get_operating_model_tw(100000, params[16,])
# sim_dat <- sampling(op_low,0.5,100)
# 
# range20 <- ggplot(op_low, aes(X, Y)) +
#   geom_raster(aes(fill = eta)) +
#   geom_point(aes(size = observed), data = sim_dat, pch = 21) +
#   geom_hline(aes(yintercept = 50)) +
#   scale_fill_viridis_c() +
#   scale_size_area() +
#   coord_cartesian(expand = FALSE)
# ggsave("Figures/SimFigs/range20_plot.png", range20)
# 
# params[20,]
# op_low <- get_operating_model_tw(100000, params[20,])
# sim_dat <- sampling(op_low,0.5,100)
# 
# range100 <- ggplot(op_low, aes(X, Y)) +
#   geom_raster(aes(fill = eta)) +
#   geom_point(aes(size = observed), data = sim_dat, pch = 21) +
#   geom_hline(aes(yintercept = 50)) +
#   scale_fill_viridis_c() +
#   scale_size_area() +
#   coord_cartesian(expand = FALSE)
# ggsave("Figures/SimFigs/range100_plot.png", range100)
# 
# #representative density values
# params[13,]
# op_low <- get_operating_model_tw(100000, params[13,])
# sim_dat <- sampling(op_low,0.5,100)
# 
# wk_grad <- ggplot(op_low, aes(X, Y)) +
#   geom_raster(aes(fill = eta)) +
#   geom_point(aes(size = observed), data = sim_dat, pch = 21) +
#   geom_hline(aes(yintercept = 50)) +
#   scale_fill_viridis_c() +
#   scale_size_area() +
#   coord_cartesian(expand = FALSE)
# ggsave("Figures/SimFigs/wk_grad_plot.png", wk_grad)
# 
# params[88,]
# op_low <- get_operating_model_tw(100000, params[88,])
# sim_dat <- sampling(op_low,0.5,100)
# 
# strng_grad <- ggplot(op_low, aes(X, Y)) +
#   geom_raster(aes(fill = eta)) +
#   geom_point(aes(size = observed), data = sim_dat, pch = 21) +
#   geom_hline(aes(yintercept = 50)) +
#   scale_fill_viridis_c() +
#   scale_size_area() +
#   coord_cartesian(expand = FALSE)
# ggsave("Figures/SimFigs/strng_grad_plot.png", strng_grad)