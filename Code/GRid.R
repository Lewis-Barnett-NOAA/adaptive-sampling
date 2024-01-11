
install.packages("INLA", repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), 
                 dependencies = TRUE)
inla.upgrade()
install.packages("mgcv", dependencies = TRUE)
install.packages("sdmTMB", dependencies = TRUE)
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE)
install.packages("gratia")

## LOAD IN

library(sdmTMB)
library(ggplot2)
library(INLA)
library(mgcv)
library(sdmTMB)
library(dplyr)
library(gratia)
fit <- sdmTMB(
  density ~ s(depth),
  data = pcod_2011, 
  mesh = pcod_mesh_2011,
  family = tweedie(link = "log")
)
fit

## DATA

dat <- readRDS("Data/wcvi-dogfish.rds")
head(dat)

ggplot(dat, aes(longitude, latitude, size = density)) + geom_point()

## UTM LINES

add_utm_columns(dat, ll_names = c("longitude", "latitude"))  

dat <- add_utm_columns(dat, utm_crs = 3156, ll_names = c("longitude", "latitude"))

ggplot(dat, aes(X, Y, size = density)) + geom_point(shape = 21) + coord_fixed()

ggplot(dat, aes(X, Y, size = density, colour = log(density + 1))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~year) + coord_fixed()

dat$log_depth <- log(dat$depth)
dat$year_factor <- as.factor(dat$year)


##MAKING MESH

mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 10)
plot(mesh)

mesh$mesh$n # number of vertices or knots

#ggplot() + inlabru::gg(mesh$mesh) + coord_fixed() +
  #geom_point(aes(X, Y), data = dat, alpha = 0.2, size = 0.5)

#TRouble with that one

#FITTING SPACIAL MODEL
fit_spatial <- sdmTMB(
  density ~ 1, # intercept only
  data = dat,  
  family = tweedie(link = "log"),
  mesh = mesh,
  spatial = "on",
  # silent = FALSE
)

sanity(fit_spatial)

fit_spatial <- run_extra_optimization(fit_spatial)
sanity(fit_spatial)

fit_spatial

##spatial + spatiotemporal fields

fit <- sdmTMB(
  density ~ 0 + year_factor + poly(log_depth, 2),
  data = dat,
  family = tweedie(link = "log"),
  mesh = mesh,
  spatial = "on",
  spatiotemporal = "iid", #< new
  time = "year", #< new
  # silent = FALSE
)

fit
sanity(fit)

# All to fit a spatial and saptiotemporal model to depth

visreg::visreg(fit, xvar = "log_depth")

visreg::visreg(fit, xvar = "log_depth", scale = "response")

g <- ggeffects::ggeffect(fit, "log_depth [3.5:6.7 by=0.05]")
plot(g)

##PREDICTION

wcvi_grid <- readRDS("Data/wcvi-grid.rds")
head(wcvi_grid)

wcvi_grid$log_depth <- log(wcvi_grid$depth)

grid <- purrr::map_dfr(unique(dat$year), ~ tibble(wcvi_grid, year = .x))
grid$year_factor <- as.factor(grid$year)
head(grid)

p0 <- predict(fit)
p0

p <- predict(fit, newdata = grid)
p

# Depth and year effect contribution:
# (Everything not a random field)
ggplot(p, aes(X, Y, fill = exp(est_non_rf))) +
  facet_wrap(~year) +
  geom_raster() +
  coord_fixed()

# Spatial random field:
ggplot(p, aes(X, Y, fill = omega_s)) +
  facet_wrap(~year) +
  geom_raster() +
  scale_fill_gradient2() +
  coord_fixed()

# Spatial-temporal random field:
ggplot(p, aes(X, Y, fill = epsilon_st)) +
  facet_wrap(~year) +
  geom_raster() +
  scale_fill_gradient2() +
  coord_fixed()

# Overall estimate of density in link (log) space:
ggplot(p, aes(X, Y, fill = est)) +
  facet_wrap(~year) +
  geom_raster() +
  coord_fixed()

# Overall estimate of density: (with log-distributed colour)
ggplot(p, aes(X, Y, fill = exp(est))) +
  facet_wrap(~year) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(trans = "log10")

##RESIDUAL CHECKING

dat$resid <- residuals(fit)

ggplot(dat, aes(X, Y, colour = resid)) +
  facet_wrap(~year) +
  geom_point(size = 0.5) +
  coord_fixed() +
  scale_colour_gradient2()

qqnorm(dat$resid)
qqline(dat$resid)

##Fitting an anisotropic model

fit_reml <- update(fit, reml = TRUE)
fit_aniso <- update(fit, reml = TRUE, anisotropy = TRUE)

sanity(fit_reml)
sanity(fit_aniso)

fit_reml
fit_aniso ##Anisotropy vs isotropy, or directional vs equal correlation decay in all directions

plot_anisotropy(fit_aniso)

AIC(fit_reml, fit_aniso)

##INDEX STANDARDIZATION

p <- predict(fit, newdata = grid, return_tmb_object = TRUE)
index <- get_index(p, area = grid$cell_area, bias_correct = FALSE)

ggplot(index, aes(year, est)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab("Year") +
  ylab("Biomass estimate (tonnes)")

index_c <- get_index(p, area = grid$cell_area, bias_correct = TRUE)
index_c$Method <- "Bias correction"

bind_rows(index, index_c) %>%
  ggplot(aes(year, est, fill = Method)) +
  geom_line(aes(colour = Method)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab("Year") +
  ylab("Biomass estimate (tonnes)")


##FIT GAM

fit_gam <- gam(
  formula = density ~ s(depth) + as.factor(year),
  family = tw(link = "log"),
  data = dat
)

summary(fit_gam)

plot(fit_gam, shade = TRUE, residuals = TRUE)

gam.check(fit_gam)

## NEW DATA

pred_gam <- predict(fit_gam, type = "response", newdata = grid)
pred_gam_df <- cbind(grid, pred_gam)

ggplot(pred_gam_df, aes(X, Y, fill = pred_gam)) + geom_raster() +
  scale_fill_viridis_c() + facet_wrap(~year) + coord_fixed() +
  labs(fill = "Log Biomass density\n(kg/km^2)")

## BIOMASS INDEX FROM GAM

sims <- gratia::fitted_samples(fit_gam, n=10, data=grid, 
                               scale="response", seed=9)
sims
sims$year <- grid$year[sims$row]
sims$biomass <- sims$fitted * 4 # expand from density to biomass, given area

level <- 0.95 # specify probability for confidence interval

# Get sum of simulated biomass (density*area) across grid cells, with CI
lwr_fn <- function(x) {as.numeric(quantile(x, probs = (1 - level) / 2))}
upr_fn <- function(x) {as.numeric(quantile(x, probs = 1 - (1 - level) / 2))}

sims_sum <-  sims %>% 
  group_by(year,draw) %>% 
  summarise_at("biomass", list(biomass = sum)) %>%
  group_by(year) %>%
  summarise_at("biomass", list(est = median, # could use mean
                               lwr = lwr_fn,
                               upr = upr_fn))
sims_sum


ggplot(sims_sum, aes(year, est)) + geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Biomass estimate (kg)')


## Fit a GAM analogous to a spatial model in sdmTMB

fit_gam_s <- gam(
  formula = density ~ s(depth) + as.factor(year) + 
    s(X,Y), #<<
  family = tw(link = "log"),
  data = dat
)

plot(fit_gam_s)

plot(fit_gam_s, scheme = 1)

plot(fit_gam_s, scheme = 2)

##Fit a GAM analogous to a spatiotemporal model in sdmTMB

fit_gam_st <- gam(
  formula = density ~ s(depth) + as.factor(year) +
    s(X,Y, by = year), #<<
  family = tw(link = "log"),
  data = dat)

plot(fit_gam_st)

this life