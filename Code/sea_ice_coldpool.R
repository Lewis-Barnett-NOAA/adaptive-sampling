library("sf")
library("dplyr")
library("akgfmaps")
library("ggplot2")

# API ----
#ice_poly <- st_read("https://nsidc.org/api/mapservices/NSIDC/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=NSIDC:g02135_polyline_n") # &filter=time@filter_from=1982-03
# filter to march of every year from ice_ext_2$timestamp

# # Data downloaded by FTP 2023-04-25: ftp://sidads.colorado.edu/pub/DATASETS/NOAA/G02135/
# # More information: https://nsidc.org/data/g02135/versions/3
# ice_ras <- read.csv("Data/Sea_ice_data/NSIDC.csv")
# #TODO: replace with same product for full range of latitudes and redo (these data are cut off halfway through EBS)
# 
# # filter data by region, time, sea ice cover (only available from 1988 onward for March in this version and need lower latitudes)
# ice_df <- ice_ras %>% filter(Month == 3, 
#                              Year > 1981) %>%
#                       na.omit()
# 
# # Minimum latitude by year
# # TODO: try Seaice >= 0.1
# ice_ext <- ice_df %>% filter(Seaice >= 0.15) %>% group_by(Year) %>% summarise(Extent = min(Latitude))
# colnames(ice_ext) <- tolower(colnames(ice_ext))

# compare to coldpool extent ----
coldpool <- coldpool:::cold_pool_index
colnames(coldpool) <- tolower(colnames(coldpool))

# # combine and fit regressions ----
# df <- left_join(ice_ext, coldpool) %>% na.omit()
# 
# m_lte2 <- lm(log(area_lte2_km2) ~ extent, df)
# m_lte1 <- lm(area_lte1_km2 ~ log(extent), df)
# m_lte0 <- lm(area_lte0_km2 ~ log(extent), df)
# m_lteminus1 <- lm(area_lteminus1_km2 ~ log(extent), df)
# 
# print(summary(m_lte2))
# print(summary(m_lte1))
# print(summary(m_lte0))
# print(summary(m_lteminus1))
# # r-square degrades with decreasing temperature threshold
# 
# #best relationship plotted
# plot(df$extent, log(df$area_lte2_km2))
# abline(m_lte2)
# 
# # TODO: separate by outer/inner/middle domain boxes from stratum shapefile and calculate extent specific to each? ----
# ice <- ice_df %>% filter(Seaice >= 0.15)
# ice_sf <- st_as_sf(ice, coords = c("Longitude", "Latitude"))
# st_crs(ice_sf) <- 4326
# 
# ebs <- akgfmaps::get_base_layers("ebs", "EPSG:4326")

# TODO: use bathymetry to define inner/outer/middle and then do spatial join between areas and ice_sf
# Note that there may be complications as bathymetry includes geometries for islands and no inner boundary for land
# ggplot() +
#   geom_sf(data = ebs$akland) +
#   geom_sf(data = ebs$bathymetry) +
#   geom_sf(data = ice_sf, color = "red") +
#   geom_sf(data = ebs$graticule, color = "grey70", alpha = 0.5) +
#   coord_sf(xlim = ebs$plot.boundary$x,
#            ylim = ebs$plot.boundary$y) +
#   scale_x_continuous(name = "Longitude",
#                      breaks = ebs$lon.breaks) +
#   scale_y_continuous(name = "Latitude",
#                      breaks = ebs$lat.breaks) +
#   theme_bw()

#st_filter(ice_sf, ebs$bathymetry) #could use this if make bathymetry into polygon with 0 depth/land

# TODO: Could also calculate probability of temp being below 2, 1, 0 etc given ice extent for each grid cell of current survey or finer spatial scale?


# Same regression for whole area, but using the proportion of EBS with sea ice, ----
# using data from ERDAPP extracted by Matt Callahan on 2023-4-27 
ice_prop <- read.csv("Data/Sea_ice_data/ebs_march_ice.csv")
colnames(ice_prop) <- tolower(colnames(ice_prop))

df2 <- left_join(ice_prop, coldpool) %>% na.omit()

m_lte2_prop <- lm(area_lte2_km2 ~ march_sea_ice, df2)
m_lte2_prop_log <- lm(log(area_lte2_km2) ~ march_sea_ice, df2)
m_lte1_prop <- lm(area_lte1_km2 ~ march_sea_ice, df2)
m_lte0_prop <- lm(area_lte0_km2 ~ march_sea_ice, df2)
m_lteminus1_prop <- lm(area_lteminus1_km2 ~ march_sea_ice, df2)

print(summary(m_lte2_prop))
print(summary(m_lte2_prop_log))
print(summary(m_lte1_prop))
print(summary(m_lte0_prop))
print(summary(m_lteminus1_prop))

#best relationship plotted
m1 <- summary(m_lte2_prop)
m1_log <- summary(m_lte2_prop_log)
saveRDS(m_lte2_prop, "Data/Sea_ice_data/ice_coldpool_lm.RDS")

pdf("Data/Sea_ice_data/ice_coldpool_regression.pdf")
plot(df2$march_sea_ice, df2$area_lte2_km2, xlim=c(0,0.7), 
     ylim=c(0,max(df2$area_lte2_km2)*1.05), xaxs="i", yaxs="i",
     xlab="Proportion of EBS with sea ice in prior March",
     ylab="Cold pool extent index (sq-km)")
abline(m_lte2_prop, col="blue")
#text(df2$march_sea_ice, df2$area_lte2_km2,labels=df2$year) #maybe a little autocorrelation?
march_sea_ice <- seq(0,0.7, by = 0.05)
ci <- predict(m_lte2_prop, newdata=data.frame(march_sea_ice), interval="confidence",
                         level = 0.95)
matlines(march_sea_ice, ci[,2:3], col = "blue", lty=2)
mtext(paste0("R-sq = ",round(m1$r.squared,2)),adj = 0)
mtext(paste0("sigma = ",round(m1$sigma,0)), adj = 0.3)
mtext(paste0("y ~ ",round(m1$coefficients[1],2)," + ", 
             round(m1$coefficients[2],2),"x"), adj = 1)
dev.off()

# log y provides higher R-square, but relationship really isn't linear,
# thus would be better fit with a concave/saturating
pdf("Data/Sea_ice_data/ice_coldpool_regression_log.pdf")
plot(df2$march_sea_ice, log(df2$area_lte2_km2), xlim=c(0,0.7), xaxs="i",
     xlab="Proportion of EBS with sea ice in prior March",
     ylab="Cold pool extent index (log sq-km)")
abline(m_lte2_prop_log, col = "blue")
ci_log <- predict(m_lte2_prop_log, newdata=data.frame(march_sea_ice), interval="confidence",
              level = 0.95)
matlines(march_sea_ice, ci_log[,2:3], col = "blue", lty=2)
mtext(paste0("R-sq = ",round(m1_log$r.squared,2)),adj = 0)
mtext(paste0("sigma = ",round(m1_log$sigma,2)), adj = 0.5)
mtext(paste0("y ~ ",round(m1_log$coefficients[1],2)," + ", 
             round(m1_log$coefficients[2],2),"x"), adj = 1)
dev.off()

# example simulation ----
# note that you'll want to set predictions < 0 to 0
sims <- simulate(m_lte2_prop, nsim = 10)
plot(m_lte2_prop$model$march_sea_ice, sims$sim_1)
# compile multiple simulations into a dataframe
cbind(m_lte2_prop$model$march_sea_ice, sims)
