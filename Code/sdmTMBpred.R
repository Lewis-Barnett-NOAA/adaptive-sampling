library(sdmTMB)
library(ggplot2)
library(INLA)
library(mgcv)
library(sdmTMB)
library(dplyr)
library(gratia)

##DATA
dat <- readRDS("Data/wcvi-dogfish.rds")
head(dat)

ggplot(dat, aes(longitude, latitude, size = density)) + geom_point()

##UTM Lines
add_utm_columns(dat, ll_names = c("longitude", "latitude"))  

dat <- add_utm_columns(dat, utm_crs = 3156, ll_names = c("longitude", "latitude"))

ggplot(dat, aes(X, Y, size = density)) + geom_point(shape = 21) + coord_fixed()

ggplot(dat, aes(X, Y, size = density, colour = log(density + 1))) +
  geom_point(alpha = 0.3) +
  facet_wrap(~year) + coord_fixed()

dat$log_depth <- log(dat$depth)
dat$year_factor <- as.factor(dat$year)

##MESH
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 10)
plot(mesh)

mesh$mesh$n # number of vertices or knots

##SPATIAL + SPATIOTEMPORAL MODEL WITH TWEEDIE DISTRIBUTION
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

##PREDICT ONTO SURVEY GRID
wcvi_grid <- readRDS("Data/wcvi-grid.rds")
head(wcvi_grid)

wcvi_grid$log_depth <- log(wcvi_grid$depth)
head(wcvi_grid)

grid <- purrr::map_dfr(unique(dat$year), ~ tibble(wcvi_grid, year = .x))
grid$year_factor <- as.factor(grid$year)
head(grid)

p0 <- predict(fit) #predict on og data
p0

fit
p <- predict(fit, newdata = grid) #predict on new data
p

##CHECK FOR RESIDUALS
dat$resid <- residuals(fit)

ggplot(dat, aes(X, Y, colour = resid)) +
  facet_wrap(~year) +
  geom_point(size = 0.5) +
  coord_fixed() +
  scale_colour_gradient2()

qqnorm(dat$resid)
qqline(dat$resid)

##CALCULATE INDEX
p <- predict(fit, newdata = grid, return_tmb_object = TRUE)
index <- get_index(p, area = grid$cell_area, bias_correct = FALSE)
index

index_columns <- index[, c('year', 'est')]
index_columns

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

############################
#Sim data from scratch

predictor_dat <- grid
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 10)
sim_dat <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  mesh = mesh,
  family = poisson(link = "log"),
  range = 0.3,
  sigma_O = 0.4,
  seed = 1,
  B = 1 # B0 = intercept
)
head(sim_dat)

#Sample 200 points for fitting
set.seed(1)
sim_dat_obs <- sim_dat[sample(seq_len(nrow(sim_dat)), 200), ]

ggplot(sim_dat, aes(X, Y)) +
  geom_raster(aes(fill = exp(eta))) + # mean without observation error
  geom_point(aes(size = observed), data = sim_dat_obs, pch = 21) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_cartesian(expand = FALSE)


##SAMPLING FUNCTIONS
dat1 <- readRDS("Data/sim_dat_8s.rds")

filledContourPlot.fn <- function(dat, valueName="eta.scaled",...) {
  X <- as.numeric(names(table(dat$X)))
  Y <- as.numeric(names(table(dat$X)))
  
  val <- matrix(NA,nrow=length(Y),ncol=length(X))
  for(xx in 1:nrow(val)) {
    ind <- dat$X==xx
    val[xx,dat$Y[ind]] <- dat[ind,valueName]
  }
  # val <- matrix(dat.yr1$eta,100,100) #alternative way to do it, but make sure it maps correctly to X and Y
  
  filled.contour(x=X,y=Y,val,...)#nlevels = 30,color.palette=colorRampPalette(c("darkseagreen1", "darkgreen")))
  
  dots <- list(...)
  if(!("xlab" %in% names(dots))) mtext("X",side=1,outer=T,line=-2)
  if(!("ylab" %in% names(dots))) mtext("Y",side=2,outer=T,line=-2)
  
  invisible(list(x=X,y=Y,z=val))
}

filledContourPlot.fn(dat1)


##CATCHABILITY FROM GAMMA DISTRIBUTION
catchability.gamma <- function(nn,mm,vv){
  shape=mm^2/vv
  scale=vv/mm
  catchability=rgamma(nn,shape=shape, scale=scale)
  return(catchability)
}

catchability.gamma(10,.8,.01)


##SAMPLE SIM DATA
observationsLN.fn <- function(n, median, cv, catchabilityParams=list(mean=1,var=0.1)) {
  logSigma <- sqrt(log(cv^2+1))
  samp <- rlnorm(n, log(median), logSigma)
  if(is.null(catchabilityParams)) {
    qq <- 1
  } else {
    if(catchabilityParams$var<=0) {
      qq <- catchabilityParams$mean
    } else {
      qq <- catchability.gamma(n, catchabilityParams$mean, catchabilityParams$var)			
    }
  }
  out <- list(observed=qq*samp, catchability=qq, samp=samp)
  return(out)
}

observationsLN.fn(3, 21, 0.2, list(mean=1.0,var=0.0000001))
observationsLN.fn(3, 21, 0.2, list(mean=1.0,var=0.01))
observationsLN.fn(3, 21, 0.2, list(mean=1.0,var=0.0))
observationsLN.fn(3, 21, 0.2, list(mean=0.8,var=-1))
observationsLN.fn(3, 21, 0.2, NULL)


##SAMPLE FROM X AND Y GRID

sampleGrid.fn <- function(xy, dat, obsCV, catchabilityPars,varName="eta.scaled") {
  gridLocs <- paste(dat$X,dat$Y)
  sampLocs <- paste(xy$X,xy$Y)
  out <- dat[gridLocs%in%sampLocs,]
  
  out$observed <- observationsLN.fn(n=length(out$eta), median=out[,varName], obsCV, catchabilityPars)$observed
  return(out)
}

#dat <- data.frame(X=c(1,1,1,2,2,2,3,3,3),Y=c(1,2,3,1,2,3,1,2,3),eta=c(1:9))
 xy <- data.frame(X=c(1,1,2,3,3,3),Y=c(2,3,2,1,2,3))
 obsCV <- 0.1
 catchabilityPars <- list(mean=1,var=0.01)
 sampleGrid.fn(xy, dat, obsCV, catchabilityPars)
sampleGrid.fn(xy,dat1,obsCV,catchabilityPars)


##USING FUNCTIONS
dat <- readRDS("Data/sim_dat_8s.rds")
dat.yr1 <- dat[dat$year==1,]

#plot locations (a grid of 100 X 100)
plot(dat.yr1$X,dat.yr1$Y,pch=20,cex=0.5,xlab="X",ylab="Y",las=1)

#plot depth contours
filledContourPlot.fn(dat.yr1, valueName="depth", nlevels = 25,color.palette=colorRampPalette(c("skyblue", "blue4")))

#plot eta.scaled
filledContourPlot.fn(dat.yr1, valueName="eta.scaled", nlevels = 30,color.palette=colorRampPalette(c("darkseagreen1", "darkgreen")))


#example using sampleGrid

#all grid locations
xyAll <- paste(dat.yr1$X,dat.yr1$Y)
xyAll

#sample from grid locations
xySample <- sample(xyAll, 10, replace=F)
tmp <-  t(as.data.frame(strsplit(xySample,"\\s+")))
rownames(tmp) <- NULL
xy <- data.frame(X=as.numeric(tmp[,1]),Y=as.numeric(tmp[,2]))

obsCV <- 0.1
catchabilityPars <- list(mean=1,var=0.0)
sampleGrid.fn(xy, dat.yr1, obsCV, catchabilityPars)

#One data set with all 6 years
file=paste0("Data/","sim_dat_8s.RDS")
dat <- readRDS(file)
for(yr in unique(dat$year)) {
  dat.yr<- dat[dat$year==yr,]
  mean.eta=mean(dat.yr$eta.scaled)
  filledContourPlot.fn(dat.yr, nlevels = 30,valueName="eta.scaled",zlim=c(0,max(dat$eta.scaled)),main=paste("Year", yr, file, "Mean.Eta", round(mean.eta,4)))
  #aa=locator(n=1)
  #if the margins are not right, close all graphics windows and try again
}

#All data sets 
for(ii in 1:20) {
  file=paste0("Data/sim_dat_",ii,"s.RDS")
  #file=paste0("./Simulated Data-20221223T205153Z-001/Simulated Data/","sim_dat_",ii,"s.RDS")
  dat <- readRDS(file)
  dat.yr1 <- dat[dat$year==2,]
  mean.eta=mean(dat.yr1$eta.scaled)
  filledContourPlot.fn(dat.yr1, nlevels = 30,valueName="eta.scaled",zlim=c(0,120),
                       main=ii)
  aa=locator(n=1)
  #if the margins are not right, close all graphics windows and try again
}


