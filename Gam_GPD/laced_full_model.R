library(fields)
library(INLA)
library(sp)
library(stringr)

load("data/german-precipitation.RData")
load("data/sites_DE_sel_elev.Rdata")
load("data/fitted_quant_75.Rdata")


### Threshold exceedances
exc.gamma     <- c(data_sel-fitted.quant.mat) #stacked site by site, i.e., exc of site 1, exc of site 2, etc...
dat.exc.gamma <- data.frame("exc"=exc.gamma,
                            "lat"=rep(sites_DE_sel$lat, each=nrow(data_sel)),
                            "lon"=rep(sites_DE_sel$lon, each=nrow(data_sel)),
                            "elev"=rep(sites_DE_sel$elevation, each=nrow(data_sel)),
                            "year"=rep(dates_sel$year, ncol(data_sel)),
                            "month"=rep(dates_sel$month, ncol(data_sel)))

############ subsample from exceedances 
set.seed(22)
subsample_tail.exc <- sample(which(dat.exc.gamma$exc>0), 
                             round(length(which(dat.exc.gamma$exc>0))/10))

dat.exc.gamma.subsample.exc <- dat.exc.gamma[subsample_tail.exc,]

library(evgam)
dat.exc.gamma_gpd <- dat.exc.gamma.subsample.exc
dat.exc.gamma_gpd$exc[dat.exc.gamma_gpd$exc <= 0] <- NA

fmla_gpd           <- list(exc ~ s(year)+s(month, bs="cc", k=12)+s(elev)+s(lat,lon), ~ s(lat, lon)+ s(elev))
m_gpd_full_sub_exc <- evgam(fmla_gpd, dat.exc.gamma_gpd, family = "gpd")

save(m_gpd_full_sub_exc, file="data/m_gpd_full_sub_exc.Rdata")

scale.full.mat <- shape.full.mat <- scale.sd.full.mat <- shape.sd.full.mat <- NULL

for(yy in 1:length(unique(dates_sel$year))){
  print(yy)
  dat2pred          <- data.frame("lat"=rep(sites_DE_sel$lat,12),
                                  "lon"=rep(sites_DE_sel$lon,12),
                                  "elev"=rep(sites_DE_sel$elevation,12),
                                  "month"=rep(1:12, each=nrow(sites_DE_sel)),
                                  "year"=rep(unique(dates_sel$year)[yy], 12*nrow(sites_DE_sel)))
  gpd_lat_elev_pred    <- predict(m_gpd_full_sub_exc, dat2pred, type = "response", se.fit = TRUE)
  gpd.fitted           <- gpd_lat_elev_pred$fitted
  gpd.fitted.se        <- gpd_lat_elev_pred$se.fit
  scale.full.mat       <- rbind(scale.full.mat, matrix(gpd.fitted[,1], byrow=TRUE, ncol=length(sites_DE_sel$lon)))
  shape.full.mat       <- rbind(shape.full.mat, matrix(gpd.fitted[,2], byrow=TRUE, ncol=length(sites_DE_sel$lon)))
  scale.sd.full.mat    <- rbind(scale.sd.full.mat, matrix(gpd.fitted.se[,1], byrow=TRUE, ncol=length(sites_DE_sel$lon)))
  shape.sd.full.mat    <- rbind(shape.sd.full.mat, matrix(gpd.fitted.se[,2], byrow=TRUE, ncol=length(sites_DE_sel$lon)))
}

dates_mat <- cbind(rep(unique(dates_sel$year), each=12), rep(1:12, length(unique(dates_sel$year)))) 

fct.merge.scale <- function(ss){
  dates_fitted_scale <- data.frame("year"=dates_mat[,1],
                                   "month"=dates_mat[,2],
                                   "scale"=scale.full.mat[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_scale, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_scale)))
  return(merge.try$scale)
}

fct.merge.scale.sd <- function(ss){
  dates_fitted_scale <- data.frame("year"=dates_mat[,1],
                                   "month"=dates_mat[,2],
                                   "scale.sd"=scale.sd.full.mat[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_scale, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_scale)))
  return(merge.try$scale.sd)
}

fct.merge.shape <- function(ss){
  dates_fitted_shape <- data.frame("year"=dates_mat[,1],
                                   "month"=dates_mat[,2],
                                   "shape"=shape.full.mat[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_shape, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_shape)))
  return(merge.try$shape)
}

fct.merge.shape.sd <- function(ss){
  dates_fitted_shape <- data.frame("year"=dates_mat[,1],
                                   "month"=dates_mat[,2],
                                   "shape.sd"=shape.sd.full.mat[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_shape, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_shape)))
  return(merge.try$shape.sd)
}

#####
fitted.scale.mat.list <- lapply(1:ncol(data_sel), fct.merge.scale) 
fitted.scale.mat      <- do.call(cbind,fitted.scale.mat.list)

fitted.scale.sd.mat.list <- lapply(1:ncol(data_sel), fct.merge.scale.sd) 
fitted.scale.sd.mat      <- do.call(cbind,fitted.scale.sd.mat.list)

fitted.shape.mat.list <- lapply(1:ncol(data_sel), fct.merge.shape) 
fitted.shape.mat      <- do.call(cbind,fitted.shape.mat.list)

fitted.shape.sd.mat.list <- lapply(1:ncol(data_sel), fct.merge.shape.sd) 
fitted.shape.sd.mat      <- do.call(cbind,fitted.shape.sd.mat.list)

fitted.gpd.full <- list("scale"=fitted.scale.mat,
                        "shape"=fitted.shape.mat,
                        "scale.sd"=fitted.scale.sd.mat,
                        "shape.sd"=fitted.shape.sd.mat)
save(fitted.gpd.full, file="data/fitted.gpd.full.Rdata")

##### Plots of the GPD parameters on January 2000
which.idx <- which((dates_sel$year==2000)&(dates_sel$month==1))[1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "scale"=fitted.gpd.full$scale[which.idx,],
                 "scale.se"=fitted.gpd.full$scale.sd[which.idx,],
                 "shape"=fitted.gpd.full$shape[which.idx,],
                 "shape.se"=fitted.gpd.full$shape.sd[which.idx,])

library(ggplot2)                 
gg.scale <- ggplot() +
  geom_point(aes(lon,lat,colour=scale), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="fitted scale")+
  theme_bw()

gg.scale.se <- ggplot() +
  geom_point(aes(lon,lat,colour=scale.se), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="se scale")+
  theme_bw()

gg.shape <- ggplot() +
  geom_point(aes(lon,lat,colour=shape), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="fitted shape")+
  theme_bw()

gg.shape.se <- ggplot() +
  geom_point(aes(lon,lat,colour=shape.se), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="se shape")+
  theme_bw()

library(gridExtra)
png("fitted_GPD_full_jan2000.png", width = 15, height = 12,units = 'in', res = 200)

grid.arrange(gg.scale,gg.scale.se,gg.shape,gg.shape.se, nrow=2,ncol=2)

dev.off()

##### Plots of the GPD parameters on August 2000
which.idx <- which((dates_sel$year==2000)&(dates_sel$month==8))[1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "scale"=fitted.gpd.full$scale[which.idx,],
                 "scale.se"=fitted.gpd.full$scale.sd[which.idx,],
                 "shape"=fitted.gpd.full$shape[which.idx,],
                 "shape.se"=fitted.gpd.full$shape.sd[which.idx,])

library(ggplot2)                 
gg.scale <- ggplot() +
  geom_point(aes(lon,lat,colour=scale), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="fitted scale")+
  theme_bw()

gg.scale.se <- ggplot() +
  geom_point(aes(lon,lat,colour=scale.se), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="se scale")+
  theme_bw()

gg.shape <- ggplot() +
  geom_point(aes(lon,lat,colour=shape), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="fitted shape")+
  theme_bw()

gg.shape.se <- ggplot() +
  geom_point(aes(lon,lat,colour=shape.se), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="se shape")+
  theme_bw()

library(gridExtra)
png("fitted_GPD_full_aug2000.png", width = 15, height = 12,units = 'in', res = 200)

grid.arrange(gg.scale,gg.scale.se,gg.shape,gg.shape.se, nrow=2,ncol=2)

dev.off()

