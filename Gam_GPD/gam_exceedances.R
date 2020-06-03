library(fields)
library(INLA)
library(sp)
library(stringr)

DATA="~/research/data/ECAD/ECA_blend_rr/"
setwd("~/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/")
OUT="~/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/data/"
load(paste0(OUT, "german-precipitation.RData"))

load("/Users/worklinda/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/sites_DE_sel_elev.Rdata")
load("/Users/worklinda/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Gam_Gamma/quantile.mat.Rdata")
load("/Users/worklinda/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Gam_Gamma/fitted_quant_75.Rdata")


### Threshold exceedances
exc.gamma     <- c(data_sel-fitted.quant.mat) #stacked site by site, i.e., exc of site 1, exc of site 2, etc...
dat.exc.gamma <- data.frame("exc"=exc.gamma,
                            "lat"=rep(sites_DE_sel$lat, each=nrow(data_sel)),
                            "lon"=rep(sites_DE_sel$lon, each=nrow(data_sel)),
                            "elev"=rep(sites_DE_sel$elevation, each=nrow(data_sel)))

############ subsample from exceedances 
set.seed(22)
subsample_tail.exc <- sample(which(dat.exc.gamma$exc>0), 
                             round(length(which(dat.exc.gamma$exc>0))/10))

dat.exc.gamma.subsample.exc <- dat.exc.gamma[subsample_tail.exc,]

library(evgam)
dat.exc.gamma_gpd <- dat.exc.gamma.subsample.exc
dat.exc.gamma_gpd$exc[dat.exc.gamma_gpd$exc <= 0] <- NA

fmla_gpd                   <- list(exc ~ s(lat, lon)+ s(elev), ~ s(lat, lon)+ s(elev))
m_gpd_lat_lon_elev_sub_exc <- evgam(fmla_gpd, dat.exc.gamma_gpd, family = "gpd")

save(m_gpd_lat_lon_elev_sub_exc, file="data/m_gpd_lat_lon_elev_sub_exc.Rdata")

#Estimated scale and shape parameters for the 4527 locations in the dataset
dat2pred          <- data.frame("lat"=sites_DE_sel$lat, "lon"=sites_DE_sel$lon, "elev"=sites_DE_sel$elevation)
gpd_lat_elev_pred <- predict(m_gpd_lat_lon_elev_sub_exc, dat2pred, type = "response", se.fit = TRUE)
gpd.fitted        <- gpd_lat_elev_pred$fitted
gpd.fitted.se     <- gpd_lat_elev_pred$se.fit
gpd_estimates     <- list("fitted_gpd"=gpd.fitted, "sd_gpd"=gpd.fitted.se)

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "scale"=gpd.fitted[,1],
                 "scale.se"=gpd.fitted.se[,1],
                 "shape"=gpd.fitted[,2],
                 "shape.se"=gpd.fitted.se[,2])

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
png("fitted_GPD_spatial.png", width = 15, height = 12,units = 'in', res = 200)

grid.arrange(gg.scale,gg.scale.se,gg.shape,gg.shape.se, nrow=2,ncol=2)

dev.off()


