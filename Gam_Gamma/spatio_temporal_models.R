library(fields)
library(INLA)
library(sp)
library(stringr)

DATA="~/research/data/ECAD/ECA_blend_rr/"
setwd("~/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/")
OUT="~/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/data/"
load(paste0(OUT, "german-precipitation.RData"))

load("/Users/worklinda/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/sites_DE_sel_elev.Rdata")
### 14030*4527 matrix
anom.training.positive <- data_sel

# replace zeros with NAs
anom.training.positive[anom.training.positive == 0] <- NA

#gam with latitude and longitude in both the location and scale
library(mgcv)

dat   <- data.frame("obs"=as.vector(anom.training.positive),
                    "site"=rep(1:ncol(anom.training.positive), each= nrow(anom.training.positive)),
                    "lat"=rep(sites_DE_sel$lat, each= nrow(anom.training.positive)),
                    "lon"=rep(sites_DE_sel$lon, each= nrow(anom.training.positive)),
                    "elev"=rep(sites_DE_sel$elevation, each= nrow(anom.training.positive)),
                    "month"=rep(dates_sel$month, ncol(anom.training.positive)),
                    "year"=rep(dates_sel$year, ncol(anom.training.positive)))

# Keep one site out of 40: selection made randomly
set.seed(22)
site2keep        <- sample(1:ncol(anom.training.positive), round(ncol(anom.training.positive)/40))

#"homogeneous" subsampling from the map
load("/Users/worklinda/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/SpatialSample.Rdata")
site2keep <- spatial.sample$id
##################################################################################################################################
###                                                    fixed shape                                                             ###
##################################################################################################################################
## fully spatial
mod.scale.spatial    <- gam(list(obs~s(lat,lon)+s(elev),
                             ~1), data = dat,
                        subset= which(dat$site %in% site2keep),
                        family=gammals, select=TRUE)
save(mod.scale.spatial, file="mod.scale.spatial.Rdata")
rm(mod.scale.spatial)

## fully temporal
mod.scale.temporal   <- gam(list(obs~s(year)+s(month, bs="cc", k=12),
                                 ~1), data = dat,
                            subset= which(dat$site %in% site2keep),
                            family=gammals, select=TRUE)
save(mod.scale.temporal, file="mod.scale.temporal.Rdata")
rm(mod.scale.temporal)

##spatio-temporal with years
mod.scale.st   <- gam(list(obs~s(year)+s(elev)+s(lat,lon),
                           ~1), data = dat,
                      subset= which(dat$site %in% site2keep),
                      family=gammals, select=TRUE)
save(mod.scale.st, file="mod.scale.st.Rdata")
rm(mod.scale.st)

##spatio-temporal with months
mod.scale.st.month   <- gam(list(obs~s(month, bs="cc", k=12)+s(elev)+s(lat,lon),
                                 ~1), data = dat,
                        subset= which(dat$site %in% site2keep),
                        family=gammals, select=TRUE)
save(mod.scale.st.month, file="mod.scale.st.month.Rdata")
rm(mod.scale.st.month)

## spatio-temporal model with year and months
mod.scale.st <- gam(list(obs~s(year)+s(month, bs="cc", k=12)+s(elev)+s(lat,lon),
                         ~1), data = dat,
                    subset= which(dat$site %in% site2keep),
                    family=gammals, select=TRUE)
save(mod.scale.st, file="mod.scale.st.year.months.Rdata")

###fitted spatio-temporal model with "fixed shape"seasons"
####### create two seasons instead of modelling a cyclic monthly effect

dat   <- data.frame("obs"=as.vector(anom.training.positive),
                    "site"=rep(1:ncol(anom.training.positive), each= nrow(anom.training.positive)),
                    "lat"=rep(sites_DE_sel$lat, each= nrow(anom.training.positive)),
                    "lon"=rep(sites_DE_sel$lon, each= nrow(anom.training.positive)),
                    "elev"=rep(sites_DE_sel$elevation, each= nrow(anom.training.positive)),
                    "month"=rep(dates_sel$month, ncol(anom.training.positive)),
                    "year"=rep(dates_sel$year, ncol(anom.training.positive)))

dat$season <- ifelse(dat$month %in% seq(5,9), "summer", "winter")


mod.scale.st.season   <- gam(list(obs~s(year)+as.factor(season)+s(elev)+s(lat,lon),
                                  ~1), data = dat,
                             subset= which(dat$site %in% site2keep),
                             family=gammals, select=TRUE)
save(mod.scale.st.season, file="mod.scale.st.season.additive.Rdata")
rm(mod.scale.st.season)

###########################
### fitted models
###########################
###fitted spatial model with fixed shape
fitted.Gamma <- predict(mod.scale.spatial, type="response", newdata = data.frame("lat"=spatial.sample$lat,
                                                                                 "lon"=spatial.sample$lon,
                                                                                 "elev"=sites_DE_sel$elevation[spatial.sample$id]))
fitted.shape <- rep(unique(fitted.Gamma[,2]), each=length(site2keep))
fitted.scale <- fitted.Gamma[,1]

# dat.select   <- dat[which(dat$site %in% site2keep),]
# lat2plot     <- unique(dat.select$lat[!is.na(dat.select$obs)]) #113 unique values
# #lon2plot     <- unique(dat.select$lon[!is.na(dat.select$obs)]) #111 unique values (instead of 113... three sites have the exact same lon)
# elev2plot    <- unique(dat.select$elev[!is.na(dat.select$obs)]) #113 unique values
# 
# idx.unique <- NULL
# for(i in 1:nrow(spatial.sample)){
#   idx.unique <- c(idx.unique, which((sites_DE_sel$lat==spatial.sample$lat[i])&(sites_DE_sel$lon==spatial.sample$lon[i])))
# }
# lon2plot <- sites_DE_sel$lon[idx.unique]
# 
# idx.unique.dat <- NULL
# for(i in 1:length(lat2plot)){
#   idx.unique.dat <- c(idx.unique.dat, which((dat.select$lat[!is.na(dat.select$obs)]==lat2plot[i])&
#                                               (dat.select$elev[!is.na(dat.select$obs)]==elev2plot[i]))[1]) #add condition on year and month if spatio-temporal model
# }
# 
# fitted.scale <- fitted.Gamma[idx.unique.dat,1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat)
df2plot <- data.frame("lon2plot"=spatial.sample$lon,
                      "lat2plot"=spatial.sample$lat,
                      "scale"=fitted.scale)

library(ggplot2)                 
ggplot() +
  geom_point(aes(lon,lat),colour="gray72", data = df, cex=.2)+
  geom_point(data = df2plot,aes(lon2plot,lat2plot,colour=scale))+
  scale_colour_gradientn(colours=rainbow(150), name="fitted scale")+
  theme_bw()

###fitted temporal model with fixed shape
library("viridis")

year2pred <- seq(min(dat$year),max(dat$year),length.out = 100)

yeffect <- NULL
for(mm in 1:12){
  yeffect <- cbind(yeffect, predict(mod.scale.temporal, type="response", newdata=data.frame("month"=mm, 
                                                                                            "year"=year2pred))[,1])
}

pdf("year_effect.pdf", width = 7, height = 5)
par(mar=c(3,3.2,1.5,0.5),mgp=c(1.6,0.5,0),font.main=1.3,cex=1.3,cex.main=1)

col2use <- plasma(12)
plot(year2pred, yeffect[,1], type="l", ylab="Fitted scale", xlab="Time", col=col2use[1], ylim=c(range(yeffect)[1],range(yeffect)[2]))
text(x=year2pred[2],y=yeffect[2,1], labels = as.character(1), cex=0.6)
for(mm in 2:12){
  lines(year2pred, yeffect[,mm], col=col2use[mm])
  text(x=year2pred[2],y=yeffect[2,mm], labels = as.character(mm), cex=0.6)
}
# legend("bottomright", legend=month.abb,
#        col=col2use, lty=1, cex=0.6)

dev.off()

###fitted spatio-temporal (with years) model with fixed shape
fitted.Gamma <- predict(mod.scale.st, type="response", newdata = data.frame("lat"=spatial.sample$lat,
                                                                            "lon"=spatial.sample$lon,
                                                                            "elev"=sites_DE_sel$elevation[spatial.sample$id],
                                                                            "year"=rep(2000, nrow(spatial.sample))))
fitted.shape <- rep(unique(fitted.Gamma[,2]), each=length(site2keep))
fitted.scale <- fitted.Gamma[,1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat)
df2plot <- data.frame("lon2plot"=spatial.sample$lon,
                      "lat2plot"=spatial.sample$lat,
                      "scale"=fitted.scale)

library(ggplot2)                 
ggplot() + ggtitle("2000") +
  geom_point(aes(lon,lat),colour="gray72", data = df, cex=.2)+
  geom_point(data = df2plot,aes(lon2plot,lat2plot,colour=scale))+
  scale_colour_gradientn(colours=rainbow(150), name="fitted scale")+
  theme_bw()

###fitted spatio-temporal (with months) model with fixed shape
#######January
fitted.Gamma <- predict(mod.scale.st.month, type="response", newdata = data.frame("lat"=spatial.sample$lat,
                                                                                  "lon"=spatial.sample$lon,
                                                                                  "elev"=sites_DE_sel$elevation[spatial.sample$id],
                                                                                  "month"=rep(1, nrow(spatial.sample))))
fitted.shape <- rep(unique(fitted.Gamma[,2]), each=length(site2keep))
fitted.scale <- fitted.Gamma[,1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat)
df2plot <- data.frame("lon2plot"=spatial.sample$lon,
                      "lat2plot"=spatial.sample$lat,
                      "scale"=fitted.scale)

library(ggplot2)                 
gg_jan <- ggplot() + ggtitle("January") +
  geom_point(aes(lon,lat),colour="gray72", data = df, cex=.2)+
  geom_point(data = df2plot,aes(lon2plot,lat2plot,colour=scale))+
  scale_colour_gradientn(colours=rainbow(150), name="fitted scale")+
  theme_bw()

#######August
fitted.Gamma <- predict(mod.scale.st.month, type="response", newdata = data.frame("lat"=spatial.sample$lat,
                                                                                  "lon"=spatial.sample$lon,
                                                                                  "elev"=sites_DE_sel$elevation[spatial.sample$id],
                                                                                  "month"=rep(8, nrow(spatial.sample))))
fitted.shape <- rep(unique(fitted.Gamma[,2]), each=length(site2keep))
fitted.scale <- fitted.Gamma[,1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat)
df2plot <- data.frame("lon2plot"=spatial.sample$lon,
                      "lat2plot"=spatial.sample$lat,
                      "scale"=fitted.scale)

library(ggplot2)                 
gg_aug <- ggplot() + ggtitle("August") +
  geom_point(aes(lon,lat),colour="gray72", data = df, cex=.2)+
  geom_point(data = df2plot,aes(lon2plot,lat2plot,colour=scale))+
  scale_colour_gradientn(colours=rainbow(150), name="fitted scale")+
  theme_bw()

library(gridExtra)
grid.arrange(gg_jan,gg_aug, nrow=1,ncol=2)

###fitted spatio-temporal (with "seasons") model with fixed shape
####### Summer&Winter
fitted.Gamma <- predict(mod.scale.st.season, type="response", newdata = data.frame("lat"=rep(spatial.sample$lat,2),
                                                                                   "lon"=rep(spatial.sample$lon,2),
                                                                                   "elev"=rep(sites_DE_sel$elevation[spatial.sample$id],2),
                                                                                   "season"=rep(c("summer","winter"), each=nrow(spatial.sample)),
                                                                                   "year"=rep(2000, 2*nrow(spatial.sample))))
fitted.shape        <- rep(unique(fitted.Gamma[,2]), each=length(site2keep))
fitted.scale.summer <- fitted.Gamma[1:nrow(spatial.sample),1]
fitted.scale.winter <- fitted.Gamma[(nrow(spatial.sample)+1):(2*nrow(spatial.sample)),1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat)
####### Summer 2000
df2plot <- data.frame("lon2plot"=spatial.sample$lon,
                      "lat2plot"=spatial.sample$lat,
                      "scale"=fitted.scale.summer)

library(ggplot2)                 
gg_summer <- ggplot() + ggtitle("Summer 2000") +
  geom_point(aes(lon,lat),colour="gray72", data = df, cex=.2)+
  geom_point(data = df2plot,aes(lon2plot,lat2plot,colour=scale))+
  scale_colour_gradientn(colours=rainbow(150), name="fitted scale")+
  theme_bw()

####### Winter 2000
df2plot <- data.frame("lon2plot"=spatial.sample$lon,
                      "lat2plot"=spatial.sample$lat,
                      "scale"=fitted.scale.winter)

library(ggplot2)                 
gg_winter <- ggplot() + ggtitle("Winter 2000") +
  geom_point(aes(lon,lat),colour="gray72", data = df, cex=.2)+
  geom_point(data = df2plot,aes(lon2plot,lat2plot,colour=scale))+
  scale_colour_gradientn(colours=rainbow(150), name="fitted scale")+
  theme_bw()

library(gridExtra)
grid.arrange(gg_summer,gg_winter, nrow=1,ncol=2)

###fitted spatio-temporal (with years and months) model with fixed shape
load("/Users/worklinda/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/mod.scale.st.year.months.Rdata")
####### January vs August
fitted.Gamma <- predict(mod.scale.st, type="response", newdata = data.frame("lat"=rep(spatial.sample$lat,2),
                                                                            "lon"=rep(spatial.sample$lon,2),
                                                                            "elev"=rep(sites_DE_sel$elevation[spatial.sample$id],2),
                                                                            "month"=rep(c(1,8), each=nrow(spatial.sample)),
                                                                            "year"=rep(2000, 2*nrow(spatial.sample))))
fitted.shape        <- rep(unique(fitted.Gamma[,2]), each=length(site2keep))
fitted.scale.summer <- fitted.Gamma[1:nrow(spatial.sample),1]
fitted.scale.winter <- fitted.Gamma[(nrow(spatial.sample)+1):(2*nrow(spatial.sample)),1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat)
####### January 2000
df2plot <- data.frame("lon2plot"=spatial.sample$lon,
                      "lat2plot"=spatial.sample$lat,
                      "scale"=fitted.scale.summer)

library(ggplot2)                 
gg_jan <- ggplot() + ggtitle("January 2000") +
  geom_point(aes(lon,lat),colour="gray72", data = df, cex=.2)+
  geom_point(data = df2plot,aes(lon2plot,lat2plot,colour=scale))+
  scale_colour_gradientn(colours=rainbow(150), name="fitted scale")+
  theme_bw()

####### August 2000
df2plot <- data.frame("lon2plot"=spatial.sample$lon,
                      "lat2plot"=spatial.sample$lat,
                      "scale"=fitted.scale.winter)

library(ggplot2)                 
gg_aug <- ggplot() + ggtitle("August 2000") +
  geom_point(aes(lon,lat),colour="gray72", data = df, cex=.2)+
  geom_point(data = df2plot,aes(lon2plot,lat2plot,colour=scale))+
  scale_colour_gradientn(colours=rainbow(150), name="fitted scale")+
  theme_bw()

library(gridExtra)
grid.arrange(gg_jan,gg_aug, nrow=1,ncol=2)

##################################################################################################################################
###                                                     vary shape                                                             ###
##################################################################################################################################
## fully temporal
mod2             <- mgcv::gam(list(obs~s(month, bs="cc", k=12)+s(year),
                             ~s(month, bs="cc", k=12)+s(year)), data = dat,
                        subset= which(dat$site %in% site2keep),family=gammals, select=TRUE)


