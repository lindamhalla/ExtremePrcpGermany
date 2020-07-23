library(fields)
library(INLA)
library(sp)
library(stringr)

load("data/german-precipitation.RData")
load("data/sites_DE_sel_elev.Rdata")
load("data/fitted_quant_90.Rdata")

### remove Zugspitze
data_sel     <- data_sel[,-13]
sites_DE_sel <- sites_DE_sel[-13,]

### Threshold exceedances
exc.gamma     <- c(data_sel-fitted.quant.mat) #stacked site by site, i.e., exc of site 1, exc of site 2, etc...
dat.exc.gamma <- data.frame("exc"=exc.gamma,
                            "lat"=rep(sites_DE_sel$lat, each=nrow(data_sel)),
                            "lon"=rep(sites_DE_sel$lon, each=nrow(data_sel)),
                            "elev"=rep(sites_DE_sel$elevation, each=nrow(data_sel)),
                            "year"=rep(dates_sel$year, ncol(data_sel)),
                            "month"=rep(dates_sel$month, ncol(data_sel)))

library(evgam)

dat.exc.gamma_gpd <- dat.exc.gamma
dat.exc.gamma_gpd$exc[dat.exc.gamma_gpd$exc <= 0] <- NA

fmla_gpd           <- list(exc ~ s(year)+s(month, bs="cc", k=12)+s(elev)+te(lat,lon), 
                           ~ s(year)+s(month, bs="cc", k=12)+s(elev)+te(lat,lon))
m_gpd_full_sub_exc <- evgam(fmla_gpd, dat.exc.gamma_gpd, family = "gpd")

save(m_gpd_full_sub_exc, file="data/m_gpd_full_ss_sub_exc.Rdata")

fitted.scale.mat <- fitted.shape.mat <- fitted.scale.sd.mat <- fitted.shape.sd.mat <- NULL

for(s in 1:nrow(sites_DE_sel)){
  gpd_lat_elev_pred <- predict(m_gpd_full_sub_exc, type="response", newdata = data.frame("lat"=rep(sites_DE_sel$lat[s],nrow(dates_sel)),
                                                                                         "lon"=rep(sites_DE_sel$lon[s],nrow(dates_sel)),
                                                                                         "elev"=rep(sites_DE_sel$elevation[s],nrow(dates_sel)),
                                                                                         "month"=dates_sel$month,
                                                                                         "year"=dates_sel$year),
                               se.fit=TRUE)
  gpd.fitted           <- gpd_lat_elev_pred$fitted
  gpd.fitted.se        <- gpd_lat_elev_pred$se.fit
  
  fitted.scale.mat <- cbind(fitted.scale.mat, gpd.fitted[,1])
  fitted.shape.mat <- cbind(fitted.shape.mat, gpd.fitted[,2])
  fitted.scale.sd.mat <- cbind(fitted.scale.sd.mat, gpd.fitted.se[,1])
  fitted.shape.sd.mat <- cbind(fitted.shape.sd.mat, gpd.fitted.se[,2])
}


fitted.gpd.full <- list("scale"=fitted.scale.mat,
                        "shape"=fitted.shape.mat,
                        "scale.sd"=fitted.scale.sd.mat,
                        "shape.sd"=fitted.shape.sd.mat)
save(fitted.gpd.full, file="data/fitted.gpd.full_ss.Rdata")

##########################################################################################
##### Plots of the GPD parameters on January 2000
##########################################################################################
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

grid.arrange(gg.scale,gg.scale.se,gg.shape,gg.shape.se, nrow=2,ncol=2)

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

grid.arrange(gg.scale,gg.scale.se,gg.shape,gg.shape.se, nrow=2,ncol=2)

##########################################################################################################################
##########################################################################################################################
df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "scale.jan"=fitted.gpd.full$scale[which((dates_sel$year==2000)&(dates_sel$month==1))[1],],
                 "scale.jan.se"=fitted.gpd.full$scale.sd[which((dates_sel$year==2000)&(dates_sel$month==1))[1],],
                 "shape.jan"=fitted.gpd.full$shape[which((dates_sel$year==2000)&(dates_sel$month==1))[1],],
                 "shape.jan.se"=fitted.gpd.full$shape.sd[which((dates_sel$year==2000)&(dates_sel$month==1))[1],],
                 "scale.aug"=fitted.gpd.full$scale[which((dates_sel$year==2000)&(dates_sel$month==8))[1],],
                 "scale.aug.se"=fitted.gpd.full$scale.sd[which((dates_sel$year==2000)&(dates_sel$month==8))[1],],
                 "shape.aug"=fitted.gpd.full$shape[which((dates_sel$year==2000)&(dates_sel$month==8))[1],],
                 "shape.aug.se"=fitted.gpd.full$shape.sd[which((dates_sel$year==2000)&(dates_sel$month==8))[1],])


library(ggplot2)   
############## Scale parameter
####### january 2000
gg_january <- ggplot() + ggtitle("january 2000") +
  geom_point(aes(lon,lat,colour=scale.jan), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Fitted scale", limits=c(min(c(df$scale.jan,df$scale.aug)),
                                                                             max(c(df$scale.jan,df$scale.aug))))+
  theme_bw()

####### august 2000
gg_august <- ggplot() + ggtitle("august 2000") +
  geom_point(aes(lon,lat,colour=scale.aug), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Fitted scale", limits=c(min(c(df$scale.jan,df$scale.aug)),
                                                                             max(c(df$scale.jan,df$scale.aug))))+
  theme_bw()

gg_january <- gg_january+ theme(legend.position = "none")
gg_august <- gg_august+ theme(legend.position = "bottom")

library(gridExtra)

grid.arrange(gg_january,gg_august, nrow=1,ncol=2)

############## Shape parameter
####### january 2000
gg_january <- ggplot() + ggtitle("january 2000") +
  geom_point(aes(lon,lat,colour=shape.jan), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Fitted shape", limits=c(min(c(df$shape.jan,df$shape.aug)),
                                                                             max(c(df$shape.jan,df$shape.aug))))+
  theme_bw()

####### august 2000
gg_august <- ggplot() + ggtitle("august 2000") +
  geom_point(aes(lon,lat,colour=shape.aug), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Fitted shape", limits=c(min(c(df$shape.jan,df$shape.aug)),
                                                                             max(c(df$shape.jan,df$shape.aug))))+
  theme_bw()

gg_january <- gg_january+ theme(legend.position = "none")
gg_august <- gg_august+ theme(legend.position = "bottom")

library(gridExtra)

grid.arrange(gg_january,gg_august, nrow=1,ncol=2)

##########################
##### Plot of the yearly and monthly effects

library("viridis")

year2pred  <- seq(min(dates_sel$year),max(dates_sel$year),length.out = 100)
yeffect.sc <- yeffect.sh <- NULL
for(mm in 1:12){
  bla    <- predict(m_gpd_full_sub_exc, newdata=data.frame("month"=mm,
                                                           "year"=year2pred,
                                                           "lat"=sites_DE_sel$lat[22],
                                                           "lon"=sites_DE_sel$lon[22],
                                                           "elev"=sites_DE_sel$elev[22]),
                    type = "response", se.fit = TRUE)
  
  yeffect.sc <- cbind(yeffect.sc, bla$fitted[,1])
  yeffect.sh <- cbind(yeffect.sh, bla$fitted[,2])
}

# year_month_effect_gpd_scale
col2use <- plasma(12)
plot(year2pred, yeffect.sc[,1], type="l", ylab="Fitted scale", xlab="Time", col=col2use[1], 
     ylim=c(range(yeffect.sc)[1],range(yeffect.sc)[2]))
text(x=year2pred[2],y=yeffect.sc[2,1], labels = as.character(1), cex=0.6)
for(mm in 2:12){
  lines(year2pred, yeffect.sc[,mm], col=col2use[mm])
  text(x=year2pred[2],y=yeffect.sc[2,mm], labels = as.character(mm), cex=0.6)
}

# year_month_effect_gpd_shape
col2use <- plasma(12)
plot(year2pred, yeffect.sh[,1], type="l", ylab="Fitted shape", xlab="Time", col=col2use[1], 
     ylim=c(range(yeffect.sh)[1],range(yeffect.sh)[2]))
text(x=year2pred[2],y=yeffect.sh[2,1], labels = as.character(1), cex=0.6)
for(mm in 2:12){
  lines(year2pred, yeffect.sh[,mm], col=col2use[mm])
  text(x=year2pred[2],y=yeffect.sh[2,mm], labels = as.character(mm), cex=0.6)
}

