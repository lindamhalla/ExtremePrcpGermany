library(fields)
library(INLA)
library(sp)
library(stringr)

load("data/german-precipitation.RData")
load("data/sites_DE_sel_elev.Rdata")

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

#"homogeneous" subsampling from the map
load("data/SpatialSample.Rdata")
site2keep <- spatial.sample$id
##################################################################################################################################
###                                                    fixed shape                                                             ###
##################################################################################################################################

## spatio-temporal model with year and months
mod.scale.st <- gam(list(obs~s(year)+s(month, bs="cc", k=12)+s(elev)+te(lat,lon),
                         ~1), data = dat,
                    subset= which(dat$site %in% site2keep),
                    family=gammals, select=TRUE)
save(mod.scale.st, file="data/mod.scale.st.year.months.Rdata")

### remove Zugspitze
data_sel     <- data_sel[,-13]
sites_DE_sel <- sites_DE_sel[-13,]

##### save the shape and scale of the Gamma (as parametrised in R)
fitted.scale.mat <- fitted.shape.mat <- NULL

for(s in 1:nrow(sites_DE_sel)){
  gamma_pred  <- predict(mod.scale.st, type="response", newdata = data.frame("lat"=rep(sites_DE_sel$lat[s],nrow(dates_sel)),
                                                                             "lon"=rep(sites_DE_sel$lon[s],nrow(dates_sel)),
                                                                             "elev"=rep(sites_DE_sel$elevation[s],nrow(dates_sel)),
                                                                             "month"=dates_sel$month,
                                                                             "year"=dates_sel$year))
  
  gamma_shape   <- exp(gamma_pred[,2])
  gamma_mean    <- gamma_pred[,1]
  
  fitted.scale.mat <- cbind(fitted.scale.mat, gamma_mean*gamma_shape)
  fitted.shape.mat <- cbind(fitted.shape.mat, 1/gamma_shape)
}

fitted.gamma.param <- list("scale"=fitted.scale.mat,
                           "shape"=fitted.shape.mat)
save(fitted.gamma.param, file="data/fitted.gamma.param.Rdata")

###fitted spatio-temporal (with years and months) model with fixed shape
load("data/mod.scale.st.year.months.Rdata")
####### January vs August
fitted.Gamma <- predict(mod.scale.st, type="response", newdata = data.frame("lat"=rep(sites_DE_sel$lat,2),
                                                                            "lon"=rep(sites_DE_sel$lon,2),
                                                                            "elev"=rep(sites_DE_sel$elevation,2),
                                                                            "month"=rep(c(1,8), 
                                                                                        each=nrow(sites_DE_sel)),
                                                                            "year"=rep(2000, 2*nrow(sites_DE_sel))))
fitted.shape          <- fitted.Gamma[1:nrow(sites_DE_sel),2] #same shape for all (january and august)
fitted.scale.january  <- fitted.Gamma[1:nrow(sites_DE_sel),1]
fitted.scale.august   <- fitted.Gamma[(nrow(sites_DE_sel)+1):(2*nrow(sites_DE_sel)),1]
fitted.median.january <- qgamma(0.5, shape=1/exp(fitted.shape), scale=fitted.scale.january*exp(fitted.shape))
fitted.median.august  <- qgamma(0.5, shape=1/exp(fitted.shape), scale=fitted.scale.august*exp(fitted.shape))

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "scale.january"=fitted.scale.january,
                 "median.january"=fitted.median.january,
                 "scale.august"=fitted.scale.august,
                 "median.august"=fitted.median.august)

library(ggplot2)   
####### january 2000
gg_january <- ggplot() + ggtitle("january 2000") +
  geom_point(aes(lon,lat,colour=median.january), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated median", limits=c(min(c(df$median.january,df$median.august)),
                                                                                 max(c(df$median.january,df$median.august))))+
  theme_bw()

####### august 2000
gg_august <- ggplot() + ggtitle("august 2000") +
  geom_point(aes(lon,lat,colour=median.august), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated median", limits=c(min(c(df$median.january,df$median.august)),
                                                                                 max(c(df$median.january,df$median.august))))+
  theme_bw()

gg_january <- gg_january+ theme(legend.position = "none")
gg_august <- gg_august+ theme(legend.position = "bottom")

library(gridExtra)

grid.arrange(gg_january,gg_august, nrow=1,ncol=2)

##########################
##### Plot of the yearly and monthly effects

library("viridis")

year2pred <- seq(min(dates_sel$year),max(dates_sel$year),length.out = 100)
yeffect   <- yeffect.sd <- NULL
for(mm in 1:12){
  bla     <- predict(mod.scale.st, type="response", newdata=data.frame("month"=mm,
                                                                         "year"=year2pred,
                                                                         "lat"=sites_DE_sel$lat[22],
                                                                         "lon"=sites_DE_sel$lon[22],
                                                                         "elev"=sites_DE_sel$elev[22]))
  
  fitted.shape <- bla[,2]
  fitted.mean  <- bla[,1]
  
  yeffect <- cbind(yeffect, qgamma(0.5, shape=1/exp(fitted.shape), scale=fitted.mean*exp(fitted.shape)))
}

col2use <- plasma(12)
plot(year2pred, yeffect[,1], type="l", ylab="Estimated median", xlab="Time", col=col2use[1], 
     ylim=c(range(yeffect)[1],range(yeffect)[2]))
text(x=year2pred[2],y=yeffect[2,1], labels = as.character(1), cex=0.6)
for(mm in 2:12){
  lines(year2pred, yeffect[,mm], col=col2use[mm])
  text(x=year2pred[2],y=yeffect[2,mm], labels = as.character(mm), cex=0.6)
}
