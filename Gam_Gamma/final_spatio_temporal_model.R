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

## spatio-temporal model with year and months
mod.scale.st <- gam(list(obs~s(year)+s(month, bs="cc", k=12)+s(elev)+s(lat,lon),
                         ~1), data = dat,
                    subset= which(dat$site %in% site2keep),
                    family=gammals, select=TRUE)
save(mod.scale.st, file="mod.scale.st.year.months.Rdata")

##### save the shape and scale of the Gamma (as parametrised in R)
scale.gamma.mat <- shape.gamma.mat <- NULL

for(yy in 1:length(unique(dates_sel$year))){
  print(yy)
  dat2pred          <- data.frame("lat"=rep(sites_DE_sel$lat,12),
                                  "lon"=rep(sites_DE_sel$lon,12),
                                  "elev"=rep(sites_DE_sel$elevation,12),
                                  "month"=rep(1:12, each=nrow(sites_DE_sel)),
                                  "year"=rep(unique(dates_sel$year)[yy], 12*nrow(sites_DE_sel)))
  gamma_pred    <- predict(mod.scale.st, dat2pred, type = "response")
  gamma_shape   <- exp(gamma_pred[,2])
  gamma_mean    <- gamma_pred[,1]
  scale.gamma.mat       <- rbind(scale.gamma.mat, matrix(gamma_mean*gamma_shape, byrow=TRUE, ncol=length(sites_DE_sel$lon)))
  shape.gamma.mat       <- rbind(shape.gamma.mat, matrix(1/gamma_shape, byrow=TRUE, ncol=length(sites_DE_sel$lon)))
}

dates_mat <- cbind(rep(unique(dates_sel$year), each=12), rep(1:12, length(unique(dates_sel$year)))) 

fct.merge.scale <- function(ss){
  dates_fitted_scale <- data.frame("year"=dates_mat[,1],
                                   "month"=dates_mat[,2],
                                   "scale"=scale.gamma.mat[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_scale, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_scale)))
  return(merge.try$scale)
}

fct.merge.shape <- function(ss){
  dates_fitted_shape <- data.frame("year"=dates_mat[,1],
                                   "month"=dates_mat[,2],
                                   "shape"=shape.gamma.mat[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_shape, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_shape)))
  return(merge.try$shape)
}

#####
fitted.scale.mat.list <- lapply(1:ncol(data_sel), fct.merge.scale) 
fitted.scale.mat      <- do.call(cbind,fitted.scale.mat.list)

fitted.shape.mat.list <- lapply(1:ncol(data_sel), fct.merge.shape) 
fitted.shape.mat      <- do.call(cbind,fitted.shape.mat.list)

fitted.gamma.param <- list("scale"=fitted.scale.mat,
                           "shape"=fitted.shape.mat)
save(fitted.gamma.param, file="data/fitted.gamma.param.Rdata")

###fitted spatio-temporal (with years and months) model with fixed shape
load("/Users/worklinda/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/mod.scale.st.year.months.Rdata")
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
png("spatial_year_month_gamma_median.png", width = 12, height = 5,units = 'in', res = 200)

grid.arrange(gg_january,gg_august, nrow=1,ncol=2)

dev.off()

##########################
##### Plot of the yearly and monthly effects

library("viridis")

year2pred <- seq(min(dat$year),max(dat$year),length.out = 100)
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

pdf("year_month_effect_median.pdf", width = 7, height = 5)
par(mar=c(3,3.2,1.5,0.5),mgp=c(1.6,0.5,0),font.main=1.3,cex=1.3,cex.main=1)

col2use <- plasma(12)
plot(year2pred, yeffect[,1], type="l", ylab="Estimated prob_0", xlab="Time", col=col2use[1], 
     ylim=c(range(yeffect)[1],range(yeffect)[2]))
text(x=year2pred[2],y=yeffect[2,1], labels = as.character(1), cex=0.6)
for(mm in 2:12){
  lines(year2pred, yeffect[,mm], col=col2use[mm])
  text(x=year2pred[2],y=yeffect[2,mm], labels = as.character(mm), cex=0.6)
}

dev.off()
