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

### 14030*4526 matrix
anom.training.positive <- data_sel

# replace zeros with NAs
# replace observations >=threshold with NAs
# replace all 0<Z<u with Z/u \in (0,1)
for(i in 1:ncol(anom.training.positive)){
  anom.training.positive[(anom.training.positive[,i] == 0)|(anom.training.positive[,i] >= fitted.quant.mat[,i]),i] <- NA
}
# anom.training.positive[(anom.training.positive == 0)|(anom.training.positive >= fitted.quant.mat)] <- NA
# anom.training.positive <- anom.training.positive[!is.na(anom.training.positive)]

#gam with latitude and longitude in both the location and scale
dat   <- data.frame("obs"=as.vector(anom.training.positive)/as.vector(fitted.quant.mat),
                    "site"=rep(1:ncol(anom.training.positive), each= nrow(anom.training.positive)),
                    "lat"=rep(sites_DE_sel$lat, each= nrow(anom.training.positive)),
                    "lon"=rep(sites_DE_sel$lon, each= nrow(anom.training.positive)),
                    "elev"=rep(sites_DE_sel$elevation, each= nrow(anom.training.positive)),
                    "month"=rep(dates_sel$month, ncol(anom.training.positive)),
                    "year"=rep(dates_sel$year, ncol(anom.training.positive)))

# # Keep one site out of 40: selection made randomly
# set.seed(22)
# site2keep        <- sample(1:ncol(anom.training.positive), round(ncol(anom.training.positive)/40))
# 
#"homogeneous" subsampling from the map
load("data/SpatialSample.Rdata")
site2keep <- spatial.sample$id
site2keep[site2keep>13] <- site2keep[site2keep>13]-1
#########################################################################################################################
###                                       covariate-dependent Beta mean                                               ###
#########################################################################################################################
library(mgcv)
library(parallel)

cl <- makeCluster(detectCores()-1)

Sys.time()
## spatio-temporal model with year and months
mod.beta.st <- bam(obs~s(year)+s(month, bs="cc", k=12)+s(elev)+te(lat,lon), data = dat,
                    family=betar(link="logit"), select=TRUE,
                    # subset= which(dat$site %in% site2keep),
                    cluster=cl)
Sys.time()
stopCluster(cl)

save(mod.beta.st, file="data/mod.beta.st.year.months.Rdata")


##### save the shapes of the beta (as parametrised in R)
fitted.shape1.mat <- fitted.shape2.mat <- NULL

for(s in 1:nrow(sites_DE_sel)){
  print(s)
  beta_mean <- predict(mod.beta.st, type="response", newdata = data.frame("lat"=rep(sites_DE_sel$lat[s],nrow(dates_sel)),
                                                                          "lon"=rep(sites_DE_sel$lon[s],nrow(dates_sel)),
                                                                          "elev"=rep(sites_DE_sel$elevation[s],nrow(dates_sel)),
                                                                          "month"=dates_sel$month,
                                                                          "year"=dates_sel$year))
  beta_phi     <- mod.beta.st$family$getTheta(trans=TRUE)
  beta_shape1  <- beta_mean*beta_phi
  beta_shape2  <- beta_phi*(1-beta_mean)
  
  fitted.shape1.mat <- cbind(fitted.shape1.mat, beta_shape1)
  fitted.shape2.mat <- cbind(fitted.shape2.mat, beta_shape2)
}

fitted.beta.param <- list("shape1"=fitted.shape1.mat,
                          "shape2"=fitted.shape2.mat)
save(fitted.beta.param, file="data/fitted.beta.param.Rdata")

###fitted spatio-temporal (with years and months) model with fixed shape
load("data/mod.beta.st.year.months.Rdata")
####### January vs August
fitted.beta <- predict(mod.beta.st, type="response", newdata = data.frame("lat"=rep(sites_DE_sel$lat,2),
                                                                            "lon"=rep(sites_DE_sel$lon,2),
                                                                            "elev"=rep(sites_DE_sel$elevation,2),
                                                                            "month"=rep(c(1,8), 
                                                                                        each=nrow(sites_DE_sel)),
                                                                            "year"=rep(2000, 2*nrow(sites_DE_sel))))
fitted.theta        <- mod.beta.st$family$getTheta(trans=TRUE)
fitted.mean.january <- fitted.beta[1:nrow(sites_DE_sel)]
fitted.mean.august  <- fitted.beta[(nrow(sites_DE_sel)+1):(2*nrow(sites_DE_sel))]
fitted.median.january <- qbeta(0.5, shape1=fitted.mean.january*fitted.theta,
                                shape2=(1-fitted.mean.january)*fitted.theta)*
  apply(fitted.quant.mat[which((dates_sel$month==1)&(dates_sel$year==2000)),],2,unique)
fitted.median.august  <- qbeta(0.5, shape1=fitted.mean.august*fitted.theta,
                               shape2=(1-fitted.mean.august)*fitted.theta)*
  apply(fitted.quant.mat[which((dates_sel$month==8)&(dates_sel$year==2000)),],2,unique)

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "median.january"=fitted.median.january,
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
png("spatial_year_month_beta_median.png", width = 12, height = 5,units = 'in', res = 200)

grid.arrange(gg_january,gg_august, nrow=1,ncol=2)

dev.off()

##########################
##### Plot of the yearly and monthly effects

library("viridis")

year2pred <- seq(min(dat$year),max(dat$year)-1,by=1)
yeffect   <- yeffect.sd <- NULL
for(mm in 1:12){
  bla     <- predict(mod.beta.st, type="response", newdata=data.frame("month"=rep(mm,length(year2pred)),
                                                                       "year"=year2pred,
                                                                       "lat"=rep(sites_DE_sel$lat[22],length(year2pred)),
                                                                       "lon"=rep(sites_DE_sel$lon[22],length(year2pred)),
                                                                       "elev"=rep(sites_DE_sel$elev[22],length(year2pred))))
  
  fitted.shape1 <- bla*mod.beta.st$family$getTheta(trans=TRUE)
  fitted.shape2 <- (1-bla)*mod.beta.st$family$getTheta(trans=TRUE)
  quant.site    <- NULL
  for(j in year2pred){
    quant.site <- c(quant.site, 
                    unique(fitted.quant.mat[which((dates_sel$month==mm)&(dates_sel$year ==j)),22])[1])
    #print(c(j,unique(fitted.quant.mat[which((dates_sel$month==mm)&(dates_sel$year ==j)),22])))
  }
  
  # yeffect <- cbind(yeffect, qbeta(0.5, shape1=fitted.shape1, shape2=fitted.shape2)*
  #                    unique(fitted.quant.mat[which((dates_sel$month==mm)&(dates_sel$year %in% year2pred)),22]))
  yeffect <- cbind(yeffect, qbeta(0.5, shape1=fitted.shape1, shape2=fitted.shape2)*
                     quant.site)
}

pdf("year_month_effect_median_beta.pdf", width = 7, height = 5)
par(mar=c(3,3.2,1.5,0.5),mgp=c(1.6,0.5,0),font.main=1.3,cex=1.3,cex.main=1)

col2use <- plasma(12)
plot(year2pred, yeffect[,1], type="l", ylab="Estimated median", xlab="Time", col=col2use[1], 
     ylim=c(range(yeffect)[1],range(yeffect)[2]))
text(x=year2pred[2],y=yeffect[2,1], labels = as.character(1), cex=0.6)
for(mm in 2:12){
  lines(year2pred, yeffect[,mm], col=col2use[mm])
  text(x=year2pred[2],y=yeffect[2,mm], labels = as.character(mm), cex=0.6)
}

dev.off()
