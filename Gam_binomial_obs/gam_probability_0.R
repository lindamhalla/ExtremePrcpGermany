
### model the probability of positive precipitation

library(fields)
library(INLA)
library(sp)
library(stringr)

load("data/german-precipitation.RData")
load("data/sites_DE_sel_elev.Rdata")
load("data/fitted_quant_75.Rdata")

### remove Zugspitze
data_sel     <- data_sel[,-13]
sites_DE_sel <- sites_DE_sel[-13,]

### Probability of being positive
obs.gamma     <- c(data_sel) #stacked site by site, i.e., obs of site 1, obs of site 2, etc...
dat.obs.binom <- data.frame("obs"=as.numeric(obs.gamma>0),
                            "lat"=rep(sites_DE_sel$lat, each=nrow(data_sel)),
                            "lon"=rep(sites_DE_sel$lon, each=nrow(data_sel)),
                            "elev"=rep(sites_DE_sel$elevation, each=nrow(data_sel)),
                            "year"=rep(dates_sel$year, ncol(data_sel)),
                            "month"=rep(dates_sel$month, ncol(data_sel)))

library(mgcv)
library(parallel)

cl <- makeCluster(detectCores()-1)
m_obs_binomial <- bam(obs~s(year)+s(month, bs="cc", k=12)+s(elev)+te(lat,lon),
                  family="binomial", data=dat.obs.binom, select=TRUE,
                  cluster=cl)
stopCluster(cl)

save(m_obs_binomial, file="data/m_obs_binomial.Rdata")

### save fitted probabilities in a matrix
fitted.pr.obs.mat <- fitted.pr.obs.sd.mat <- NULL

for(s in 1:nrow(sites_DE_sel)){
  pr.exc_pred  <- predict(m_obs_binomial, type="response", newdata = data.frame("lat"=rep(sites_DE_sel$lat[s],nrow(dates_sel)),
                                                                                "lon"=rep(sites_DE_sel$lon[s],nrow(dates_sel)),
                                                                                "elev"=rep(sites_DE_sel$elevation[s],nrow(dates_sel)),
                                                                                "month"=dates_sel$month,
                                                                                "year"=dates_sel$year),
                          se.fit=TRUE)
  
  fitted.pr.obs.mat    <- cbind(fitted.pr.obs.mat, pr.exc_pred$fit)
  fitted.pr.obs.sd.mat <- cbind(fitted.pr.obs.sd.mat, pr.exc_pred$se.fit)
}


#####
fitted.pr.obs <- list("fitted.pr.obs.mat"=fitted.pr.obs.mat,
                      "fitted.pr.obs.sd.mat"=fitted.pr.obs.sd.mat)
save(fitted.pr.obs, file="data/fitted.pr.obs.Rdata")

##### Plots of the prob of positive precip on January 2000
which.idx <- which((dates_sel$year==2000)&(dates_sel$month==1))[1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "pr.obs"=fitted.pr.obs.mat[which.idx,],
                 "pr.obs.se"=fitted.pr.obs.sd.mat[which.idx,])

library(ggplot2)                 
gg.pr.obs <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.obs), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="fitted prob_0")+
  theme_bw()

gg.pr.obs.se <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.obs.se), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="se prob_0")+
  theme_bw()

library(gridExtra)

grid.arrange(gg.pr.obs,gg.pr.obs.se, nrow=1,ncol=2)

##### Plots of the prob of positive precip on August 2000
which.idx <- which((dates_sel$year==2000)&(dates_sel$month==8))[1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "pr.obs"=fitted.pr.obs.mat[which.idx,],
                 "pr.obs.se"=fitted.pr.obs.sd.mat[which.idx,])

library(ggplot2)                 
gg.pr.obs.aug <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.obs), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="fitted prob_0")+
  theme_bw()

gg.pr.obs.se.aug <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.obs.se), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="se prob_0")+
  theme_bw()

grid.arrange(gg.pr.obs.aug,gg.pr.obs.se.aug, nrow=1,ncol=2)

#######
df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "pr.january"=fitted.pr.obs.mat[which((dates_sel$year==2000)&(dates_sel$month==1))[1],],
                 "pr.august"=fitted.pr.obs.mat[which((dates_sel$year==2000)&(dates_sel$month==8))[1],])

library(ggplot2)   
####### january 2000
gg_january <- ggplot() + ggtitle("january 2000") +
  geom_point(aes(lon,lat,colour=pr.january), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated prob_0", limits=c(min(c(df$pr.january,df$pr.august)),
                                                                                 max(c(df$pr.january,df$pr.august))))+
  theme_bw()

####### august 2000
gg_august <- ggplot() + ggtitle("august 2000") +
  geom_point(aes(lon,lat,colour=pr.august), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated prob_0", limits=c(min(c(df$pr.january,df$pr.august)),
                                                                                 max(c(df$pr.january,df$pr.august))))+
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
  bla     <- predict(m_obs_binomial, type="response", newdata=data.frame("month"=mm,
                                                                         "year"=year2pred,
                                                                         "lat"=sites_DE_sel$lat[22],
                                                                         "lon"=sites_DE_sel$lon[22],
                                                                         "elev"=sites_DE_sel$elev[22]), se.fit = TRUE)
  yeffect    <- cbind(yeffect, bla$fit)
  yeffect.sd <- cbind(yeffect.sd, bla$se.fit)
}

col2use <- plasma(12)
plot(year2pred, yeffect[,1], type="l", ylab="Estimated prob_0", xlab="Time", col=col2use[1], 
     ylim=c(range(yeffect)[1],range(yeffect)[2]))
text(x=year2pred[2],y=yeffect[2,1], labels = as.character(1), cex=0.6)
for(mm in 2:12){
  lines(year2pred, yeffect[,mm], col=col2use[mm])
  text(x=year2pred[2],y=yeffect[2,mm], labels = as.character(mm), cex=0.6)
}
