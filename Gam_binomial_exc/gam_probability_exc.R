### model the probability of exceeding the fixed threshold (the 90% quantile of the fitted Gamma distribution)

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

### Probability of exceeding the fixed threshold
exc.gamma     <- c(data_sel-fitted.quant.mat) #stacked site by site, i.e., exc of site 1, exc of site 2, etc...
dat.exc.binom <- data.frame("exc"=as.numeric(exc.gamma>0),
                            "lat"=rep(sites_DE_sel$lat, each=nrow(data_sel)),
                            "lon"=rep(sites_DE_sel$lon, each=nrow(data_sel)),
                            "elev"=rep(sites_DE_sel$elevation, each=nrow(data_sel)),
                            "year"=rep(dates_sel$year, ncol(data_sel)),
                            "month"=rep(dates_sel$month, ncol(data_sel)))

library(mgcv)
library(parallel)

cl <- makeCluster(detectCores()-1)
m_binomial <- bam(exc~s(year)+s(month, bs="cc", k=12)+s(elev)+te(lat,lon),
                  family="binomial", data=dat.exc.binom, select=TRUE,
                  cluster=cl)
stopCluster(cl)

save(m_binomial, file="data/m_binomial.Rdata")

######### Obtain matrix of probabilities of exceedance
fitted.pr.exc.mat <- fitted.pr.exc.sd.mat <- NULL

for(s in 1:nrow(sites_DE_sel)){
  pr.exc_pred  <- predict(m_binomial, type="response", newdata = data.frame("lat"=rep(sites_DE_sel$lat[s],nrow(dates_sel)),
                                                                              "lon"=rep(sites_DE_sel$lon[s],nrow(dates_sel)),
                                                                              "elev"=rep(sites_DE_sel$elevation[s],nrow(dates_sel)),
                                                                              "month"=dates_sel$month,
                                                                              "year"=dates_sel$year),
                          se.fit=TRUE)

  fitted.pr.exc.mat    <- cbind(fitted.pr.exc.mat, pr.exc_pred$fit)
  fitted.pr.exc.sd.mat <- cbind(fitted.pr.exc.sd.mat, pr.exc_pred$se.fit)
}

fitted.pr.exc <- list("fitted.pr.exc.mat"=fitted.pr.exc.mat,
                      "fitted.pr.exc.sd.mat"=fitted.pr.exc.sd.mat)
save(fitted.pr.exc, file="data/fitted.pr.exc.Rdata")

##### Plots of the prob of positive precip on January 2000
which.idx <- which((dates_sel$year==2000)&(dates_sel$month==1))[1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "pr.exc"=fitted.pr.exc.mat[which.idx,],
                 "pr.exc.se"=fitted.pr.exc.sd.mat[which.idx,])

library(ggplot2)                 
gg.pr.exc <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.exc), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(9), name="fitted prob_exc")+
  theme_bw()

gg.pr.exc.se <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.exc.se), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(9), name="se prob_exc")+
  theme_bw()

library(gridExtra)

grid.arrange(gg.pr.exc,gg.pr.exc.se, nrow=1,ncol=2)

##### Plots of the prob of positive precip on August 2000
which.idx <- which((dates_sel$year==2000)&(dates_sel$month==8))[1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "pr.exc"=fitted.pr.exc.mat[which.idx,],
                 "pr.exc.se"=fitted.pr.exc.sd.mat[which.idx,])

library(ggplot2)                 
gg.pr.exc.aug <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.exc), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(9), name="fitted prob_exc")+
  theme_bw()

gg.pr.exc.se.aug <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.exc.se), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(9), name="se prob_exc")+
  theme_bw()

grid.arrange(gg.pr.exc.aug,gg.pr.exc.se.aug, nrow=1,ncol=2)

#######
df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "pr.january"=100*fitted.pr.exc.mat[which((dates_sel$year==2000)&(dates_sel$month==1))[1],],
                 "pr.august"=100*fitted.pr.exc.mat[which((dates_sel$year==2000)&(dates_sel$month==8))[1],])

library(ggplot2)   
####### january 2000
gg_january <- ggplot() + ggtitle("january 2000") +
  geom_point(aes(lon,lat,colour=pr.january), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(9), name="Estimated prob_exc (%)", limits=c(min(c(df$pr.january,df$pr.august)),
                                                                                 max(c(df$pr.january,df$pr.august))))+
  theme_bw()

####### august 2000
gg_august <- ggplot() + ggtitle("august 2000") +
  geom_point(aes(lon,lat,colour=pr.august), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(9), name="Estimated prob_exc (%)", limits=c(min(c(df$pr.january,df$pr.august)),
                                                                                 max(c(df$pr.january,df$pr.august))))+
  theme_bw()

gg_january <- gg_january+ theme(legend.position = "none")
gg_august <- gg_august+ theme(legend.position = "bottom")

library(gridExtra)

grid.arrange(gg_january,gg_august, nrow=1,ncol=2)

##########################
##### Plot of the yearly and monthly effects

library("viridis")

year2pred <- seq(min(dat.exc.binom$year),max(dat.exc.binom$year),length.out = 100)
yeffect   <- yeffect.sd <- NULL
for(mm in 1:12){
  bla     <- predict(m_binomial, type="response", newdata=data.frame("month"=mm,
                                                                     "year"=year2pred,
                                                                     "lat"=sites_DE_sel$lat[22],
                                                                     "lon"=sites_DE_sel$lon[22],
                                                                     "elev"=sites_DE_sel$elev[22]), se.fit = TRUE)
  yeffect    <- cbind(yeffect, bla$fit)
  yeffect.sd <- cbind(yeffect.sd, bla$se.fit)
}

col2use <- plasma(12)
plot(year2pred, yeffect[,1], type="l", ylab="Estimated prob_exc", xlab="Time", col=col2use[1], 
     ylim=c(range(yeffect)[1],range(yeffect)[2]))
text(x=year2pred[2],y=yeffect[2,1], labels = as.character(1), cex=0.6)
for(mm in 2:12){
  lines(year2pred, yeffect[,mm], col=col2use[mm])
  text(x=year2pred[2],y=yeffect[2,mm], labels = as.character(mm), cex=0.6)
}
