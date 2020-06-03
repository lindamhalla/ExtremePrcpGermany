
### Probability of obseeding the fixed threshold
obs.gamma     <- c(data_sel) #stacked site by site, i.e., obs of site 1, obs of site 2, etc...
dat.obs.binom <- data.frame("obs"=as.numeric(obs.gamma>0),
                            "lat"=rep(sites_DE_sel$lat, each=nrow(data_sel)),
                            "lon"=rep(sites_DE_sel$lon, each=nrow(data_sel)),
                            "elev"=rep(sites_DE_sel$elevation, each=nrow(data_sel)),
                            "year"=rep(dates_sel$year, ncol(data_sel)),
                            "month"=rep(dates_sel$month, ncol(data_sel)))

library(mgcv)
library(parallel)

cl <- makeCluster(detectCores()-2)
m_obs_binomial <- bam(obs~s(year)+s(month, bs="cc", k=12)+s(elev)+s(lat,lon),
                  family="binomial", data=dat.obs.binom, select=TRUE,
                  cluster=cl)
stopCluster(cl)

save(m_obs_binomial, file="data/m_obs_binomial.Rdata")

pr.obs.full.mat <- pr.obs.sd.full.mat <- NULL

for(yy in 1:length(unique(dates_sel$year))){
  print(yy)
  dat2pred          <- data.frame("lat"=rep(sites_DE_sel$lat,12),
                                  "lon"=rep(sites_DE_sel$lon,12),
                                  "elev"=rep(sites_DE_sel$elevation,12),
                                  "month"=rep(1:12, each=nrow(sites_DE_sel)),
                                  "year"=rep(unique(dates_sel$year)[yy], 12*nrow(sites_DE_sel)))
  pr.obs_pred        <- predict(m_obs_binomial, dat2pred, type = "response", se.fit = TRUE)
  pr.obs.full.mat    <- rbind(pr.obs.full.mat, matrix(pr.obs_pred$fit, byrow=TRUE, ncol=length(sites_DE_sel$lon)))
  pr.obs.sd.full.mat <- rbind(pr.obs.sd.full.mat, matrix(pr.obs_pred$se.fit, byrow=TRUE, ncol=length(sites_DE_sel$lon)))
}

dates_mat <- cbind(rep(unique(dates_sel$year), each=12), rep(1:12, length(unique(dates_sel$year)))) 

fct.merge.pr.obs <- function(ss){
  dates_fitted_pr.obs <- data.frame("year"=dates_mat[,1],
                                    "month"=dates_mat[,2],
                                    "pr.obs"=pr.obs.full.mat[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_pr.obs, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_pr.obs)))
  return(merge.try$pr.obs)
}

fct.merge.pr.obs.sd <- function(ss){
  dates_fitted_pr.obs <- data.frame("year"=dates_mat[,1],
                                    "month"=dates_mat[,2],
                                    "pr.obs.sd"=pr.obs.sd.full.mat[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_pr.obs, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_pr.obs)))
  return(merge.try$pr.obs.sd)
}


#####
fitted.pr.obs.mat.list <- lapply(1:ncol(data_sel), fct.merge.pr.obs) 
fitted.pr.obs.mat      <- do.call(cbind,fitted.pr.obs.mat.list)

fitted.pr.obs.sd.mat.list <- lapply(1:ncol(data_sel), fct.merge.pr.obs.sd) 
fitted.pr.obs.sd.mat      <- do.call(cbind,fitted.pr.obs.sd.mat.list)

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
png("fitted_pr_obs_jan2000.png", width = 15, height = 12,units = 'in', res = 200)

grid.arrange(gg.pr.obs,gg.pr.obs.se, nrow=1,ncol=2)

dev.off()

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
png("fitted_pr_obs_spatial.png", width = 12, height = 5,units = 'in', res = 200)

grid.arrange(gg_january,gg_august, nrow=1,ncol=2)

dev.off()

##########################
##### Plot of the yearly and monthly effects

library("viridis")

year2pred <- seq(min(dat.obs.binom$year),max(dat.obs.binom$year),length.out = 100)
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


pdf("year_month_effect_pr_obs.pdf", width = 7, height = 5)
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

