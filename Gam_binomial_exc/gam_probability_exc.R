### model the probability of exceeding the fixed threshold (the 75% quantile of the fitted Gamma distribution)

library(fields)
library(INLA)
library(sp)
library(stringr)

load("data/german-precipitation.RData")
load("data/sites_DE_sel_elev.Rdata")
load("data/fitted_quant_75.Rdata")


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

cl <- makeCluster(detectCores()-2)
m_binomial <- bam(exc~s(year)+s(month, bs="cc", k=12)+s(elev)+s(lat,lon),
                  family="binomial", data=dat.exc.binom, select=TRUE,
                  cluster=cl)
stopCluster(cl)

save(m_binomial, file="data/m_binomial.Rdata")

######### Obtain matrix of probabilities of exceedance
pr.exc.full.mat <- pr.exc.sd.full.mat <- NULL

for(yy in 1:length(unique(dates_sel$year))){
  print(yy)
  dat2pred          <- data.frame("lat"=rep(sites_DE_sel$lat,12),
                                  "lon"=rep(sites_DE_sel$lon,12),
                                  "elev"=rep(sites_DE_sel$elevation,12),
                                  "month"=rep(1:12, each=nrow(sites_DE_sel)),
                                  "year"=rep(unique(dates_sel$year)[yy], 12*nrow(sites_DE_sel)))
  pr.exc_pred        <- predict(m_binomial, dat2pred, type = "response", se.fit = TRUE)
  pr.exc.full.mat    <- rbind(pr.exc.full.mat, matrix(pr.exc_pred$fit, byrow=TRUE, ncol=length(sites_DE_sel$lon)))
  pr.exc.sd.full.mat <- rbind(pr.exc.sd.full.mat, matrix(pr.exc_pred$se.fit, byrow=TRUE, ncol=length(sites_DE_sel$lon)))
}

dates_mat <- cbind(rep(unique(dates_sel$year), each=12), rep(1:12, length(unique(dates_sel$year)))) 

fct.merge.pr.exc <- function(ss){
  dates_fitted_pr.exc <- data.frame("year"=dates_mat[,1],
                                    "month"=dates_mat[,2],
                                    "pr.exc"=pr.exc.full.mat[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_pr.exc, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_pr.exc)))
  return(merge.try$pr.exc)
}

fct.merge.pr.exc.sd <- function(ss){
  dates_fitted_pr.exc <- data.frame("year"=dates_mat[,1],
                                    "month"=dates_mat[,2],
                                    "pr.exc.sd"=pr.exc.sd.full.mat[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_pr.exc, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_pr.exc)))
  return(merge.try$pr.exc.sd)
}


#####
fitted.pr.exc.mat.list <- lapply(1:ncol(data_sel), fct.merge.pr.exc) 
fitted.pr.exc.mat      <- do.call(cbind,fitted.pr.exc.mat.list)

fitted.pr.exc.sd.mat.list <- lapply(1:ncol(data_sel), fct.merge.pr.exc.sd) 
fitted.pr.exc.sd.mat      <- do.call(cbind,fitted.pr.exc.sd.mat.list)

##### Plots of the prob of positive precip on January 2000
which.idx <- which((dates_sel$year==2000)&(dates_sel$month==1))[1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "pr.exc"=fitted.pr.exc.mat[which.idx,],
                 "pr.exc.se"=fitted.pr.exc.sd.mat[which.idx,])

library(ggplot2)                 
gg.pr.exc <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.exc), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="fitted prob_exc")+
  theme_bw()

gg.pr.exc.se <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.exc.se), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="se prob_exc")+
  theme_bw()

library(gridExtra)
png("fitted_pr_exc_jan2000.png", width = 15, height = 12,units = 'in', res = 200)

grid.arrange(gg.pr.exc,gg.pr.exc.se, nrow=1,ncol=2)

dev.off()

##### Plots of the prob of positive precip on August 2000
which.idx <- which((dates_sel$year==2000)&(dates_sel$month==8))[1]

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "pr.exc"=fitted.pr.exc.mat[which.idx,],
                 "pr.exc.se"=fitted.pr.exc.sd.mat[which.idx,])

library(ggplot2)                 
gg.pr.exc.aug <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.exc), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="fitted prob_exc")+
  theme_bw()

gg.pr.exc.se.aug <- ggplot() +
  geom_point(aes(lon,lat,colour=pr.exc.se), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="se prob_exc")+
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
  scale_colour_gradientn(colours=rainbow(150), name="Estimated prob_exc (%)", limits=c(min(c(df$pr.january,df$pr.august)),
                                                                                 max(c(df$pr.january,df$pr.august))))+
  theme_bw()

####### august 2000
gg_august <- ggplot() + ggtitle("august 2000") +
  geom_point(aes(lon,lat,colour=pr.august), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated prob_exc (%)", limits=c(min(c(df$pr.january,df$pr.august)),
                                                                                 max(c(df$pr.january,df$pr.august))))+
  theme_bw()

gg_january <- gg_january+ theme(legend.position = "none")
gg_august <- gg_august+ theme(legend.position = "bottom")

library(gridExtra)
png("fitted_pr_exc_spatial.png", width = 12, height = 5,units = 'in', res = 200)

grid.arrange(gg_january,gg_august, nrow=1,ncol=2)

dev.off()

##########################
##### Plot of the yearly and monthly effects

library("viridis")

year2pred <- seq(min(dat.obs.binom$year),max(dat.obs.binom$year),length.out = 100)
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

pdf("year_month_effect_pr_exc.pdf", width = 7, height = 5)
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
