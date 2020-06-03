#########################################################################################################
#################### plots of estimated median
#########################################################################################################

#fully spatial
load("mod.scale.spatial.Rdata")
fitted.Gamma <- predict(mod.scale.spatial, type="response", newdata = data.frame("lat"=sites_DE_sel$lat,
                                                                                 "lon"=sites_DE_sel$lon,
                                                                                 "elev"=sites_DE_sel$elevation))


fitted.shape  <- exp(fitted.Gamma[,2])
fitted.mean   <- fitted.Gamma[,1]
fitted.median <- qgamma(0.5, shape=1/fitted.shape, scale=fitted.mean*fitted.shape)

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "scale"=fitted.mean*fitted.shape,
                 "median"=fitted.median)

library(ggplot2)                 
gg.spatial <- ggplot() +
  geom_point(aes(lon,lat,colour=median), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="estimated median")+
  theme_bw()

png("spatial_gamma.png", width = 12, height = 9,units = 'in', res = 200)

gg.spatial

dev.off()

#fully temporal
load("mod.scale.temporal.Rdata")
library("viridis")

year2pred <- seq(min(dat$year),max(dat$year),length.out = 100)
yeffect   <- NULL
for(mm in 1:12){
  yeffect <- cbind(yeffect, predict(mod.scale.temporal, type="response", newdata=data.frame("month"=mm, 
                                                                                            "year"=year2pred))[,1])
}

fitted.shape  <- predict(mod.scale.temporal, type="response", newdata=data.frame("month"=mm, "year"=year2pred))[,2]
fitted.median <- NULL
for(mm in 1:12){
  fitted.median <- cbind(fitted.median, qgamma(0.5, shape=1/exp(fitted.shape), scale=yeffect[,mm]*exp(fitted.shape)))
}

pdf("year_effect_median.pdf", width = 7, height = 5)
par(mar=c(3,3.2,1.5,0.5),mgp=c(1.6,0.5,0),font.main=1.3,cex=1.3,cex.main=1)

col2use <- plasma(12)
plot(year2pred, fitted.median[,1], type="l", ylab="Estimated median", xlab="Time", col=col2use[1], 
     ylim=c(range(fitted.median)[1],range(fitted.median)[2]))
text(x=year2pred[2],y=fitted.median[2,1], labels = as.character(1), cex=0.6)
for(mm in 2:12){
  lines(year2pred, fitted.median[,mm], col=col2use[mm])
  text(x=year2pred[2],y=fitted.median[2,mm], labels = as.character(mm), cex=0.6)
}

dev.off()

###fitted spatio-temporal (with years) model
load("mod.scale.st.Rdata")
fitted.Gamma <- predict(mod.scale.st, type="response", newdata = data.frame("lat"=sites_DE_sel$lat,
                                                                            "lon"=sites_DE_sel$lon,
                                                                            "elev"=sites_DE_sel$elevation,
                                                                            "year"=2000))
fitted.shape  <- rep(unique(fitted.Gamma[,2]), each=length(site2keep))
fitted.scale  <- fitted.Gamma[,1]
fitted.median <- qgamma(0.5, shape=fitted.shape, scale=fitted.scale)

df.2000 <- data.frame("lon"=sites_DE_sel$lon,
                      "lat"=sites_DE_sel$lat,
                      "scale"=fitted.scale,
                      "median"=fitted.median)

fitted.Gamma <- predict(mod.scale.st, type="response", newdata = data.frame("lat"=sites_DE_sel$lat,
                                                                            "lon"=sites_DE_sel$lon,
                                                                            "elev"=sites_DE_sel$elevation,
                                                                            "year"=2018))
fitted.shape  <- rep(unique(fitted.Gamma[,2]), each=length(site2keep))
fitted.scale  <- fitted.Gamma[,1]
fitted.median <- qgamma(0.5, shape=fitted.shape, scale=fitted.scale)

df.2018 <- data.frame("lon"=sites_DE_sel$lon,
                      "lat"=sites_DE_sel$lat,
                      "scale"=fitted.scale,
                      "median"=fitted.median)

library(ggplot2)                 
gg2000 <- ggplot() + ggtitle("2000") +
  geom_point(aes(lon,lat,colour=median), data = df.2000, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated median", limits=c(min(c(df.2000$median,df.2018$median)),
                                                                                 max(c(df.2000$median,df.2018$median))))+
  theme_bw()

library(ggplot2)                 
gg2018 <- ggplot() + ggtitle("2018") +
  geom_point(aes(lon,lat,colour=median), data = df.2018, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated median", limits=c(min(c(df.2000$median,df.2018$median)),
                                                                                 max(c(df.2000$median,df.2018$median))))+
  theme_bw()

gg2000 <- gg2000+ theme(legend.position = "none")
gg2018 <- gg2018+ theme(legend.position = "bottom")

png("spatial_years_gamma_median.png", width = 12, height = 5,units = 'in', res = 200)

grid.arrange(gg2000,gg2018,nrow=1, ncol=2)

dev.off()

###fitted spatio-temporal (with months) model with fixed shape
load("mod.scale.st.month.Rdata")
#######January
fitted.Gamma <- predict(mod.scale.st.month, type="response", newdata = data.frame("lat"=sites_DE_sel$lat,
                                                                                  "lon"=sites_DE_sel$lon,
                                                                                  "elev"=sites_DE_sel$elevation,
                                                                                  "month"=1))
fitted.shape  <- rep(unique(fitted.Gamma[,2]), each=length(site2keep))
fitted.scale  <- fitted.Gamma[,1]
fitted.median <- qgamma(0.5, shape=fitted.shape, scale=fitted.scale)

df.jan <- data.frame("lon"=sites_DE_sel$lon,
                     "lat"=sites_DE_sel$lat,
                     "scale"=fitted.scale,
                     "median"=fitted.median)

#######August
fitted.Gamma <- predict(mod.scale.st.month, type="response", newdata = data.frame("lat"=sites_DE_sel$lat,
                                                                                  "lon"=sites_DE_sel$lon,
                                                                                  "elev"=sites_DE_sel$elevation,
                                                                                  "month"=8))
fitted.shape  <- rep(unique(fitted.Gamma[,2]), each=length(site2keep))
fitted.scale  <- fitted.Gamma[,1]
fitted.median <- qgamma(0.5, shape=fitted.shape, scale=fitted.scale)

df.aug <- data.frame("lon"=sites_DE_sel$lon,
                     "lat"=sites_DE_sel$lat,
                     "scale"=fitted.scale,
                     "median"=fitted.median)

library(ggplot2)                 
gg_jan <- ggplot() + ggtitle("January") +
  geom_point(aes(lon,lat,colour=median), data = df.jan, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated median", limits=c(min(c(df.jan$median,df.aug$median)),
                                                                                 max(c(df.jan$median,df.aug$median))))+
  theme_bw()

gg_aug <- ggplot() + ggtitle("August") +
  geom_point(aes(lon,lat,colour=median), data = df.aug, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated median", limits=c(min(c(df.jan$median,df.aug$median)),
                                                                                 max(c(df.jan$median,df.aug$median))))+
  theme_bw()

gg_jan <- gg_jan+ theme(legend.position = "none")
gg_aug <- gg_aug+ theme(legend.position = "bottom")

library(gridExtra)
png("spatial_months_gamma_median.png", width = 12, height = 5,units = 'in', res = 200)

grid.arrange(gg_jan,gg_aug, nrow=1,ncol=2)

dev.off()

###fitted spatio-temporal (with "seasons") model with fixed shape
load("mod.scale.st.season.additive.Rdata")
####### Summer&Winter
fitted.Gamma <- predict(mod.scale.st.season, type="response", newdata = data.frame("lat"=rep(sites_DE_sel$lat,2),
                                                                                   "lon"=rep(sites_DE_sel$lon,2),
                                                                                   "elev"=rep(sites_DE_sel$elevation,2),
                                                                                   "season"=rep(c("summer","winter"), 
                                                                                                each=nrow(sites_DE_sel)),
                                                                                   "year"=rep(2000, 2*nrow(sites_DE_sel))))
fitted.shape         <- unique(fitted.Gamma[,2])
fitted.scale.summer  <- fitted.Gamma[1:nrow(sites_DE_sel),1]
fitted.scale.winter  <- fitted.Gamma[(nrow(sites_DE_sel)+1):(2*nrow(sites_DE_sel)),1]
fitted.median.summer <- qgamma(0.5, shape=fitted.shape, scale=fitted.scale.summer)
fitted.median.winter <- qgamma(0.5, shape=fitted.shape, scale=fitted.scale.winter)

df <- data.frame("lon"=sites_DE_sel$lon,
                 "lat"=sites_DE_sel$lat,
                 "scale.summer"=fitted.scale.summer,
                 "median.summer"=fitted.median.summer,
                 "scale.winter"=fitted.scale.winter,
                 "median.winter"=fitted.median.winter)

library(ggplot2)   
####### Summer 2000
gg_summer <- ggplot() + ggtitle("Summer 2000") +
  geom_point(aes(lon,lat,colour=median.summer), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated median", limits=c(min(c(df$median.summer,df$median.winter)),
                                                                                 max(c(df$median.summer,df$median.winter))))+
  theme_bw()

####### Winter 2000
gg_winter <- ggplot() + ggtitle("Winter 2000") +
  geom_point(aes(lon,lat,colour=median.winter), data = df, cex=.8)+
  scale_colour_gradientn(colours=rainbow(150), name="Estimated median", limits=c(min(c(df$median.summer,df$median.winter)),
                                                                                 max(c(df$median.summer,df$median.winter))))+
  theme_bw()

gg_summer <- gg_summer+ theme(legend.position = "none")
gg_winter <- gg_winter+ theme(legend.position = "bottom")

library(gridExtra)
png("spatial_seasons_gamma_median.png", width = 12, height = 5,units = 'in', res = 200)

grid.arrange(gg_summer,gg_winter, nrow=1,ncol=2)

dev.off()

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