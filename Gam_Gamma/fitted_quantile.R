library(sp)
library(stringr)
library(mgcv)

load("data/german-precipitation.RData")

load("data/sites_DE_sel_elev.Rdata")

### Spatio-temporal quantile (threshold) (years+months+space)
load("data/mod.scale.st.year.months.Rdata")

### remove Zugspitze
data_sel     <- data_sel[,-13]
sites_DE_sel <- sites_DE_sel[-13,]

####### 
thd.level        <- 0.90
fitted.quant.mat <- NULL

for(s in 1:nrow(sites_DE_sel)){
  fitted.Gamma <- predict(mod.scale.st, type="response", newdata = data.frame("lat"=rep(sites_DE_sel$lat[s],nrow(dates_sel)),
                                                                              "lon"=rep(sites_DE_sel$lon[s],nrow(dates_sel)),
                                                                              "elev"=rep(sites_DE_sel$elevation[s],nrow(dates_sel)),
                                                                              "month"=dates_sel$month,
                                                                              "year"=dates_sel$year))
  fitted.shape <- exp(unique(fitted.Gamma[,2]))
  fitted.mean  <- fitted.Gamma[,1]
  fitted.quant <- qgamma(thd.level, shape=1/fitted.shape, scale=fitted.mean*fitted.shape)
  
  fitted.quant.mat <- cbind(fitted.quant.mat, fitted.quant)
}

save(fitted.quant.mat, file="data/fitted_quant_90.Rdata")
