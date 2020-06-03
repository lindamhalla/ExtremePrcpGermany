library(sp)
library(stringr)

DATA="~/research/data/ECAD/ECA_blend_rr/"
setwd("~/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/")
OUT="~/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/data/"
load(paste0(OUT, "german-precipitation.RData"))

load("/Users/worklinda/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/sites_DE_sel_elev.Rdata")

### Spatio-temporal quantile (threshold) (years+months+space)
load("/Users/worklinda/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Thomas/mod.scale.st.year.months.Rdata")
####### 
thd.level        <- 0.75
fitted.quant.mat <- NULL
for(yy in 1:length(unique(dates_sel$year))){
  fitted.Gamma <- predict(mod.scale.st, type="response", newdata = data.frame("lat"=rep(sites_DE_sel$lat,12),
                                                                              "lon"=rep(sites_DE_sel$lon,12),
                                                                              "elev"=rep(sites_DE_sel$elevation,12),
                                                                              "month"=rep(1:12, 
                                                                                          each=nrow(sites_DE_sel)),
                                                                              "year"=rep(unique(dates_sel$year)[yy], 
                                                                                         12*nrow(sites_DE_sel))))
  fitted.shape <- exp(unique(fitted.Gamma[,2]))
  fitted.mean  <- fitted.Gamma[,1]
  fitted.quant <- qgamma(thd.level, shape=1/fitted.shape, scale=fitted.mean*fitted.shape)
  
  fitted.quant.mat <- rbind(fitted.quant.mat,matrix(fitted.quant, byrow=TRUE, ncol=length(sites_DE_sel$lon)))
}

dates_quant_mat <- cbind(rep(unique(dates_sel$year), each=12), rep(1:12, length(unique(dates_sel$year))))
quantile.mat    <- list("dates"=dates_quant_mat,"quant"=fitted.quant.mat)
save(quantile.mat, file="quantile.mat.Rdata")

#merge quantiles and dates

fct.merge <- function(ss){
  dates_fitted_quant <- data.frame("year"=quantile.mat$dates[,1],
                                   "month"=quantile.mat$dates[,2],
                                   "quant"=quantile.mat$quant[,ss])
  merge.try <- merge.data.frame(dates_sel, dates_fitted_quant, all.x=TRUE, 
                                by = intersect(names(dates_sel), names(dates_fitted_quant)))
  return(merge.try$quant)
}

fitted.quant.mat.list <- lapply(1:ncol(data_sel), fct.merge) 
fitted.quant.mat      <- do.call(cbind,fitted.quant.mat.list)

save(fitted.quant.mat, file="/Users/worklinda/Desktop/HEC_Montreal/Postdoc/Daniela_Thomas_precipitation/Gam_Gamma/fitted_quant_75.Rdata")
