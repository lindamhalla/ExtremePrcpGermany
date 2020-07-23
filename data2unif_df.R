library(evd)

# install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")
library(TruncatedDistributions)
load("data/german-precipitation.RData")
load("data/sites_DE_sel_elev.Rdata")
load("data/fitted_quant_90.Rdata")
load("data/fitted.gpd.full_ss.Rdata")
load("data/fitted.pr.exc.Rdata")
load("data/fitted.gamma.param.Rdata")
load("data/fitted.beta.param.Rdata")
load("data/fitted.pr.obs.Rdata")
load("data/m_gpd_full_ss_sub_exc.Rdata")

### remove Zugspitze
data_sel     <- data_sel[,-13]
sites_DE_sel <- sites_DE_sel[-13,]

############### Transform data to uniform using the beta distribution for observations between zero and the threshold
####### 
anom.training.unif <- data_sel

### Transform data to the uniform scale
fct2unif <- function(y,u,sh1.beta,sh2.beta,sh.gpd,sc.gpd,pu,p0){
  pu <- pu/p0 #we need the probability of exceeding the thd conditional on being positive
  if(is.na(y)) ret <- NA
  else{
    if(y==0)
      ret <- 1-p0
    else if((y>0)&(y <= u)) 
      ret <- (1-p0)+p0*(1-pu)*pbeta(y/u, shape1=sh1.beta, shape2=sh1.beta)
    if(y>u) 
      ret <- (1-p0)+p0*(1-pu*(1+sh.gpd*(y-u)/sc.gpd)^(-1/sh.gpd))
  }
  return(ret)
}

for(i in 1:nrow(data_sel)){
  if((i%%500)==0) print(paste0("i=",i))
  for(j in 1:ncol(data_sel)){
    if((j%%500)==0) print(paste0("j=",j))
    anom.training.unif[i,j] <- fct2unif(data_sel[i,j],fitted.quant.mat[i,j],
                                        sh1.beta=fitted.beta.param$shape1[i,j],
                                        sh2.beta=fitted.beta.param$shape2[i,j],
                                        sh.gpd=fitted.gpd.full$shape[i,j], 
                                        sc.gpd=fitted.gpd.full$scale[i,j],
                                        pu=fitted.pr.exc$fitted.pr.exc.mat[i,j],
                                        p0=fitted.pr.obs$fitted.pr.obs.mat[i,j])
  }
}
save(anom.training.unif, file="data/anom.training.unif.beta.Rdata")

############### Transform data to uniform using the truncated Gamma distribution for observations between zero and the threshold
####### 
anom.training.unif <- data_sel

### Transform data to the uniform scale
fct2unif <- function(y,u,sh.gamma,sc.gamma,sh.gpd,sc.gpd,pu,p0){
  pu <- pu/p0 #we need the probability of exceeding the thd conditional on being positive
  if(is.na(y)) ret <- NA
  else{
    if(y==0)
      ret <- 1-p0
    else if((y>0)&(y <= u)) 
      ret <- (1-p0)+p0*(1-pu)*ptgamma(y, shape=sh.gamma, scale=sc.gamma, a=0, b = u)
    if(y>u) 
      ret <- (1-p0)+p0*(1-pu*(1+sh.gpd*(y-u)/sc.gpd)^(-1/sh.gpd))
  }
  return(ret)
}

for(i in 1:nrow(data_sel)){
  if((i%%500)==0) print(paste0("i=",i))
  for(j in 1:ncol(data_sel)){
    if((j%%500)==0) print(paste0("j=",j))
    anom.training.unif[i,j] <- fct2unif(data_sel[i,j],fitted.quant.mat[i,j],
                                        sh.gamma=fitted.gamma.param$shape[i,j], sc.gamma=fitted.gamma.param$scale[i,j],
                                        sh.gpd=fitted.gpd.full$shape[i,j], sc.gpd=fitted.gpd.full$scale[i,j],
                                        pu=fitted.pr.exc$fitted.pr.exc.mat[i,j],
                                        p0=fitted.pr.obs$fitted.pr.obs.mat[i,j])
  }
}
save(anom.training.unif, file="data/anom.training.unif.Rdata")

#####################################################################################################################################################################
#####################################################################################################################################################################
############### Save data in a list

DE_precip_data_fit <- list("obs"=data_sel,
                          "sites"=sites_DE_sel,
                          "date"=dates_sel,
                          "gpd.sc"=fitted.gpd.full$scale,
                          "gpd.sh"=fitted.gpd.full$shape,
                          "quant90"=fitted.quant.mat,
                          "pr.exc"=fitted.pr.exc$fitted.pr.exc.mat,
                          "gamma.sc"=fitted.gamma.param$scale,
                          "gamma.sh"=fitted.gamma.param$shape,
                          "beta.sh1"=fitted.beta.param$shape1,
                          "beta.sh2"=fitted.beta.param$shape2,
                          "p0"=fitted.pr.obs$fitted.pr.obs.mat,
                          "data.unif.p0.gamma"=anom.training.unif,
                          "data.unif.p0.beta"=anom.training.unif.beta)

save(DE_precip_data_fit, file="data/DE_precip_data_fit_unif.Rdata")

