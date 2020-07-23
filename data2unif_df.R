load("data/anom.training.unif.beta.Rdata")
anom.training.unif.beta <- anom.training.unif
load("data/anom.training.unif.Rdata")
load("data/fitted.gpd.full_ss.Rdata")
load("data/fitted_quant_90.Rdata")
load("data/fitted.pr.exc.Rdata")
load("data/fitted.gamma.param.Rdata")
load("data/fitted.pr.obs.Rdata")
load("data/fitted.beta.param.Rdata")

load("data/german-precipitation.RData")
load("data/sites_DE_sel_elev.Rdata")
### remove Zugspitze
data_sel     <- data_sel[,-13]
sites_DE_sel <- sites_DE_sel[-13,]

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

