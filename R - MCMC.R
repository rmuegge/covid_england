library(gsubfn)
library(ggplot2)
library(readxl)
library(reshape2)
library(rgdal)
library(leaflet)
library(spdep)
library(dplyr)
library(CARBayesST)
library(coda)
library(geosphere)

covid <- read.csv("coviddata.csv")
Bdry <- readOGR(dsn="Local_Authority_Districts__May_2020__Boundaries_UK_BFC.shp")
week1 <- covid[which(covid$week=="21_01"),]
Covid.Bdry <- merge(x=Bdry,y=week1,by.x="LAD20CD",by.y="area",all.x=FALSE)
Covid.Bdry.ll <- spTransform(Covid.Bdry,CRS("+proj=longlat +datum=WGS84 +no_defs"))
lookup <- data.frame(Code=Covid.Bdry@data$LAD20CD,spatialorder=1:nrow(Covid.Bdry@data))
covid.temp <- merge(x=covid,y=lookup,by.x="area",by.y="Code")
covid.ordered <- arrange(covid.temp,week,spatialorder)
W.nb <- poly2nb(Covid.Bdry.ll,row.names=Covid.Bdry.ll@data$LAD20CD)
W <- nb2mat(W.nb,style="B")

fit1 <- ST.CARar(deaths~offset(log(expecteddeaths)),family="poisson",
  data=covid.ordered,W=W,AR=1,MALA=TRUE,burnin=200000,n.sample=2200000,thin=1000,verbose=TRUE)

#saveRDS(fit1,"ST_ar1.rds")

fit2 <- ST.CARar(deaths~offset(log(expecteddeaths)),family="poisson",
                 data=covid.ordered,W=W,AR=2,MALA=TRUE,burnin=200000,n.sample=2200000,thin=1000,verbose=TRUE)

#saveRDS(fit2,"ST_ar2.rds")

# Convergence check: traceplots
sims <- fit2$samples
mc_beta <- mcmc(sims$beta)
mc_beta <- ggs(mc_beta)
mc_beta$Parameter <- c("beta0")
beta_plot <- ggs_traceplot(mc_beta) + ylab("Value") + geom_line(colour="#003865") + xlab("")

mc_rho3 <- mcmc(sims$rho)
mc_rho3 <- ggs(mc_rho3)
mc_alpha1 <- mc_rho3[which(mc_rho3$Parameter=="rho1.T"),]
mc_alpha2 <- mc_rho3[which(mc_rho3$Parameter=="rho2.T"),]
mc_rho <- mc_rho3[which(mc_rho3$Parameter=="rho.S"),]

mc_rho$Parameter <- c("rho")
rho_plot <- ggs_traceplot(mc_rho) + ylab("Value") + geom_line(colour="#003865") + xlab("")

mc_alpha1$Parameter <- c("alpha1")
alpha1_plot <- ggs_traceplot(mc_alpha1) + ylab("Value") + geom_line(colour="#003865") + xlab("")

mc_alpha2$Parameter <- c("alpha2")
alpha2_plot <- ggs_traceplot(mc_alpha2) + ylab("Value") + geom_line(colour="#003865") + xlab("")

mc_tau2 <- mcmc(sims$tau2)
mc_tau2 <- ggs(mc_tau2)
mc_tau2$Parameter <- c("tau^2")
tau2_plot <- ggs_traceplot(mc_tau2) + ylab("Value") + geom_line(colour="#003865")

mc_phi <- mcmc(sims$phi)
mc_phi <- ggs(mc_phi)
mc_phi <- mc_phi[which(mc_phi$Parameter==1),]
mc_phi$Parameter <- c("phi_11")
phi_plot <- ggs_traceplot(mc_phi) + ylab("Value") + geom_line(colour="#003865")

grid.arrange(beta_plot,rho_plot,alpha1_plot,alpha2_plot,tau2_plot,phi_plot,nrow=3)

#DIC is lower for ar2. Work with ar2 results.
res_ar1$modelfit
res_ar2$modelfit

