rm(list = ls())
library(reshape2)
library(FLXSA)
library(ggplot2)
library(ggplotFL)
library(FLCore)
library(mixdist)
library(plyr)
library(data.table)
library(dplyr)
library(FLa4a)
library(FLEDA)
library(FLAssess)
library(FLash)
library(FLBRP) 
library(kobe)
library(ggthemes);theme_set(theme_bw())
setwd("~/pCloudDrive/ICCAT/runs/a4a")
load("~/pCloudDrive/ICCAT/data/swom_mixed_data_bioIndexes_corrected_2008.rda")

### Stock object
stk <- swom

### Index object
idx <- u

### Stock object with discards
# load("~/pCloudDrive/ICCAT/data/swom_mixed_data_bioIndexes_dsc.rda")
stk_disc <- swom

### set the units
units(landings(stk)) <- units(landings(stk_disc)) <-"t"
units(catch(stk)) <- units(catch(stk_disc)) <- "t"
units(catch.n(stk)) <- units(catch.n(stk_disc)) <- ""
units(catch.wt(stk)) <- units(catch.wt(stk_disc)) <- "t"
units(stock.n(stk)) <- units(stock.n(stk_disc)) <- ""
units(stock.wt(stk)) <- units(stock.wt(stk_disc)) <- "t"
units(mat(stk)) <- units(mat(stk_disc)) <- ""
units(m(stk)) <- units(m(stk_disc)) <- "m"
units(harvest(stk)) <- units(harvest(stk_disc)) <- "f"

# catch.n(stk)[catch.n(stk_disc) == 0] <- 0.1

range(stk, c("minfbar", "maxfbar")) = c(2, 4)
range(stk_disc, c("minfbar", "maxfbar")) = c(2, 4)

### Biomass index
idx <- u
summary(idx)

idx_bio <- FLIndices()

for (i in 1:length(idx)) {
  dms <- list(age="all", year=range(idx[[i]])["minyear"]:range(idx[[i]])["maxyear"])
  idx_bio[[i]] <- FLIndexBiomass(FLQuant(NA, dimnames=dms))
  index(idx_bio[[i]]) <- quantSums(index(idx[[i]]))
}

names(idx_bio) <- names(idx)


for(i in 1:length(idx_bio)) {
  range(idx_bio[[i]],c("startf","endf")) = c(0,1)
  range(idx_bio[[i]],c("min","max")) = c(2,4)
  units(index(idx_bio[[i]])) <- ""
}


### Trim the stocks
stk <- trim(stk, year = 1985:2018)
stk_disc <- trim(stk_disc, year = 1985:2018)


### Lorenzen M
stk_L <- stk
par=FLPar(Linf=238.59,K=0.185,t0=-1.404,a=9.62e-3,b=3.06,l50=132,
          m1=20/c(mean(stock.wt(stk_L)["3"])),m2=FLPar(-0.309))

m=FLife:::lorenzen(stock.wt(stk_L),par[c("m1","m2")])
par["m1"]=0.2*par["m1"]/mean(m["3"])  

m(stk_L)=FLife:::lorenzen(stock.wt(stk_L),par[c("m1","m2")])

# fix units:
units(landings(stk_L)) <- "t"
units(catch(stk_L)) <- "t"
units(catch.n(stk_L)) <- ""
units(catch.wt(stk_L)) <- "t"
units(stock.n(stk_L)) <- ""
units(stock.wt(stk_L)) <- "t"
units(mat(stk_L)) <- ""
units(m(stk_L)) <- "m"
units(harvest(stk_L)) <- "f"


#############################################################################################
################################# Candidate runs ############################################
#############################################################################################
stk_L <- setPlusGroup(stk_L, 5)

fit9_L <- sca(stk_L, idx_bio, fmodel = fmod2, srmodel = srmod1, qmodel = qmod2)
a4a.stk9_L <- stk_L + fit9_L

res9_L <- residuals(fit9_L, stk_L, idx_bio)

fits9_L <- simulate(fit9_L, 1000)
a4a.stk9_L_sim <- stk_L + fits9_L
plot(a4a.stk9_L_sim)


######################################################################################################
################################# MCMC fit ###########################################################
######################################################################################################
fits9_L <- simulate(fit9_L, 1000)

mc <- SCAMCMC(mcmc=10000, mcsave=100)
fit9_L_mc <- sca(stk_L, idx_bio, fmodel = fmod2, srmodel = srmod1, qmodel = qmod2, fit="MCMC", mcmc=mc)
fitSumm(fit9_L_mc)

fit9_L_mc.mc <- as.mcmc(fit9_L_mc)
plot(fit9_L_mc.mc)

mc2 <- SCAMCMC(mcmc=10000, mcsave=100, mcprobe = 0.4)
fit9_L_mc2 <- sca(stk_L, idx_bio, fmodel = fmod2, srmodel = srmod1, qmodel = qmod2, fit="MCMC", mcmc=mc2)
fitSumm(fit9_L_mc2)

fit9_L_mc.mc2 <- as.mcmc(fit9_L_mc2)
plot(fit9_L_mc.mc2)


plot(FLStocks(mcmc = stk_L+fit9_L_mc.mc, mcmc2 = stk_L + fit9_L_mc.mc2, lorenM = a4a.stk9_L_sim))
