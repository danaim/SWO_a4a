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
load("~/pCloudDrive/ICCAT/data/swom_mixed_data_bioIndexes_dsc.rda")
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
########################### Various exploratory plots #######################################
#############################################################################################
plot(catch.wt(stk))

ggplot(aes(year,data, color = as.factor(age)), data = as.data.frame(catch.wt(stk_disc))) +
  geom_line(stat = 'identity')

#### discards included
catch_comp <- as.data.frame(catch.n(stk))[,c(1,2,7)]
catch_comp$stk <- 'no_disc'

temp <- as.data.frame(catch.n(stk_disc))[,c(1,2,7)]
temp$stk <- 'disc'
catch_comp <- rbind(catch_comp, temp); rm(temp)

ggplot(aes(year,data, color = as.factor(stk)), data = catch_comp)+theme(legend.position="none")+
  geom_line() + facet_grid(age ~. )
#############################################################################################
################################# Candidate runs ############################################
#############################################################################################
stk <- setPlusGroup(stk, 5)
stk_L <- setPlusGroup(stk_L, 5)
stk_disc <- setPlusGroup(stk_disc, 5)

# fmod1 <- ~s(year, k = 17) + s(age, k = 3) + s(year, k = 17, by = as.numeric(age == 0))

fmod2 <- ~s(year, k = 17) + s(age, k = 3)
fmod3 <- ~s(year, k = 17) + s(replace(age, age>4, 4), k = 3)
# fmod4 <- ~s(year, k = 17) + s(age, k = 4)
# fmod5 <- ~s(year, k = 17) + s(replace(age, age>4, 4), k = 4)


srmod1 <- ~s(year, k = 15) ### this is the final srmodel
srmod2 <- ~s(year, k = 18)

# there is a threashold in 18 knots after which the model 
# gives a very high recruitment
srmod3 <- ~s(year, k = 19)

#vmod <- list(~s(age, k = 4), ~1, ~1,~1,~1,~1)
qmod1 <- list(~1,~s(year, k = 5), ~s(year, k = 3), ~s(year, k = 3),~1)
qmod2 <- list(~1,~1, ~s(year, k = 3), ~s(year, k = 3),~1)

### results in crazy recruitment that also drives away the fit of the catch just before 2010
# fit1 <- sca(stk, idx_bio, fmodel = fmod1, srmodel = srmod1)
# a4a.stk1 <- stk + fit1

fit2 <- sca(stk, idx_bio, fmodel = fmod2, srmodel = srmod1)
a4a.stk2 <- stk + fit2

fit3 <- sca(stk, idx_bio, fmodel = fmod3, srmodel = srmod1)
a4a.stk3 <- stk + fit3

## lower SSB and unstable in retro
# fit4 <- sca(stk, idx_bio, fmodel = fmod4, srmodel = srmod1)
# a4a.stk4 <- stk + fit4

## The same as fit3
# fit5 <- sca(stk, idx_bio, fmodel = fmod5, srmodel = srmod1)
# a4a.stk5 <- stk + fit5

fit6 <- sca(stk, idx_bio, fmodel = fmod3, srmodel = srmod1, qmodel = qmod1)
a4a.stk6 <- stk + fit6

# fit 7 fit 8 rejected due to recruitment model
# fit7 <- sca(stk, idx_bio, fmodel = fmod2, srmodel = srmod3)
# a4a.stk7 <- stk + fit7
# 
# fit8 <- sca(stk, idx_bio, fmodel = fmod3, srmodel = srmod2)
# a4a.stk8 <- stk + fit8

### Two candidate ones:

fit9 <- sca(stk, idx_bio, fmodel = fmod2, srmodel = srmod1, qmodel = qmod2)
a4a.stk9 <- stk + fit9

# fit10 <- sca(stk, idx_bio, fmodel = fmod2, srmodel = srmod1, qmodel = qmod1)
# a4a.stk10 <- stk + fit10
# 
# plot(FLStocks(fit10 = a4a.stk10, fit2 = a4a.stk2, 
#               fit3 = a4a.stk3, fit9 = a4a.stk9, fit6 = a4a.stk6))

AIC(fit9);AIC(fit10)
BIC(fit9);BIC(fit10)
fitSumm(fit9);fitSumm(fit10)

## Check Lorenz M and discards
fit9_L <- sca(stk_L, idx_bio, fmodel = fmod2, srmodel = srmod1, qmodel = qmod2)
a4a.stk9_L <- stk_L + fit9_L

# fit10_L <- sca(stk_L, idx_bio, fmodel = fmod2, srmodel = srmod1, qmodel = qmod1)
# a4a.stk10_L <- stk_L + fit10_L

fit9_disc <- sca(stk_disc, idx_bio, fmodel = fmod2, srmodel = srmod1, qmodel = qmod1)
a4a.stk9_disc <- stk_disc + fit9_disc

# fit10_disc <- sca(stk_disc, idx_bio, fmodel = fmod2, srmodel = srmod1, qmodel = qmod1)
# a4a.stk10_disc <- stk_disc + fit10_disc

plot(FLStocks(fit9 = a4a.stk10, fit10_L = a4a.stk10_L, fit10_disc = a4a.stk10_disc))

plot(FLStocks(fit9 = a4a.stk9, fit9_L = a4a.stk9_L, fit9_disc = a4a.stk9_disc))

# plot(a4a.stk9)

# res1 <- residuals(fit1, stk, idx_bio)
# res2 <- residuals(fit2, stk, idx_bio)
# res3 <- residuals(fit3, stk, idx_bio)
# res6 <- residuals(fit6, stk, idx_bio)
res9 <- residuals(fit9, stk, idx_bio)
res10 <- residuals(fit10, stk, idx_bio)

res_p <- residuals(fit9, stk, idx_bio, type = "pearson")
plot(res_p, main = NULL)
# res_sp <- residuals(fit1, stk, idx_bio, type = "deviances")

plot(res9)


plot(res_p)
plot(res_sp)

ggplot(aes(year,data), data = as.data.frame(harvest(a4a.stk9))) +
  geom_point() + geom_line() + geom_smooth() +
  facet_wrap(~age, scales = 'free')

############################################################################
########################## retro ###########################################
############################################################################

retro <- function(stk, idxs, retro=3, kfrac="missing", k, ...){
  args <- list(...)
  if(missing(kfrac)) kfrac <- 0.5
  lst0 <- split(0:retro, 0:retro)
  lst0 <- lapply(lst0, function(x){
    yr <- range(stk)["maxyear"] - x
    args$stock <- window(stk, end=yr)
    args$indices <- FLIndices(window(idxs, end=yr))
    KY <- unname(k - floor(x*kfrac))
    # NOTE THIS HAS TO BE ADAPTED FOR THE MODEL USED  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fmod <- substitute(~s(year, k = KY) + s(age, k = 3), list(KY=KY))
    # fmod <- substitute(~s(year, k = KY) + s(replace(age, age>4, 4), k = 3), list(KY=KY))
    args$fmodel <- as.formula(fmod)
    fit <- do.call("sca", args)
    args$stock + fit
  })
  FLStocks(lst0)
}
retro.lst <- retro(stk_disc, idx_bio, k=17, srmodel= srmod1, qmodel = qmod2)
plot(retro.lst, main="Retrospective analysis")
############################################################################
############################################################################
############################################################################


wireframe(data ~ age + year, data = as.data.frame(harvest(a4a.stk10)), drape = TRUE,
          screen = list(x = -90, y=-45), zlab = 'F',
          col.regions = colorRampPalette(c("blue", "pink",'yellow'))(100))

landings.n(a4a.stk9)<-catch.n(a4a.stk9)
landings.wt(a4a.stk9_L)<-catch.wt(a4a.stk9)
landings(a4a.stk9)<-catch(a4a.stk9)
discards(a4a.stk9)<-0
discards.n(a4a.stk9)<-0
discards.wt(a4a.stk9)<-0

model = "bevholt"
srr <- fmle(as.FLSR(a4a.stk9, model = model),control=list(trace=0),fixed=list(b=100))
brp <- brp(FLBRP(a4a.stk9, sr = srr))
ref_points<- refpts(brp)

ref_points<- refpts(spe_brp)
plot(spe_brp)

landings.sel(spe_brp)

ref_points_df <- as.data.frame(ref_points[c('f0.1','msy','virgin'),
                                          c('harvest', 'yield','rec','ssb','biomass')])
ref_points_df <- ref_points_df[,c(1,2,4)]
ref_points_df <- reshape2::dcast(ref_points_df, refpt~quant)
class(ref_points_df)
plot(landings.sel(spe_brp))

sel <- as.data.frame(landings.sel(spe_brp))
sel$sel_1 <- 1 - exp(-sel$data)

ggplot(aes(age,data), data = sel) + geom_point() + geom_line()

(catch.n(a4a.stk9)/stock.n(a4a.stk9))



fits9_L <- simulate(fit9_L, 1000)
a4a.stk9_L_sim <- stk_L + fits9_L
plot(a4a.stk9_L_sim, probs=c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95))

######################################################################################################
################################# MCMC fit ###########################################################
######################################################################################################
fits9_L <- simulate(fit9_L, 10)
fitss <- simulate(fit9_L,10)

plot(FLStocks(a = a4a.stk9_L + fits9_L, b = a4a.stk9_L+fitss))
mc <- SCAMCMC(mcmc=10000, mcsave=100)
fit9_L_mc <- sca(stk_L, idx_bio, fmodel = fmod2, srmodel = srmod1, qmodel = qmod2, fit="MCMC", mcmc=mc)

fitSumm(fit9_L_mc)

fit9_L_mc.mc <- as.mcmc(fit9_L_mc)
plot(fit9_L_mc.mc)



plot(FLStocks(mcmc = stk_L+fit9_L_mc, lorenM = a4a.stk9_L_sim))

catch_diagn <- computeCatchDiagnostics(fit9_L,stk_L,idx_bio)
plot(catch_diagn)
########################################################################################
############################# kobe plot ################################################
########################################################################################

data(yft)
kobePhase(subset(yft,year==2014))+
  geom_point(aes(stock,harvest))

subset(yft,year==2014 & scenario == 1)
x<-as.numeric(fbar(a4a.stk9)/ref_points['f0.1','harvest'])
y<-as.numeric(ssb(a4a.stk9))/ref_points['f0.1','ssb']
df_kobe<-data.frame(stock = c(y,1.5),harvest = c(x,0.5))

kobePhase(df_kobe)+
  geom_point(aes(stock,harvest))

stock(a4a.stk9)
