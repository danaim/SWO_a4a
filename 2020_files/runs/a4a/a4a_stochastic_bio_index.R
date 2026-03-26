### a4a run with stochastic age slicing and biomass indices ###
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
load("~/pCloudDrive/ICCAT/data/swom_mixed_data_bioIndexes.rda")

###################################################################################################
################################ Prepare stock and Index ##########################################
###################################################################################################
### Stock object
stk <- swom
summary(stk)

### set the units
units(landings(stk)) <- "t"
units(catch(stk)) <- "t"
units(catch.n(stk)) <- ""
units(catch.wt(stk)) <- "t"
units(stock.n(stk)) <- ""
units(stock.wt(stk)) <- "t"
units(mat(stk)) <- ""
units(m(stk)) <- "m"
units(harvest(stk)) <- "f"

catch.n(stk)[catch.n(stk) == 0] <- 0.1

range(stk, c("minfbar", "maxfbar")) = c(1, 3)

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
  range(idx_bio[[i]],c("min","max")) = c(1,3)
  units(index(idx_bio[[i]])) <- "t"
}
###################################################################################################
################################################# Explore Plots ###################################
###################################################################################################
### catch @ age
df<-as.data.frame(catch.n(stk))
df<-na.omit(df)
ggplot(df,aes(x = year, y = data))+ 
  geom_line()+ facet_grid(age ~.)

ggplot(as.data.frame(catch.n(stk)))+
  geom_point(aes(year,age,size=data))+
  scale_size_area(max_size=8,  guide="none") +
  xlab("Year")+ylab("Age")+
  scale_y_continuous(breaks=seq(0,8,2))+
  theme_bw()+theme(legend.position="none")

df<-as.data.frame(catch.n(stk))
df<-na.omit(df)
ggplot(df,aes(x = year, y = data,color= as.factor(age)))+ 
  geom_line(size = 1)+ggtitle("Catch at age MUT GSA 1")

## index plots
plot_idx <- as.data.frame(idx_bio)

ggplot(aes(year, data), data = plot_idx[plot_idx$slot == 'index',c(1,3,8,9)]) +
  geom_line(stat = 'identity')+ facet_grid(cname ~.)

###################################################################################################
################################################# Assessment ######################################
###################################################################################################
## Full time series and index assessment

## Explore the f model
fit <- sca(stk, idx_bio) #NOT strange fitting to catch
a4a.stk <- stk + fit
plot(a4a.stk)

## put fmodel
## simple separable factor
fmod <- ~factor(age) + factor(year)
fit2 <- sca(stk, idx_bio, fmodel = fmod) # NOT crazy recruitment
a4a.stk2 <- stk + fit2
plot(a4a.stk2)

# change fmodel to splines and play with some ks
# play also with age 

fmod <- ~s(year, k = 23) + s(age, k = 3)
fit3 <- sca(stk, idx_bio, fmodel = fmod)
a4a.stk3 <- stk + fit3
plot(a4a.stk3)

fmod <- ~s(year, k = 20) + s(age, k = 4)
fit4 <- sca(stk, idx_bio, fmodel = fmod)
a4a.stk4 <- stk + fit4
plot(a4a.stk4)

fmod <- ~s(year, k = 15) + s(age, k = 4)
fit5 <- sca(stk, idx_bio, fmodel = fmod)
a4a.stk5 <- stk + fit5

# add sr model to the very low fishing mortality
fmod <- ~s(year, k = 20) + s(age, k = 4)
srmod <- ~s(year,k = 20)
fit6 <- sca(stk, idx_bio, fmodel = fmod, srmodel = srmod)
a4a.stk6 <- stk + fit6
plot(a4a.stk6)


fmod <-  ~ s(replace(age, age>4, 4), k = ) + s(year, k=20)
srmod <- ~s(year,k = 20)
fit7 <- sca(stk, idx_bio, fmodel = fmod)
a4a.stk7 <- stk + fit7
plot(a4a.stk7)

stks <- FLStocks(catch_obs = stk, fit3 = a4a.stk3, fit4 = a4a.stk6, fit5 = a4a.stk7)
plot(stks)

stks[[1]]

# Fishing mortality 3D
wireframe(data ~ age + year, data = as.data.frame(harvest(a4a.stk10)), drape = TRUE,
          main="Fishing mortality", screen = list(x = -90, y=-45),
          col.regions = colorRampPalette(c("blue", "pink",'yellow'))(100))

###
stk2 <- trim(stk, year = 1985:2018) 

fmod <-  ~ s(replace(age, age>4, 4), k = 3) + s(year, k=20)
fit8 <- sca(stk2, idx_bio, fmodel = fmod)
a4a.stk8 <- stk2 + fit8
plot(a4a.stk8)

# Check smoother at age sensitivity
fmod <- ~s(year, k = 23) + s(age, k = 4)
fit9 <- sca(stk, idx_bio, fmodel = fmod)
a4a.stk9 <- stk + fit9

fmod <- ~s(year, k = 23) + s(age, k = 3)
fit10 <- sca(stk, idx_bio, fmodel = fmod)
a4a.stk10 <- stk + fit10

stks <- FLStocks(catch_obs = stk, fit10 = a4a.stk10, fit9 = a4a.stk9)
plot(stks)


# ad srmodel
fmod <-  ~ s(replace(age, age>4, 4), k = 4) + s(year, k=23)
fit3 <- sca(stk, idx_bio, fmodel = fmod)
a4a.stk3 <- stk + fit3
plot(a4a.stk3)


fmod <-  ~ s(replace(age, age>4, 4), k = 4) + s(year, k=23)
srmod <- ~s(year, k = 20)
fit4 <- sca(stk, idx_bio, fmodel = fmod, srmodel = srmod)
a4a.stk4 <- stk + fit4
plot(a4a.stk4)


stks <- FLStocks(catch_obs = stk, fit3 = a4a.stk3, fit4 = a4a.stk4)
plot(stks)

# fmod <- ~s(age, k = 3, by = breakpts(year, 1990))
# # add n1model
# 
# n1mod <- ~s(age, k =4)

wireframe(data ~ age + year, data = as.data.frame(harvest(a4a.stk)), drape = TRUE,
          main="Fishing mortality", screen = list(x = -90, y=-45),
          col.regions = colorRampPalette(c("blue", "pink",'yellow'))(100))


## Built some extra diagnostics
qqmath(res)
res <- residuals(fit1, stk, idx_bio)
plot(res)

plot(fit, stk)

res


fit.pred <- predict(fit)
lapply(fit.pred, names)

###################################################################################################
############################################## Retrospective ######################################
###################################################################################################

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
    fmod <- substitute(~ s(replace(age, age>4, 4), k = 4) + s(year, k=KY), list(KY=KY))
    args$fmodel <- as.formula(fmod)
    fit <- do.call("sca", args)
    args$stock + fit
  })
  FLStocks(lst0)
}

retro.lst <- retro(stk, idx_bio, k=23, srmodel= srmod)
plot(retro.lst, main="Retrospective analysis")

fmod <-  ~ s(replace(age, age>4, 4), k = 4) + s(year, k=23)
###################################################################################################
############################################## Simulations #######################################
###################################################################################################

fits <- simulate(fit4, 1000)
flqs <- FLQuants(sim=iterMedians(stock.n(fits)), det=stock.n(fit4))

sims <- as.data.frame(flqs)
ggplot(aes(year,data), data = sims[sims$qname == 'sim',c(1,2,7,8)]) +
  geom_line() + facet_wrap(age~.)

sim_stk <- stk + fits
plot(sim_stk)
