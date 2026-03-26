### a4a run with stochastic age slicing and biomass indices ###
### trimmed time series #######################################
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

range(stk, c("minfbar", "maxfbar")) = c(2, 4)
# stk <- setPlusGroup(stk,5)
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
  units(index(idx_bio[[i]])) <- "t"
}

### Trim the time series
stk <- trim(stk, year = 1985:2018)

###################################################################################################
########################################### Explore Plots #########################################
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
  geom_line(size = 1)+ggtitle("Catch at age")

## index plots
plot_idx <- as.data.frame(idx_bio)

ggplot(aes(year, data), data = plot_idx[plot_idx$slot == 'index',c(1,3,8,9)]) +
  geom_line(stat = 'identity')+ facet_grid(cname ~.)

temp <- FLIndices()
dms <- list(age=range(stk)["min"]:range(stk)["max"], 
            year=range(stk)["minyear"]:range(stk)["maxyear"])
temp[[1]] <- FLIndex(FLQuant(NA, dimnames=dms))
index(temp[[1]]) <- catch.n(stk)
plot(temp[[1]], type="internal")
###################################################################################################
####################################### Assessment ################################################
###################################################################################################


###################################################################################################
###################################################################################################
# change fmodel to splines and play with some ks
# play also with age 

# Comments on separable splines model:
## it is very sensitive to k on age (unresonable results for k>3)
## stick with k = 3 = floor(nage/2)

fmod <- ~s(year, k = 17) + s(age, k = 3)
fit3 <- sca(stk, idx_bio, fmodel = fmod)
a4a.stk3 <- stk + fit3
plot(a4a.stk3)

res <- residuals(fit3, stk, idx_bio)
res_p <- residuals(fit3, stk, idx_bio, type = "pearson")
res_sp <- residuals(fit3, stk, idx_bio, type = "standard.pearson")
plot(res)

# Fishing mortality 3D
wireframe(data ~ age + year, data = as.data.frame(harvest(a4a.stk3)), drape = TRUE,
          main="Fishing mortality", screen = list(x = -90, y=-45),
          col.regions = colorRampPalette(c("blue", "pink",'yellow'))(100))


####################################################################################
####################################################################################
#### Best runs bad residuals plus group 6

fmod <- ~s(year, k = 15) + s(age, k = 3) + te(age, year,k=c(4,7)) ####BEST residuals

# add sr model
fmod <- ~s(year, k = 15) + s(age, k = 3) + te(age, year,k=c(4,10))
srmod <- ~s(year, k = 15)
fit7 <- sca(stk, idx_bio, fmodel = fmod, srmodel = srmod)
a4a.stk7 <- stk + fit7
plot(a4a.stk7)
res <- residuals(fit7, stk, idx_bio)
plot(res)
AIC(fit7)
rec(a4a.stk7)

# Fishing mortality 3D
wireframe(data ~ age + year, data = as.data.frame(harvest(a4a.stk7)), drape = TRUE,
          main="Fishing mortality", screen = list(x = -90, y=-45),
          col.regions = colorRampPalette(c("blue", "pink",'yellow'))(100))




fmod <- ~s(year, k = 17) + s(age, k = 3)
srmod <- ~s(year, k = 15)
qmod <- list(~s(year, k = 15),~1,~s(year, k =5),~s(year, k =5),~1)
fit10 <- sca(stk, idx_bio, fmodel = fmod, qmodel = qmod, srmodel = srmod)
a4a.stk10 <- stk + fit10
plot(a4a.stk10)
res <- residuals(fit10, stk, idx_bio)
plot(res)
AIC(fit10)

fmod <- ~s(year, k = 17) + s(replace(age, age>4, 4), k = 4)
srmod <- ~s(year, k = 15)
fit8 <- sca(stk, idx_bio, fmodel = fmod, srmodel = srmod)
a4a.stk8 <- stk + fit8
plot(a4a.stk8)
res <- residuals(fit8, stk, idx_bio)
plot(res)
AIC(fit8)

stks <- FLStocks(catch_obs = stk, fit6 = a4a.stk6, fit7 = a4a.stk7, fit8 = a4a.stk8)
plot(stks)

# Increase degrees of freedom of the smoother splines to get better fit to age 6
# as a result we get very low ssb

fmod <- ~s(year, k = 17) + s(age, k = 4)
srmod <- ~s(year, k = 15) 
fit9 <- sca(stk, idx_bio, fmodel = fmod, srmodel = srmod)
a4a.stk9 <- stk + fit9
plot(a4a.stk9)
res <- residuals(fit9, stk, idx_bio)
plot(res)

fmod2<-  ~s(replace(age, age > 2, 2), k = 3) + s(year, k = 6) + te(age, year,k=c(3,7))

# # Take out spanish lls if it makes any difference (it doesn't !!)
# idx_bio2 <- FLIndices(idx_bio[[1]],idx_bio[[2]],idx_bio[[3]],idx_bio[[4]])
# 
# fmod <- ~s(year, k = 17) + s(age, k = 4)
# srmod <- ~s(year, k = 15)
# fit10 <- sca(stk, idx_bio2, fmodel = fmod, srmodel = srmod)
# a4a.stk10 <- stk + fit10
# plot(a4a.stk10)
# res <- residuals(fit10, stk, idx_bio2)
# plot(res)
# 
# stks <- FLStocks(catch_obs = stk, fit9 = a4a.stk9, fit10 = a4a.stk10)
# plot(stks)

## Set plus group 5:
stk2 <- setPlusGroup(stk, 5)

# add sr model
fmod <- ~s(year, k = 17) + s(age, k = 3)
srmod <- ~s(year, k = 15)
fit11 <- sca(stk2, idx_bio, fmodel = fmod, srmodel = srmod)
a4a.stk11 <- stk2 + fit11
plot(a4a.stk11)
res <- residuals(fit11, stk2, idx_bio)
plot(res)
AIC(fit11)

fmod <- ~s(year, k = 17) + s(replace(age, age>4, 4), k = 4)
srmod <- ~s(year, k = 15)
fit12 <- sca(stk2, idx_bio, fmodel = fmod, srmodel = srmod)
a4a.stk12 <- stk2 + fit12
plot(a4a.stk12)
res <- residuals(fit12, stk2, idx_bio)
plot(res)
bubbles(res)
AIC(fit12)

stks <- FLStocks(catch_obs = stk2, fit11 = a4a.stk11, fit12 = a4a.stk12)
plot(stks)

brp4 <- FLBRP(a4a.stk12)
summary(brp4)
harvest(brp4)

#Reference points
model = "geomean"
spe_srr <- fmle(as.FLSR(a4a.stk10, model = model))
spe_brp <- brp(FLBRP(a4a.stk10, sr = spe_srr))
ref_points<- refpts(spe_brp)


retro <- function(stk, idxs, retro=5, kfrac="missing", k, ...){
  args <- list(...)
  if(missing(kfrac)) kfrac <- 0.5
  lst0 <- split(0:retro, 0:retro)
  lst0 <- lapply(lst0, function(x){
    yr <- range(stk)["maxyear"] - x
    args$stock <- window(stk, end=yr)
    args$indices <- FLIndices(window(idxs, end=yr))
    KY <- unname(k - floor(x*kfrac))
    # NOTE THIS HAS TO BE ADAPTED FOR THE MODEL USED  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # fmod <- substitute(~ s(replace(age, age>4, 4), k = 4) + s(year, k=KY), list(KY=KY))
    fmod <- substitute(~ s(age, k = 4) + s(year, k=KY), list(KY=KY))
    args$fmodel <- as.formula(fmod)
    fit <- do.call("sca", args)
    args$stock + fit
  })
  FLStocks(lst0)
}
retro.lst <- retro(stk, idx_bio, k=17, srmodel= srmod)
plot(retro.lst, main="Retrospective analysis")
