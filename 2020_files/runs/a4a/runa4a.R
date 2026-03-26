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
load("~/pCloudDrive/ICCAT/data/swom_ageit_data_bioIndexes.rda")
load("~/pCloudDrive/ICCAT/data/swom_continuityStrict_data.rda")
# Prepare stk
stk <- swom
summary(stk)

## set the units
units(landings(stk)) <- "t"
units(catch(stk)) <- "t"
units(catch.n(stk)) <- "numbers"
units(catch.wt(stk)) <- "t"
units(stock.n(stk)) <- "numbers"
units(stock.wt(stk)) <- "t"
units(mat(stk)) <- ""
units(m(stk)) <- "m"
units(harvest(stk)) <- "f"

catch.n(stk)[catch.n(stk) == 0] <- 0.1

range(stk, c("minfbar", "maxfbar")) = c(1, 3)

stk <- trim(stk, year = 1987:2018)

# Biomass index
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

for(i in 1:length(idx)) {
  range(idx[[i]],c("startf","endf")) = c(0,1)
  range(idx[[i]],c("min","max")) = c(1,3)
  units(index(idx[[i]])) <- ""
}

## Plots

df<-as.data.frame(catch.n(stk))
df<-na.omit(df)
ggplot(df,aes(x = year, y = data))+ 
  geom_line()+ facet_grid(age ~.)

plot_idx <- as.data.frame(idx_bio)

ggplot(aes(year, data), data = plot_idx[plot_idx$slot == 'index',c(1,3,8,9)]) +
  geom_line(stat = 'identity')+ facet_grid(cname ~.)

plot_idx2 <- as.data.frame(idx_bio)

ggplot(aes(year, data), data = plot_idx[plot_idx$slot == 'index',c(1,3,8,9)]) +
  geom_line(stat = 'identity')+ facet_grid(cname ~.)
## Assessment

# fmodel
fmod <-  ~s(year, k = 20) + s(age, k = 5) 
srmod <- ~ bevholt(CV=0.2)

fit <- sca(stock = stk, indices = idx_bio, fmodel = fmod)
a4a.stk <- stk + fit
plot(stk + fit)
stks <- FLStocks(catch_obs = stk, assess = a4a.stk)
plot(stks)

# Fishing mortality 3D
wireframe(data ~ age + year, data = as.data.frame(harvest(a4a.stk)), drape = TRUE,
          main="Fishing mortality", screen = list(x = -90, y=-45),
          col.regions = colorRampPalette(c("blue", "pink",'yellow'))(100))

fit2 <- sca(stock = stk, indices = idx_bio, fmodel = fmod, srmodel = srmod)
a4a.stk2 <- stk + fit2
plot(stk + fit2)
stks <- FLStocks(catch_obs = stk, assess = a4a.stk, bev_ass = a4a.stk2)
plot(stks)



## residuals
# diagnostics
res <- residuals(fit, stk, idx_bio)
plot(res)

plot(fit, stk)
plot(fit, idx_bio) 
