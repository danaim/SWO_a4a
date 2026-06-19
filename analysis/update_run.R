# Update run with the same settings
# Data are not even close to what they used to be
library(FLCore)
library(FLa4a)
library(FLBRP)
library(ggplot2);theme_set(theme_bw())
rm(list = ls())

swo.stk <- readRDS("Robj/swo_stk.rds")
idx <- readRDS("Robj/swo_bio_idx.rds")

stk <- window(swo.stk, start = 1987)

stk <- setPlusGroup(stk, 5)
range(stk)[c('minfbar','maxfbar')] <- c(2,4) 

# 2020 settings
names(idx)
fmod2 <- ~s(year, k = 19) + s(age, k = 3)
srmod1 <- ~s(year, k = 17)
qmod <- list(~1,~s(year, k = 9), ~s(year, k = 3), ~1, ~s(year, k = 3))

fit <- sca(stk, idx, fmodel = fmod2, srmodel = srmod1, qmodel = qmod) 
a4a.stk <- stk + simulate(fit, nsim = 1000) # numbers are crazy

plot(a4a.stk)

fitdef <- sca(stk, idx)
a4a.fitd <- stk + simulate(fitdef, nsim = 1000)

plot(a4a.fitd)
