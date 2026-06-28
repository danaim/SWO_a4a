# Update run with the same settings
# Data are not even close to what they used to be
library(FLCore)
library(FLa4a)
library(FLBRP)
library(ggplot2);theme_set(theme_bw())
rm(list = ls())

# Load old ones for checking and play with the indices
load("2020_files/runs/a4a/MCMC/input4MCMC.RData")
stk_2020 <- stk
stk_disc_2020 <- stk_disc
idx_bio
rm(stk);rm(stk_disc)

swo.stk <- readRDS("Robj/swo_stk.rds")
idx <- readRDS("Robj/swo_bio_idx.rds")
idx[['SIC_LL']] <- idx_bio[['SI_LL']]

stk <- window(swo.stk, start = 1985)

stk <- setPlusGroup(stk, 5) # from 2020 the distribution is different now
range(stk)[c('minfbar','maxfbar')] <- c(2,4) 

# 2020 settings
names(idx)
fmod2 <- ~s(year, k = 19) + s(age, k = 3)
srmod1 <- ~s(year, k = 17)
qmod <- list(~1,~s(year, k = 9), ~s(year, k = 3), ~1, ~s(year, k = 3),~1)

fit <- sca(stk, idx, fmodel = fmod2, srmodel = srmod1, qmodel = qmod) 
a4a.stk <- stk + simulate(fit, nsim = 1000) # numbers are crazy

plot(a4a.stk)

res <- residuals(fit, stk, idx)
plot(res)
