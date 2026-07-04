# Update run with the same settings
# Data are not even close to what they used to be
library(FLCore)
library(FLa4a)
library(FLBRP)
library(ggplot2);theme_set(theme_bw())
rm(list = ls())

### ------------------------------------------------------ ###
# load data
stk_2020 <- readRDS(file = "Robj/swo_stk_2020.rds")
idx_2020 <- readRDS(file = "Robj/swo_bio_idx_2020.rds")
swo.stk <- readRDS("Robj/swo_stk.rds")
swo.idx <- readRDS("Robj/swo_bio_idx.rds")

### ------------------------------------------------------ ###
# TODO: add plots and comparisons of data 2020 vs 2026
### ------------------------------------------------------ ###

### ------------------------------------------------------ ###
# Run of 2020 (for comparison)
### ------------------------------------------------------ ###

# NOTE: final a4a had discards included
fmod2020 <- ~s(year, k = 17) + s(age, k = 3)
srmod2020 <- ~s(year, k = 15)
qmod2020 <- list(~1,~1, ~s(year, k = 3), ~s(year, k = 3),~1,~s(year, k = 9))
fit2020 <- sca(stk_2020, idx_2020, fmodel = fmod2020, srmodel = srmod2020, qmodel = qmod2020)
a4a.2020 <- stk_2020 + simulate(fit2020, nsim = 1000)

### ------------------------------------------------------ ###
# Update run 1985 - 2024 plus group 5
### ------------------------------------------------------ ###
stk <- window(swo.stk, start = 1985) # truncate early years as in 2020
stk <- replaceZeros(stk)
stk <- setPlusGroup(stk, 5) # same plus group
range(stk)[c('minfbar','maxfbar')] <- c(2,4) 

idx <- swo.idx
names(idx)
fmod2 <- ~s(year, k = 19) + s(age, k = 3) # added ks
srmod1 <- ~s(year, k = 17) # added ks
qmod <- list(~1,~s(year, k = 9), ~s(year, k = 3), ~1, ~s(year, k = 3),~1)
fit <- sca(stk, idx, fmodel = fmod2, srmodel = srmod1, qmodel = qmod) 
a4a.stk <- stk + simulate(fit, nsim = 1000)

plot(FLStocks(update_run = a4a.stk, run_2020 = a4a.2020))
res <- residuals(fit, stk, idx)
plot(res)

### ------------------------------------------------------ ###
# Sensitivity runs
# Truncate data up to 2018
### ------------------------------------------------------ ###
stk <- window(stk, end = 2018)
idx <- window(idx, end = 2018)
fit1 <- sca(stk, idx, fmodel = fmod2020, srmodel = srmod2020, qmodel = qmod) 
a4a.sens1 <- stk + simulate(fit1, nsim = 1000)
plot(FLStocks(run_2020 = a4a.2020, run_sens1 = a4a.sens1))

### ------------------------------------------------------ ###
# Sensitivity runs
# Truncate data up to 2018 - use old indices
### ------------------------------------------------------ ###
# GT: check this one, it fails completely
fit2 <- sca(stk, idx_2020, fmodel = fmod2020, srmodel = srmod2020, qmodel = qmod2020) 
a4a.sens2 <- stk + simulate(fit2, nsim = 1000)
plot(FLStocks(run_2020 = a4a.2020, run_sens2 = a4a.sens2))

