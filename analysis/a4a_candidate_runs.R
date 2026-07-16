library(FLCore)
library(FLa4a)
library(FLBRP)
library(ggplot2);theme_set(theme_bw())
rm(list = ls())

### ------------------------------------------------------ ###
# load data
swo.stk <- readRDS("Robj/swo_stk.rds")
swo.stk.l <- readRDS("Robj/swo_stk_loren.rds")
swo.idx <- readRDS("Robj/swo_bio_idx_std.rds")
swo.idx.ext <- readRDS("Robj/swo_bio_idx_std_ext.rds")

# Only for use with the swo.idx with the swo.idx.ext we are ok
idx_use <- swo.idx[-2]
for(i in 1:length(swo.idx)) {
  range(swo.idx[[i]],c("min","max")) = c(1,4)
  units(index(swo.idx[[i]])) <- "t"
}

## Run 1: 1985 Constant M
stk85_c <- replaceZeros(window(swo.stk, start = 1985))
range(stk85_c)[c('minfbar', 'maxfbar')] <- c(1, 4)
fmod_85 <- ~ s(year, k = 20) + s(age, k = 8) + ti(year, age, k = c(6, 5))
srmod   <- ~ s(year, k = 20)
qmod    <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3), ~1)

fit85_c <- sca(stk85_c, idx_use, fmodel = fmod_85, srmodel = srmod, qmodel = qmod)
a4a85_c <- stk85_c + simulate(fit85_c, nsim = 1000)

## Run 2: 1985 Lorenzen M
stk85_l <- replaceZeros(window(swo.stk.l, start = 1985))
range(stk85_l)[c('minfbar', 'maxfbar')] <- c(1, 4)
fit85_l <- sca(stk85_l, idx_use, fmodel = fmod_85, srmodel = srmod, qmodel = qmod)
a4a85_l <- stk85_l + simulate(fit85_l, nsim = 1000)

## Sensitivity 1: 1978 Constant M
stk78_c <- replaceZeros(window(swo.stk, start = 1978))
range(stk78_c)[c('minfbar', 'maxfbar')] <- c(1, 4)
fmod_78 <- ~ s(year, k = 25) + s(age, k = 8) + ti(year, age, k = c(6, 5))

fit78_c <- sca(stk78_c, idx_use, fmodel = fmod_78, srmodel = srmod, qmodel = qmod)
a4a78_c <- stk78_c + simulate(fit78_c, nsim = 1000)

## Sensitivity 2: 1985 Constant M - extra indices
stk_s2 <- replaceZeros(window(swo.stk.l, start = 1985))
idx_s2 <- readRDS("Robj/swo_bio_idx_std_ext.rds")
range(stk_s2)[c('minfbar', 'maxfbar')] <- c(1, 4)
qmod_s2    <- list(~1, ~s(year, k = 5), ~s(year, k = 3), ~1, ~s(year, k = 3), ~1, ~1)
fit_s2 <- sca(stk_s2, idx_s2, fmodel = fmod_78, srmodel = srmod, qmodel = qmod_s2)
a4a_s2 <- stk_s2 + simulate(fit_s2, nsim = 1000)

plot(FLStocks(
  `85 Const M`  = a4a85_c,
  `85 Loren M`  = a4a85_l,
  `Observed catch` = stk85_c
)) + ylim(0.01, NA) + theme(legend.position="top")

# for constant M
res_df <- as.data.frame(residuals(fit85_c, stk85_c, idx_use))
plot(res_df) 

## Proposal of the group: Exclude Sicilian LL - use only four indices
idx_gp <- idx_use[-5]
stk85_c <- replaceZeros(window(swo.stk, start = 1985))
range(stk85_c)[c('minfbar', 'maxfbar')] <- c(1, 4)
fmod_85 <- ~ s(year, k = 20) + s(age, k = 8) + ti(year, age, k = c(6, 5))
srmod   <- ~ s(year, k = 20)
qmod_gp    <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3))
fit_gp <- sca(stk85_c, idx_gp, fmodel = fmod_85, srmodel = srmod, qmodel = qmod_gp)
a4a_gp <- stk85_c + simulate(fit_gp, nsim = 1000)

plot(FLStocks(`4 indices` = a4a_gp, `Constant M` = a4a85_c))

## Try f mode with k = 17 n1 model
range(stk85_c)[c('minfbar', 'maxfbar')] <- c(0, 9)
fmod_85_17 <- ~s(year, k = 17) + s(age, k = 8) + ti(year, age, k = c(6, 5))
fit85_c_17 <- sca(stk85_c, idx_use, fmodel = fmod_85_17, srmodel = srmod, qmodel = qmod)
a4a85_c_17 <- stk85_c + simulate(fit85_c, nsim = 1000)

wireframe(harvest(fit85_c_17))
plot_status(stk85_c, fit85_c_17, a4a85_c_17)+ geom_vline(xintercept = 1988, linetype = 'dashed')

plot_status(stk85_c, fit85_c, a4a85_c) + geom_vline(xintercept = 1988, linetype = 'dashed')

plot(FLStocks(`Constant M 17` = a4a85_c_17, `Constant M` = a4a85_c))
res_df_17 <- as.data.frame(residuals(fit85_c_17, stk85_c, idx_use))
res_df


## Sensitivity 1: 1972 Constant M
stk72_c <- replaceZeros(window(swo.stk, start = 1972))
catch.n(stk72_c)[1,ac(c(1972:1977))] <- 0
stk72_c <- replaceZeros(stk72_c)

range(stk72_c)[c('minfbar', 'maxfbar')] <- c(1, 4)
fmod_72 <- ~ s(year, k = 25) + s(age, k = 8) + ti(year, age, k = c(6, 5))
fit72_c <- sca(stk72_c, idx_use, fmodel = fmod_72, srmodel = srmod, qmodel = qmod)
a4a72_c <- stk72_c + simulate(fit72_c, nsim = 1000)

source("analysis/plot_functions.R")
plot(FLStocks(`72` = a4a72_c, `78` = a4a78_c))
plot_status(stk72_c, fit72_c, a4a72_c)

### ===================================================================== ###
### Final run ----------------------------------------------------------- ###
### ===================================================================== ###

## Proposal of the group: Exclude Sicilian LL - use only four indices
idx_gp <- idx_use[-5]
stk85_c <- replaceZeros(window(swo.stk, start = 1985))
range(stk85_c)[c('minfbar', 'maxfbar')] <- c(1, 4)
fmod_fn <- ~ s(year, k = 18) + s(age, k = 8) + ti(year, age, k = c(6, 5))
srmod   <- ~ s(year, k = 20)
qmod_gp    <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3))
fit_fn <- sca(stk85_c, idx_gp, fmodel = fmod_fn, srmodel = srmod, qmodel = qmod_gp)
a4a_fn <- stk85_c + simulate(fit_fn, nsim = 1000)
plot(a4a_fn)
