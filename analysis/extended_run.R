# Update run with the same settings
# Data are not even close to what they used to be
library(FLCore)
library(FLa4a)
library(FLBRP)
library(ggplot2);theme_set(theme_bw())
rm(list = ls())

ca <- read.csv("data/SWOMed_JABBA_catch_10April2026.csv")
ca.flq <- FLQuant(as.vector(ca[,-1]), dimnames = list(age = 'all',year = 1950:2024))

# Load old ones for checking and play with the indices
load("2020_files/runs/a4a/MCMC/input4MCMC.RData")
stk_2020 <- stk
stk_disc_2020 <- stk_disc
idx_2020 <- idx_bio
rm(stk);rm(stk_disc)
fmod2020 <- ~s(year, k = 17) + s(age, k = 3)
srmod2020 <- ~s(year, k = 15)
qmod2020 <- list(~1,~1, ~s(year, k = 3), ~s(year, k = 3),~1,~s(year, k = 9))
fit2020 <- sca(stk_disc_2020, idx_2020, fmodel = fmod2020, srmodel = srmod2020, qmodel = qmod2020)
a4a.stk2020 <- stk_disc_2020 + simulate(fit2020, nsim = 1000)


swo.stk <- readRDS("Robj/swo_stk.rds")
idx <- readRDS("Robj/swo_bio_idx.rds")
idx[['SIC_LL']] <- idx_bio[['SI_LL']]
rm(idx_bio)

stk <- swo.stk
stk <- replaceZeros(stk)
fmod2 <- ~s(year, k = 20) + s(age, k = 4)
fmod3 <- ~s(year, k = 25) + s(age, k = 4) + te(year, age, k = c(8,4))
fmod4 <- ~s(year, k = 30) + s(age, k = 5)


## different sr models
srmod1 <- ~s(year, k = 17)
srmod3 <- ~I(as.numeric(year<2000)) - 1 + s(year, k = 15,by = as.numeric(year>=1987))
srmod7 <- ~ I(as.numeric(year<1987)) + 1 + s(year, k=15, by=as.numeric(year>=1987))
srmod8 <- ~I(as.numeric(year>=1987)) + 1 + s(year, k = 15,by = as.numeric(year>1987))
srmod6 <- ~s(year, k = 15,by = as.numeric(year>=1987))
srmod4 <- ~s(year, k = 15, by = as.numeric(year, year>1987))
srmod5 <- ~bevholt(CV = 0.2)

qmod <- list(~1,~s(year, k = 9), ~s(year, k = 3), ~1, ~s(year, k = 3),~1)

# srmods <- list(srmod1, srmod2, srmod3, srmod4 , srmod5)

# extended run up to 1972
stk2 <- window(stk, start = 1978)
# fits <- scas(FLStocks(stk), list(idx), fmodel = fmod2, srmodel = srmods, qmodel = list (qmod))

fit1 <- sca(stk2, idx, fmodel = fmod2, srmodel = srmod1, qmodel = qmod) 
fit3 <- sca(stk2, idx, fmodel = fmod2, srmodel = srmod3, qmodel = qmod)
fit7 <- sca(stk2, idx, fmodel = fmod2, srmodel = srmod7, qmodel = qmod)
fit8 <- sca(stk2, idx, fmodel = fmod2, srmodel = srmod8, qmodel = qmod)
fit6 <- sca(stk2, idx, fmodel = fmod2, srmodel = srmod6, qmodel = qmod)
fit4 <- sca(stk2, idx, fmodel = fmod2, srmodel = srmod4, qmodel = qmod)
fit5 <- sca(stk2, idx, fmodel = fmod2, srmodel = srmod5, qmodel = qmod)

a4a.stk1 <- stk2 + simulate(fit1, nsim = 1000)
a4a.stk3 <- stk2 + simulate(fit3, nsim = 1000)
a4a.stk7 <- stk2 + simulate(fit7, nsim = 1000)
a4a.stk8 <- stk2 + simulate(fit8, nsim = 1000)
a4a.stk6 <- stk2 + simulate(fit6, nsim = 1000)
a4a.stk4 <- stk2 + simulate(fit4, nsim = 1000)
a4a.stk5 <- stk2 + simulate(fit5, nsim = 1000)

plot(FLStocks(fit5 = a4a.stk5, fit4 = a4a.stk4, fit7 = a4a.stk7, fit6 = a4a.stk6,
fit3 = a4a.stk3))
plot(FLStocks(fit3 = a4a.stk3, fit6 = a4a.stk6))
plot(FLStocks(baserun = a4a.stk2, run2020 = a4a.stk2020, 
    catch = stk, fit3 = a4a.stk3, fit4 = a4a.stk4))

plot(FLStocks(baserun = a4a.stk2, run2020 = a4a.stk2020, 
    catch = stk))

### --------------------------------------------------- ###
### This should be done properly with FLife and simulating 
### better a population in an almost virgin state
### for now we just taking a random CFD to put in the 
### beginning of the time series and we scale it to the
### catch. 
### -------------------------------------------------- ###

# Trial with an assumed initial population and total catches
stk_ext <- window(stk, start = 1950)
stock.wt(stk_ext) <- catch.wt(stk_ext) <- as.numeric(catch.wt(stk)[,1])
harvest.spwn(stk_ext) <- m.spwn(stk_ext) <- 0.5
catch.n(stk_ext)[,1:10] <- catch.n(stk_ext)[,ac(1980)]
catch(stk_ext) <- ca.flq

# Scale everything
xc <- catch(stk_ext)[,1:10]
xca <- catch.n(stk_ext)[,1:10]
xwt <- catch.wt(stk_ext)[,1:10]

SOP <- as.matrix(as.data.frame(xc)[7]/colSums(xca*xwt))
SOP.flq <- FLQuant(rep(SOP,each=10),dimnames=list(age=0:9, year = (1950:1959)))
catch.n(stk_ext)[,1:10] <- catch.n(stk_ext)[,1:10]*SOP.flq
# -- end of scaling


stk <- stk_ext
stk <- replaceZeros(stk_ext)
fmod2 <- ~s(year, k = 20) + s(age, k = 4)
fmod3 <- ~s(year, k = 25) + s(age, k = 4) + te(year, age, k = c(8,4))
fmod4 <- ~s(year, k = 30) + s(age, k = 5)


## different sr models
srmod1 <- ~s(year, k = 17)
srmod2 <- ~factor(replace(year, year<1988, 1988))
srmod3 <- ~s(year, k = 15, by = as.numeric(year, year>1987))

qmod <- list(~1,~s(year, k = 9), ~s(year, k = 3), ~1, ~s(year, k = 3),~1)

# extended run up to 1972
fit1 <- sca(stk, idx, fmodel = fmod2, srmodel = srmod1, qmodel = qmod)
fit2 <- sca(stk_ext, idx, fmodel = fmod2, srmodel = srmod2, qmodel = qmod) 
fit3 <- sca(stk_ext, idx, fmodel = fmod3, srmodel = srmod1, qmodel = qmod)
fit4 <- sca(stk, idx, fmodel = fmod4, srmodel = srmod1, qmodel = qmod)
plot(stk_ext)
plot(catch.n(stk_ext))

xx <- as.data.frame(idx)
ggplot(data = xx[xx$slot == 'index',]) + geom_line(aes(x = year, y = data)) + facet_wrap(~cname)
