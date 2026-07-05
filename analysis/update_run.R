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
swo.idx <- readRDS("Robj/swo_bio_idx_std.rds")

### ------------------------------------------------------ ###
# Comparisons of data 2020 vs 2026
### ------------------------------------------------------ ###
# Total catch:
ggplot(as.data.frame(FLQuants(ca2020 = catch(stk_2020),
                              ca2026 = catch(swo.stk)))) +
  geom_line(aes(x = year, y = data, color = qname))

ggplot(as.data.frame(FLQuants(caa2020 = catch.n(stk_2020), 
                              caa2026 = catch.n(swo.stk)))) +
  geom_line(aes(x = year, y = data, group = qname, color = qname))+
  facet_wrap(~age, ncol = 3)

df <- do.call(rbind, lapply(names(swo.idx), function(nm) {
  y <- as.data.frame(swo.idx[[nm]])
  y <- y[y$slot == 'index', ]
  y$index <- nm
  y
}))
df$slot <- '2026 indices'

df2 <- do.call(rbind, lapply(names(idx_2020), function(nm) {
  y <- as.data.frame(idx_2020[[nm]])
  y <- y[y$slot == 'index', ]
  y$index <- nm
  y
}))
df2$slot <- '2020 indices'

df <- rbind(df,df2)
ggplot(df, aes(x = year, y = data, group = slot, color = slot)) + geom_line()+
  facet_wrap(~index, ncol = 3)

library(dplyr)
# Standardize again just for the plotting 
df <- df %>%
  group_by(index, slot) %>%
  mutate(data_std = data / mean(data, na.rm = TRUE)) %>%
  ungroup()

ggplot(df, aes(x = year, y = data_std, group = slot, color = slot)) + geom_line()+
  facet_wrap(~index, ncol = 3)

### ------------------------------------------------------ ###
# Run of 2020 (for comparison)
### ------------------------------------------------------ ###
# Dropping the Ligurian surface (as we did in 2020)
idx_2020 <- idx_2020[-6]
# NOTE: final a4a had discards included
fmod2020 <- ~s(year, k = 17) + s(age, k = 3)
srmod2020 <- ~s(year, k = 15)
# qmod2020 <- list(~1,~1, ~s(year, k = 3), ~s(year, k = 3),~1,~s(year, k = 9))
qmod2020 <- list(~1,~1, ~s(year, k = 3), ~s(year, k = 3),~1)
fit2020 <- sca(stk_2020, idx_2020, fmodel = fmod2020, srmodel = srmod2020, qmodel = qmod2020)
a4a.2020 <- stk_2020 + simulate(fit2020, nsim = 1000)

### ------------------------------------------------------ ###
# Update run 1985 - 2024 plus group 5
### ------------------------------------------------------ ###
stk <- window(swo.stk, start = 1985) # truncate early years as in 2020
stk <- replaceZeros(stk)
stk <- setPlusGroup(stk, 5) # same plus group
range(stk)[c('minfbar','maxfbar')] <- c(2,4) 

# Also dropping the Ligurian from here
idx <- swo.idx[-2]
names(idx)

fmod2 <- ~s(year, k = 20) + s(age, k = 3) # added ks - keep ratio
srmod1 <- ~s(year, k = 17) # added ks
# qmod <- list(~1,~s(year, k = 9), ~s(year, k = 3), ~1, ~s(year, k = 3),~1)
qmod <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3),~1)
fit <- sca(stk, idx, fmodel = fmod2, srmodel = srmod1, qmodel = qmod) 
a4a.stk <- stk + simulate(fit, nsim = 1000)

plot(FLStocks(update_run = a4a.stk, run_2020 = a4a.2020))
res <- residuals(fit, stk, idx)
res2020 <- residuals(fit2020, stk_2020, idx_2020)
plot(res2020)
plot(res, by = 'age')

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
# The assessment fails
### ------------------------------------------------------ ###
fit2 <- sca(stk, idx_2020, fmodel = fmod2020, srmodel = srmod2020, qmodel = qmod2020) 
a4a.sens2 <- stk + simulate(fit2, nsim = 1000)
plot(FLStocks(run_2020 = a4a.2020, run_sens2 = a4a.sens2))

# Plot all together
plot(FLStocks(run_2020 = a4a.2020, run_sens1 = a4a.sens1, run_sens2 = a4a.sens2, catch = stk))

### ------------------------------------------------------ ###
# Status - Ref points
### ------------------------------------------------------ ###

stk <- window(swo.stk, start = 1985)
stk <- replaceZeros(stk)
stk <- setPlusGroup(stk, 5) # same plus group
aa.stk <- stk + fit

landings.n(aa.stk)<-catch.n(aa.stk)
landings.wt(aa.stk)<-catch.wt(aa.stk)
landings(aa.stk)<-catch(aa.stk)
discards(aa.stk)<-0
discards.n(aa.stk)<-0
discards.wt(aa.stk)<-0

model = "bevholt"
srr_upd <- fmle(as.FLSR(aa.stk, model = model),control=list(trace=0),fixed=list(b=100))
brp_upd<- brp(FLBRP(aa.stk, sr = srr_upd))
ref_points_upd <- refpts(brp_upd)

aa.stk <- stk_2020 + fit2020
landings.n(aa.stk)<-catch.n(aa.stk)
landings.wt(aa.stk)<-catch.wt(aa.stk)
landings(aa.stk)<-catch(aa.stk)
discards(aa.stk)<-0
discards.n(aa.stk)<-0
discards.wt(aa.stk)<-0
srr_2020 <- fmle(as.FLSR(aa.stk, model = model),control=list(trace=0),fixed=list(b=100))
brp_2020 <- brp(FLBRP(aa.stk, sr = srr_2020))
ref_points_2020 <- refpts(brp_2020)


fmsy_upd <- as.numeric(ref_points_upd['msy','harvest'])
ssbmsy_upd <- as.numeric(ref_points_upd['msy','ssb'])
p1 <- plot(a4a.stk, metrics=list(SSB=ssb, F=fbar)) +
  geom_flpar(data=FLPars(SSB=FLPar(SSBmsy = ssbmsy_upd),
  F=FLPar(FMSY=fmsy_upd)), x=c(1990, 2015))

fmsy_upd <- as.numeric(ref_points_2020['msy','harvest'])
ssbmsy_upd <- as.numeric(ref_points_2020['msy','ssb'])
p2 <- plot(a4a.2020, metrics=list(SSB=ssb, F=fbar)) +
  geom_flpar(data=FLPars(SSB=FLPar(SSBmsy = ssbmsy_upd),
  F=FLPar(FMSY=fmsy_upd)), x=c(1990, 2015))

p1+p2
