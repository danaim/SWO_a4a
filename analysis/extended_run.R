# Update run with the same settings
# Data are not even close to what they used to be
library(FLCore)
library(FLa4a)
library(FLBRP)
library(ggplot2);theme_set(theme_bw())
rm(list = ls())

# Load the old one to import the Sicilian index
load("2020_files/runs/a4a/MCMC/input4MCMC.RData")
rm(stk);rm(stk_disc)

swo.stk <- readRDS("Robj/swo_stk.rds")
idx <- readRDS("Robj/swo_bio_idx.rds")
idx[['SIC_LL']] <- idx_bio[['SI_LL']]
rm(idx_bio)

stk <- swo.stk
stk <- replaceZeros(stk)


### ------------------------------------------------------------------ ###
# fmodels
fmod2 <- ~s(year, k = 20) + s(age, k = 4)
fmod3 <- ~s(year, k = 25) + s(age, k = 4) + te(year, age, k = c(8,4))
fmod4 <- ~s(year, k = 30) + s(age, k = 5)

fmods <- list(fmod2, fmod3)
### ------------------------------------------------------------------ ###

### ------------------------------------------------------------------ ###
## different sr models
srmod1 <- ~s(year, k = 17)
srmod3 <- ~I(as.numeric(year<2000)) - 1 + s(year, k = 15, by = as.numeric(year>=1987))
srmod9 <- ~I(as.numeric(year<2000)) + 1 + s(year, k = 15, by = as.numeric(year>=1987))
srmod10 <- ~I(as.numeric(year<2000)) + 1 + s(year, k = 15, by = as.numeric(year>=1987),bs='cr' )
srmod7 <- ~I(as.numeric(year<1987)) + 1 + s(year, k = 15, by = as.numeric(year>=1987))
srmod8 <- ~I(as.numeric(year=<2024)) + 1 + s(year, k = 15,by = as.numeric(year>=1987))
srmod6 <- ~s(year, k = 15,by = as.numeric(year>=1987))
srmod4 <- ~s(year, k = 15, by = as.numeric(year, year>1987))
srmod5 <- ~bevholt(CV = 0.2)

srmods <- list(srmod1, srmod3, srmod7)
### ------------------------------------------------------------------ ###

### ------------------------------------------------------------------ ###
qmod1 <- list(~1,~s(year, k = 9), ~s(year, k = 3), ~1, ~s(year, k = 3),~1)
qmod2 <- list( ~1, ~s(year, k = 4), ~s(year, k = 3), ~1, ~s(year, k = 3), ~1 )

qmods <- list(qmod1, qmod2)
### ------------------------------------------------------------------ ###


# extended run up to 1978
stk2 <- window(stk, start = 1978)  # NAs the years 1972 - 1978

fits <- scas(FLStocks(stk2), list(idx), fmodel = fmods, 
    srmodel = srmods, qmodel = qmods,
    combination.all = TRUE)

plot(fits) # fails for now

stks <- lapply(fits, function(x) stk2 + simulate(x, nsim = 500))
ress <- lapply(fits, function(x){
    residuals(x, stk2, idx)
})

plot(ress[[6]])
wireframe(harvest(fits[[1]]))

### ------------------------------------------------- ###
### Random things
### ------------------------------------------------- ###

df <- data.frame(year = 1978:2024)

dm3 <- getX(srmod3, df)
dm9 <- getX(srmod9, df)
dm10 <- getX(srmod10, df)
ncol(dm3);ncol(dm9)


par(mfrow = c(1, 2))

image(x = df$year, y = 1:ncol(dm3), z = dm3, 
      main = "Design Matrix: srmod3 (- 1)", 
      xlab = "Year", ylab = "Basis Function",
      col = hcl.colors(50, "Blue-Red"))


image(x = df$year, y = 1:ncol(dm9), z = dm9, 
      main = "Design Matrix: srmod9 (+ 1)", 
      xlab = "Year", ylab = "Basis Function",
      col = hcl.colors(50, "Blue-Red"))
par(mfrow = c(1, 1))


par(mfrow = c(1, 2))
matplot(x = df$year, y = dm3[,1:2], type = "l", lty = 1, lwd = 2,
        col = hcl.colors(ncol(dm3), "Spectral"),
        main = "Basis Functions: srmod3 (- 1)",
        xlab = "Year", ylab = "Basis Function")

abline(h = 0, col = "black", lwd = 2, lty = 2)
matplot(x = df$year, y = dm9[,1:2], type = "l", lty = 1, lwd = 2,
        col = hcl.colors(ncol(dm9), "Spectral"),
        main = "Basis Functions: srmod9 (+ 1)",
        xlab = "Year", ylab = "Basis Function")
abline(h = 0, col = "black", lwd = 2, lty = 2)
par(mfrow = c(1, 1))

par(mfrow = c(1, 2))
matplot(x = df$year, y = dm3, type = "l", lty = 1, lwd = 2,
        col = hcl.colors(ncol(dm3), "Spectral"),
        main = "Basis Functions: srmod3 (- 1)",
        xlab = "Year", ylab = "Basis Function")

abline(h = 0, col = "black", lwd = 2, lty = 2)
matplot(x = df$year, y = dm9, type = "l", lty = 1, lwd = 2,
        col = hcl.colors(ncol(dm9), "Spectral"),
        main = "Basis Functions: srmod9 (+ 1)",
        xlab = "Year", ylab = "Basis Function")
abline(h = 0, col = "black", lwd = 2, lty = 2)
par(mfrow = c(1, 1))


matplot(x = df$year, y = dm10, type = "l", lty = 1, lwd = 2,
        col = hcl.colors(ncol(dm10), "Spectral"),
        main = "Basis Functions: srmod9 (+ 1)",
        xlab = "Year", ylab = "Parameter Value")

image(x = df$year, y = 1:ncol(dm10), z = dm10, 
      main = "Design Matrix: srmod9 (+ 1)", 
      xlab = "Year", ylab = "Parameter / Basis Function",
      col = hcl.colors(50, "Blue-Red"))


### --------------------------------------------------- ###
### This should be done properly with FLife and simulating 
### better a population in an almost virgin state
### for now we just taking a random CFD to put in the 
### beginning of the time series and we scale it to the
### catch. 
### -------------------------------------------------- ###

ca <- read.csv("data/SWOMed_JABBA_catch_10April2026.csv")
ca.flq <- FLQuant(as.vector(ca[,-1]), dimnames = list(age = 'all',year = 1950:2024))

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


