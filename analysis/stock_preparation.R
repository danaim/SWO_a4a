library(FLCore)
library(FLa4a)
library(FLBRP)
library(ggplot2);theme_set(theme_bw())
rm(list = ls())

# Load old ones just for checking
load("2020_files/runs/a4a/MCMC/input4MCMC.RData")
stk_2020 <- stk
stk_disc_2020 <- stk_disc
rm(stk);rm(stk_disc)

caa <- read.csv("data/caa.csv")
caa.flq <- FLQuant(t(as.matrix(caa[,-1])), dimnames = list(age = 0:9,year = 1972:2024))
caa.flq <- caa.flq

ca <- read.csv("data/SWOMed_JABBA_catch_10April2026.csv")
ca.flq <- FLQuant(as.vector(ca[,-1]), dimnames = list(age = 'all',year = 1950:2024))

mat_0_9  <- c(0.00, 0.00, 0.15, 0.65, 1.00, 1.00, rep(1.00, 4))
mat.flq <- FLQuant(mat_0_9, dimnames = list(age = 0:9,year = 1972:2024))

m.flq <- FLQuant(0.2, dimnames = list(age = 0:9,year = 1972:2024))

cw <- read.csv("data/waa.csv")

# GT: Check if we need something else other than mean
cw[cw$Age == 9,]$Weight <- mean(cw[cw$Age %in% c(9,10),]$Weight) 
cw.v <- cw[cw$Age %in% c(0:9),]$Weight
cw.flq <- FLQuant(cw.v, dimnames = list(age = 0:9,year = 1972:2024))
cw.flq <- cw.flq/1000

stk <- FLStock(catch.n  = caa.flq)
stk@catch.wt = cw.flq
stk@catch <- ca.flq[, ac(1972:2024)]

stk@landings <- stk@catch
stk@landings.wt = cw.flq
stk@landings <- ca.flq[, ac(1972:2024)]

stk@stock.wt <- stk@catch.wt

stk@mat <- mat.flq
stk@m <- m.flq
stk@harvest.spwn <- stk@m.spwn <-  FLQuant(0.5, dimnames = list(age = 0:9,year = 1972:2024))

min_yr <- range(stk)['minyear']
max_yr <- range(stk)['maxyear']
min_age <- range(stk)['min']
max_age <- range(stk)['max']

# SoP GT: I put 1 in the NAs in catch.n just to do the SoP
catch.n(stk)[1, ac(1972:1977)] <- 1
SOP <- as.matrix(as.data.frame(catch(stk))[7]/colSums(catch.n(stk)*catch.wt(stk)))
SOP.flq <- FLQuant(rep(SOP,each=10),dimnames=list(age=0:max_age+1, year = (min_yr:max_yr)))
catch.n(stk) <- catch.n(stk)*SOP.flq
(SOP=as.matrix(as.data.frame(catch(stk))[7]/colSums(catch.n(stk)*catch.wt(stk))))
catch(stk)/computeCatch(stk)

catch.n(stk)[1, ac(1972:1977)] <- NA

units(stk) <- standardUnits(stk)
units(catch.wt(stk)) <- units(landings.wt(stk)) <- units(discards.wt(stk)) <- "t"
units(stock.wt(stk)) <- 't'
units(stk)

ggplot(data = catch.n(setPlusGroup(stk,5))) + geom_line(aes(x = age, y = data, group = year, color = year))
saveRDS(stk,"Robj/swo_stk.rds")

### Random stuff
# GT : Don't bother to look from now on
plot(stk)

catch.wt(setPlusGroup(stk,5))[,ac(1985:2018)]%/%catch.wt(stk_disc_2020)
stk.ext <- window(stk, start = 1950, end = 2024)
plot(catch(stk))

p1 <- ggplot(data = catch.n(setPlusGroup(stk,5))) + geom_line(aes(x = age, y = data, group = year, color = year))

p2 <- ggplot(data = catch.n(stk_disc_2020)) + geom_line(aes(x = age, y = data, group = year, color = year))

p1 + p2

