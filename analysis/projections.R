library(FLCore)
library(FLa4a)
library(FLasher)
library(FLBRP)

rm(list = ls())

swo.stk   <- readRDS("Robj/swo_stk.rds")
swo.stk.l <- readRDS("Robj/swo_stk_loren.rds")
swo.idx   <- readRDS("Robj/swo_bio_idx_std.rds")

## Use standardized longline indices (excluding index 2)
idx_use <- swo.idx[-2]
for (i in seq_along(idx_use)) {
  range(idx_use[[i]], c("min", "max")) <- c(1, 4)
  units(index(idx_use[[i]])) <- "t"
}

stk85_c <- replaceZeros(window(swo.stk, start = 1985))
range(stk85_c)[c('minfbar', 'maxfbar')] <- c(1, 4)
fmod_85 <- ~ s(year, k = 18) + s(age, k = 8) + ti(year, age, k = c(6, 5))
srmod   <- ~ s(year, k = 20)
qmod    <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3), ~1)

fit85_c <- sca(stk85_c, idx_use, fmodel = fmod_85, srmodel = srmod, qmodel = qmod)
a4a85_c <- stk85_c + simulate(fit85_c, nsim = 100)

model = "bevholt"

landings.n(a4a85_c)<-catch.n(a4a85_c)
landings.wt(a4a85_c)<-catch.wt(a4a85_c)
landings(a4a85_c)<-catch(a4a85_c)
discards(a4a85_c)<-0
discards.n(a4a85_c)<-0
discards.wt(a4a85_c)<-0

a4a_sr <- fmle(as.FLSR(a4a85_c, model = model),control=list(trace=0),fixed=list(b=100))

a4a_mtf <- fwdWindow(a4a85_c, end = 2040)

# Fix TACs 8kt, 7kt and 6kt
catch_df <- data.frame(year = rep(2024:2040,3), quant = "catch",
                       value = c(10434,10434,rep(6500,15),
                       10434,10434,rep(7000,15),
                       10434,10434,rep(7500,15),
                      10434,10434,rep(8000,15),
                      10434,10434,rep(8500,15),
                      10434,10434,rep(9000,15),
                      10434,10434,rep(9500,15),
                      10434,10434,rep(10000,15),
                      10434,10434,rep(10500,15),
                      10434,10434,rep(11000,15),
                      10434,10434,rep(11500,15),
                      10434,10434,rep(12000,15)
                      ),Scen = rep(1:12, each = 17))

ctrl_catch <- lapply(1:12, function(x){
  fwdControl(catch_df[catch_df$Scen == x,-4])
})

res00 <- residuals(a4a_sr)
rec00 <- window(rec(a4a_mtf), 2024, 2040)
rec00 <- rlnorm(rec00, mean(res00), sqrt(var(res00)))

catch_prj <- lapply(1:12, function(x){
  fwd(a4a_mtf, control = ctrl_catch[[x]], sr = a4a_sr, residuals = rec00)})


a4a_catch_8 <- fwd(a4a_mtf, control = ctrl_catch[[1]], sr = a4a_sr, residuals = rec00)
a4a_catch_7 <- fwd(a4a_mtf, control = ctrl_catch[[2]], sr = a4a_sr, residuals = rec00)
a4a_catch_6 <- fwd(a4a_mtf, control = ctrl_catch[[3]], sr = a4a_sr, residuals = rec00)

plot(FLStocks(`6000` = a4a_catch_6, `7000` = a4a_catch_7, `8000` = a4a_catch_8))
catch(a4a85_c)

as.data.frame(stock(a4a_catch_8))
