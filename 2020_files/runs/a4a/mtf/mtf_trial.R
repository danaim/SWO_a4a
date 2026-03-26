rm(list = ls())
load("~/pCloudDrive/ICCAT/runs/a4a/mtf/a4afit9.RData")

a4a_mtf <- stf(a4a.stk9, nyears = 12)

model = "bevholt"
a4a_sr <- fmle(as.FLSR(a4a.stk9, model = model),control=list(trace=0),fixed=list(b=100))

model = "bevholt"
a4a_sr <- fmle(as.FLSR(a4a.stk9, model = model),control=list(trace=0))

model = "geomean"
a4a_sr <- fmle(as.FLSR(a4a.stk9, model = model),control=list(trace=0))


# Fix TACs 8kt, 7kt and 6kt
catch_df <- data.frame(year = rep(2019:2030,3), quantity = "catch",
                       val = c(rep(8000,12),rep(7000,12),rep(6000,12)))
ctrl_catch <- list()

ctrl_catch[[1]] <- fwdControl(catch_df[catch_df$val == 8000,])
ctrl_catch[[2]] <- fwdControl(catch_df[catch_df$val == 7000,])
ctrl_catch[[3]] <- fwdControl(catch_df[catch_df$val == 6000,])

a4a_catch_8 <- fwd(a4a_mtf, ctrl_catch[[1]], sr = a4a_sr)
a4a_catch_7 <- fwd(a4a_mtf, ctrl_catch[[2]], sr = a4a_sr)
a4a_catch_6 <- fwd(a4a_mtf, ctrl_catch[[3]], sr = a4a_sr)

plot(FLStocks(catch_8kt = a4a_catch_8, catch_7kt = a4a_catch_7, catch_6kt = a4a_catch_6))
plot(window(a4a_catch_8, start=2000))

# stochastuc recruitment
niters <- 1000
a4a_mtf <- stf(a4a.stk9, nyears = 12)
a4a_mtf <- propagate(a4a_mtf, niters)

multi_rec_residuals <- FLQuant(NA, dimnames = list(year=2019:2030, iter=1:niters))
sample_years <- sample(dimnames(residuals(a4a_sr))$year, niters * 12, replace = TRUE)

multi_rec_residuals[] <- exp(residuals(a4a_sr)[,sample_years])

a4a_stoch_rec <- fwd(a4a_mtf, ctrl = ctrl_catch[[1]], sr = a4a_sr, 
                     sr.residuals = multi_rec_residuals, sr.residuals.mult = TRUE)

plot(window(a4a_stoch_rec, start = 2000))
