library(FLa4a)
library(FLBRP)
library(mse)
swo.stk <- readRDS("datasets/swo_stk.rds")
idx <- readRDS("datasets/swo_bio_idx.rds")
stk <- swo.stk

stk <- replaceZeros(stk)
stk <- window(stk, start = 1978)

# Stock assessment
fmod4 <- ~s(year, k = 30) + s(age, k = 5)
srmod3 <- ~I(as.numeric(year<1987)) + 1 + s(year, k = 15,by = as.numeric(year>=1987))
qmod <- list(~1,~s(year, k = 9), ~s(year, k = 3), ~1, ~s(year, k = 3),~1)
fit3 <- sca(stk, idx, fmodel = fmod4, srmodel = srmod3, qmodel = qmod)

a4aom.ctrl <- list(stkname = "Stock", dl = 1, npy = 15, sr = "segreg", 
blim_proxy = "spr.10", bsafe_proxy = "spr.30", 
fmsy_proxy = "f0.1", conditioning_ny = 3)

stk <- stk + simulate(fit3, nsim = 25)

a4aom <- function(stk, ctrl, ...){

  # setup
  spread(ctrl)
  it <- dim(stk)[6]
  iy <- range(stk)["maxyear"]
  dy <- iy - dl
  fy <- iy + npy

  # S/R
  om.sr <- fmle(as.FLSR(qapply(stk, iterMedians), model=sr), control = list(trace = 0))
  om.srdevs <- rlnormar1(it, sdlog=sd(residuals(om.sr)), years=seq(dy, fy))

  # BRP
  om.brp <- FLBRP(stk, sr=om.sr)
  rp <- refpts(om.brp)
  rpnms <- unique(c(dimnames(rp)$refpt, blim_proxy, bsafe_proxy, fmsy_proxy))

  if(length(rpnms)>7){
    rp <- rp[c(1:nrow(rp),rep(1, length(rpnms)-7))]
    dimnames(rp)$refpt <- rpnms
    refpts(om.brp) <- rp
  }
  om.brp <- brp(om.brp)
  rp <- remap(refpts(om.brp), map=list(BLIM=c(blim_proxy, "ssb"), BSAFE = c(bsafe_proxy, "ssb"), FMSY = c(fmsy_proxy, "harvest"), SBMSY = c("msy", "ssb"), BMSY = c("msy", "biomass"), B0 = c("virgin", "biomass"), SB0 = c("virgin", "ssb")))

  # OM
  om <- FLom(stock=stk, refpts=rp, model=sr, params=params(om.sr), deviances=om.srdevs, name=stkname)
  om <- fwdWindow(om, end=fy, nsq=conditioning_ny)
  om

}

swoom <- a4aom(stk, a4aom.ctrl)


