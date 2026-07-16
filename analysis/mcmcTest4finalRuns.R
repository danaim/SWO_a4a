rm(list = ls())
library(FLa4a)
library(ggplotFL)
load("Robj/input4mcmc.Rdata")
ls()
aa.stk <- stk85_c
aa.idx <- idx_use
#====================================================================
# VARIABLES
#====================================================================

it <- 2000 # iterations
brn <- 500
mcsave <- 1000
mcmc <- mcsave*it
fmod <- ~s(year, k = 20) + s(age, k = 8) + ti(year, age, k = c(6, 5))
srmod <- ~s(year, k = 20)
qmod <- list(~1, ~s(year, k = 3),~1,~s(year, k = 3), ~1)
# n1mod <-  ~factor(age)
# vmod <- list(~1,~s(age,k=3))
alpha <- 0.2
vv <- seq(0.5,0, -0.1)
ssbcov <- fcov <- matrix(NA, ncol=length(vv), nrow=it-brn)
colnames(ssbcov) <- colnames(fcov) <- vv

#====================================================================
# runs
#====================================================================
#--------------------------------------------------------------------
# mrprobe 0.5
#--------------------------------------------------------------------
fit50 <- sca(aa.stk, aa.idx, fmodel=fmod, qmodel=qmod,
             srmodel=srmod,
             fit="MCMC", mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.499))
fit <- burnin(fit50, brn)
stk <- aa.stk + fit

fitmc <- fit
fitmc@pars@stkmodel@vcov[] <- cov(t(coef(stkmodel(pars(fit)))@.Data))
fitmc@pars@qmodel[[1]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[1]])@.Data))
fitmc@pars@qmodel[[2]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[2]])@.Data))
fitmc@pars@qmodel[[3]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[3]])@.Data))
fitmc@pars@qmodel[[4]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[4]])@.Data))
fitmc@pars@qmodel[[5]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[5]])@.Data))

fitmc@pars@vmodel[[1]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[1]])@.Data))
fitmc@pars@vmodel[[2]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[2]])@.Data))
fitmc@pars@vmodel[[3]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[3]])@.Data))
fitmc@pars@vmodel[[4]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[4]])@.Data))
fitmc@pars@vmodel[[5]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[5]])@.Data))
fitmc@pars@vmodel[[6]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[6]])@.Data))

mus <- median(ssb(stk)[,"2024"])
muf <- median(fbar(stk)[,"2024"])

for(i in 1:(it-brn)){
	fit. <- iter(fitmc, i)
	stk. <- aa.stk + simulate(fit., (it-brn))
	p <- ssb(stk.)[,"2024"]
	ssbcov[i,1] <- mus >= quantile(c(p), alpha/2) & mus <= quantile(c(p), 1-alpha/2)
	p <- fbar(stk.)[,"2024"]
	fcov[i,1] <- muf >= quantile(c(p), alpha/2) & muf <= quantile(c(p), 1-alpha/2)
}

#--------------------------------------------------------------------
# mrprobe 0.4
#--------------------------------------------------------------------
fit40 <- sca(aa.stk, aa.idx, fmodel=fmod, qmodel=qmod,
	srmodel=srmod,
	fit="MCMC", mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.4))
fit <- burnin(fit40, brn)
stk <- aa.stk + fit

fitmc <- fit
fitmc@pars@stkmodel@vcov[] <- cov(t(coef(stkmodel(pars(fit)))@.Data))
fitmc@pars@qmodel[[1]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[1]])@.Data))
fitmc@pars@qmodel[[2]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[2]])@.Data))
fitmc@pars@qmodel[[3]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[3]])@.Data))
fitmc@pars@qmodel[[4]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[4]])@.Data))
fitmc@pars@qmodel[[5]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[5]])@.Data))

fitmc@pars@vmodel[[1]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[1]])@.Data))
fitmc@pars@vmodel[[2]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[2]])@.Data))
fitmc@pars@vmodel[[3]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[3]])@.Data))
fitmc@pars@vmodel[[4]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[4]])@.Data))
fitmc@pars@vmodel[[5]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[5]])@.Data))
fitmc@pars@vmodel[[6]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[6]])@.Data))

mus <- median(ssb(stk)[,"2024"])
muf <- median(fbar(stk)[,"2024"])

for(i in 1:(it-brn)){
	fit. <- iter(fitmc, i)
	stk. <- aa.stk + simulate(fit., (it-brn))
	p <- ssb(stk.)[,"2024"]
	ssbcov[i,2] <- mus >= quantile(c(p), alpha/2) & mus <= quantile(c(p), 1-alpha/2)
	p <- fbar(stk.)[,"2024"]
	fcov[i,2] <- muf >= quantile(c(p), alpha/2) & muf <= quantile(c(p), 1-alpha/2)
}

#--------------------------------------------------------------------
# mrprobe 0.3
#--------------------------------------------------------------------
fit30 <- sca(aa.stk, aa.idx, fmodel=fmod, qmodel=qmod,
	srmodel=srmod,
	fit="MCMC", mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.3))
fit <- burnin(fit30, brn)
stk <- aa.stk + fit

fitmc <- fit
fitmc@pars@stkmodel@vcov[] <- cov(t(coef(stkmodel(pars(fit)))@.Data))
fitmc@pars@qmodel[[1]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[1]])@.Data))
fitmc@pars@qmodel[[2]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[2]])@.Data))
fitmc@pars@qmodel[[3]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[3]])@.Data))
fitmc@pars@qmodel[[4]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[4]])@.Data))
fitmc@pars@qmodel[[5]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[5]])@.Data))

fitmc@pars@vmodel[[1]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[1]])@.Data))
fitmc@pars@vmodel[[2]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[2]])@.Data))
fitmc@pars@vmodel[[3]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[3]])@.Data))
fitmc@pars@vmodel[[4]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[4]])@.Data))
fitmc@pars@vmodel[[5]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[5]])@.Data))
fitmc@pars@vmodel[[6]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[6]])@.Data))

mus <- median(ssb(stk)[,"2024"])
muf <- median(fbar(stk)[,"2024"])

for(i in 1:(it-brn)){
	fit. <- iter(fitmc, i)
	stk. <- aa.stk + simulate(fit., (it-brn))
	p <- ssb(stk.)[,"2024"]
	ssbcov[i,3] <- mus >= quantile(c(p), alpha/2) & mus <= quantile(c(p), 1-alpha/2)
	p <- fbar(stk.)[,"2024"]
	fcov[i,3] <- muf >= quantile(c(p), alpha/2) & muf <= quantile(c(p), 1-alpha/2)
}

#--------------------------------------------------------------------
# mrprobe 0.2
#--------------------------------------------------------------------
fit20 <- sca(aa.stk, aa.idx, fmodel=fmod, qmodel=qmod,
	srmodel=srmod,
	fit="MCMC", mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.2))
fit <- burnin(fit20, brn)
stk <- aa.stk + fit

fitmc <- fit
fitmc@pars@stkmodel@vcov[] <- cov(t(coef(stkmodel(pars(fit)))@.Data))
fitmc@pars@qmodel[[1]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[1]])@.Data))
fitmc@pars@qmodel[[2]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[2]])@.Data))
fitmc@pars@qmodel[[3]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[3]])@.Data))
fitmc@pars@qmodel[[4]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[4]])@.Data))
fitmc@pars@qmodel[[5]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[5]])@.Data))

fitmc@pars@vmodel[[1]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[1]])@.Data))
fitmc@pars@vmodel[[2]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[2]])@.Data))
fitmc@pars@vmodel[[3]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[3]])@.Data))
fitmc@pars@vmodel[[4]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[4]])@.Data))
fitmc@pars@vmodel[[5]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[5]])@.Data))
fitmc@pars@vmodel[[6]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[6]])@.Data))

mus <- median(ssb(stk)[,"2024"])
muf <- median(fbar(stk)[,"2024"])

for(i in 1:(it-brn)){
	fit. <- iter(fitmc, i)
	stk. <- aa.stk + simulate(fit., (it-brn))
	p <- ssb(stk.)[,"2024"]
	ssbcov[i,4] <- mus >= quantile(c(p), alpha/2) & mus <= quantile(c(p), 1-alpha/2)
	p <- fbar(stk.)[,"2024"]
	fcov[i,4] <- muf >= quantile(c(p), alpha/2) & muf <= quantile(c(p), 1-alpha/2)
}

#--------------------------------------------------------------------
# mrprobe 0.1
#--------------------------------------------------------------------
fit10 <- sca(aa.stk, aa.idx, fmodel=fmod, qmodel=qmod,
	srmodel=srmod,
	fit="MCMC", mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.1))
fit <- burnin(fit10, brn)
stk <- aa.stk + fit

fitmc <- fit
fitmc@pars@stkmodel@vcov[] <- cov(t(coef(stkmodel(pars(fit)))@.Data))
fitmc@pars@qmodel[[1]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[1]])@.Data))
fitmc@pars@qmodel[[2]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[2]])@.Data))
fitmc@pars@qmodel[[3]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[3]])@.Data))
fitmc@pars@qmodel[[4]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[4]])@.Data))
fitmc@pars@qmodel[[5]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[5]])@.Data))

fitmc@pars@vmodel[[1]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[1]])@.Data))
fitmc@pars@vmodel[[2]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[2]])@.Data))
fitmc@pars@vmodel[[3]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[3]])@.Data))
fitmc@pars@vmodel[[4]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[4]])@.Data))
fitmc@pars@vmodel[[5]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[5]])@.Data))
fitmc@pars@vmodel[[6]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[6]])@.Data))

mus <- median(ssb(stk)[,"2024"])
muf <- median(fbar(stk)[,"2024"])

for(i in 1:(it-brn)){
	fit. <- iter(fitmc, i)
	stk. <- aa.stk + simulate(fit., (it-brn))
	p <- ssb(stk.)[,"2024"]
	ssbcov[i,5] <- mus >= quantile(c(p), alpha/2) & mus <= quantile(c(p), 1-alpha/2)
	p <- fbar(stk.)[,"2024"]
	fcov[i,5] <- muf >= quantile(c(p), alpha/2) & muf <= quantile(c(p), 1-alpha/2)
}

#--------------------------------------------------------------------
# mrprobe 0
#--------------------------------------------------------------------
fit00 <- sca(aa.stk, aa.idx, fmodel=fmod, qmodel=qmod,
	srmodel=srmod,
	fit="MCMC", mcmc = SCAMCMC(mcmc = mcmc, mcsave = mcsave, mcprobe = 0.0001))
fit <- burnin(fit00, brn)
stk <- aa.stk + fit

fitmc <- fit
fitmc@pars@stkmodel@vcov[] <- cov(t(coef(stkmodel(pars(fit)))@.Data))
fitmc@pars@qmodel[[1]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[1]])@.Data))
fitmc@pars@qmodel[[2]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[2]])@.Data))
fitmc@pars@qmodel[[3]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[3]])@.Data))
fitmc@pars@qmodel[[4]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[4]])@.Data))
fitmc@pars@qmodel[[5]]@vcov[] <- cov(t(coef(qmodel(pars(fit))[[5]])@.Data))

fitmc@pars@vmodel[[1]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[1]])@.Data))
fitmc@pars@vmodel[[2]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[2]])@.Data))
fitmc@pars@vmodel[[3]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[3]])@.Data))
fitmc@pars@vmodel[[4]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[4]])@.Data))
fitmc@pars@vmodel[[5]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[5]])@.Data))
fitmc@pars@vmodel[[6]]@vcov[] <- cov(t(coef(vmodel(pars(fit))[[6]])@.Data))

mus <- median(ssb(stk)[,"2024"])
muf <- median(fbar(stk)[,"2024"])

for(i in 1:(it-brn)){
	fit. <- iter(fitmc, i)
	stk. <- aa.stk + simulate(fit., (it-brn))
	p <- ssb(stk.)[,"2024"]
	ssbcov[i,6] <- mus >= quantile(c(p), alpha/2) & mus <= quantile(c(p), 1-alpha/2)
	p <- fbar(stk.)[,"2024"]
	fcov[i,6] <- muf >= quantile(c(p), alpha/2) & muf <= quantile(c(p), 1-alpha/2)
}

#====================================================================
# process results
#====================================================================
fits <- a4aFitSAs(f50=fit50, f40=fit40, f30=fit30, f20=fit20, f10=fit10, f00=fit00)
mcmctest.df <- data.frame(mcprobe=vv, ssbcov=colMeans(ssbcov), fcov=colMeans(fcov), accrate = unlist(lapply(fits, function(x) fitSumm(x)["accrate",])))

# NOTE: 0.4 has the best compromise between CI coverage for both f 
#	and ssb and acceptance rate (0.3-0.4)

stks <- FLStocks('0.5'=aa.stk + fit50, '0.4'=aa.stk + fit40, '0.3'=aa.stk + fit30, '0.2'=aa.stk + fit20, '0.1'=aa.stk + fit10, '0.0'=aa.stk + fit00)

plot(stks)
plot(aa.stk + fit40)

## best for constM final mcprob = 0.4, accrate= 0.31
save(aa.stk,fit50,fit40,fit30,fit20,fit10,fit00,mcmctest.df, 
     file = "Robj/final_constM_mcmc.RData")












