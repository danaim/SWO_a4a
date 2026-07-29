rm(list = ls())
library(FLa4a)
library(ggplotFL)
library(coda)
load("Robj/input4mcmc.Rdata")

#====================================================================
# Setup - use same model as candidate run
#====================================================================
stk   <- stk85_c
idx   <- idx_use
fmod  <- ~ s(year, k = 20) + s(age, k = 8) + ti(year, age, k = c(6, 5))
srmod <- ~ s(year, k = 20)
qmod  <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3), ~1)
yr    <- ac(range(stk)["maxyear"])

#====================================================================
# 1. MLE fit
#====================================================================
fit_mle <- sca(stk, idx, fmodel = fmod, srmodel = srmod, qmodel = qmod)
stk_mle <- stk + fit_mle
stk_mle_sim <- stk + simulate(fit_mle, nsim = 1000)

#====================================================================
# 2. MCMC fit (single chain, tuned mcprobe)
#====================================================================
mcprobe <- 0.4   # use value from mcprobe tuning
it       <- 2000
brn      <- 500
mcsave   <- 1000

fit_mcmc <- sca(stk, idx, fmodel = fmod, srmodel = srmod, qmodel = qmod,
                fit = "MCMC",
                mcmc = SCAMCMC(mcmc = mcsave * it, mcsave = mcsave, mcprobe = mcprobe))
fit_b    <- burnin(fit_mcmc, brn)
stk_mcmc <- stk + fit_b

#====================================================================
# 3. Convergence check
#====================================================================
pars_mat <- t(coef(stkmodel(pars(fit_b)))@.Data)
chain    <- coda::mcmc(pars_mat)
cat("Acceptance rate:", fitSumm(fit_mcmc)["accrate", ], "\n")
cat("Min ESS:", min(effectiveSize(chain)), "\n")
plot(chain[, 1:min(4, ncol(chain))], main = "Parameter traces")

## Terminal year traces
par(mfrow = c(1, 2))
plot(c(ssb(stk_mcmc)[, yr]), type = "l", main = paste("SSB", yr), ylab = "SSB")
abline(h = median(c(ssb(stk_mcmc)[, yr])), col = "red", lty = 2)
plot(c(fbar(stk_mcmc)[, yr]), type = "l", main = paste("Fbar", yr), ylab = "F")
abline(h = median(c(fbar(stk_mcmc)[, yr])), col = "red", lty = 2)
par(mfrow = c(1, 1))

#====================================================================
# 4. MLE vs MCMC trajectory comparison
#====================================================================
plot(FLStocks(MLE = stk_mle_sim, MCMC = stk_mcmc))

#====================================================================
# 5. Terminal year distributions
#====================================================================
ssb_mle  <- c(ssb(stk_mle_sim)[, yr])
ssb_mcmc <- c(ssb(stk_mcmc)[, yr])
f_mle    <- c(fbar(stk_mle_sim)[, yr])
f_mcmc   <- c(fbar(stk_mcmc)[, yr])

par(mfrow = c(1, 2))
hist(ssb_mcmc, breaks = 30, freq = FALSE, col = rgb(0, 0, 1, 0.3),
     main = paste("SSB", yr), xlab = "SSB")
hist(ssb_mle, breaks = 30, freq = FALSE, col = rgb(1, 0, 0, 0.3), add = TRUE)
abline(v = c(ssb(stk_mle)[, yr]), col = "red", lwd = 2)
legend("topright", c("MCMC", "MLE sim"), fill = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.3)), bty = "n")

hist(f_mcmc, breaks = 30, freq = FALSE, col = rgb(0, 0, 1, 0.3),
     main = paste("Fbar", yr), xlab = "F")
hist(f_mle, breaks = 30, freq = FALSE, col = rgb(1, 0, 0, 0.3), add = TRUE)
abline(v = c(fbar(stk_mle)[, yr]), col = "red", lwd = 2)
legend("topright", c("MCMC", "MLE sim"), fill = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.3)), bty = "n")
par(mfrow = c(1, 1))

#====================================================================
# 6. Summary table
#====================================================================
summary_df <- data.frame(
  Metric = c("SSB", "Fbar"),
  MLE = round(c(c(ssb(stk_mle)[, yr]), c(fbar(stk_mle)[, yr])), 3),
  MLE_CI_lo = round(c(quantile(ssb_mle, 0.1), quantile(f_mle, 0.1)), 3),
  MLE_CI_hi = round(c(quantile(ssb_mle, 0.9), quantile(f_mle, 0.9)), 3),
  MCMC_med = round(c(median(ssb_mcmc), median(f_mcmc)), 3),
  MCMC_CI_lo = round(c(quantile(ssb_mcmc, 0.1), quantile(f_mcmc, 0.1)), 3),
  MCMC_CI_hi = round(c(quantile(ssb_mcmc, 0.9), quantile(f_mcmc, 0.9)), 3)
)
print(summary_df)

    ## Compare terminal year SSB distributions                                                        
boxplot(list(                                                                                     
      MLE  = c(ssb(stk_mle_sim)[, yr]),                                                               
      MCMC = c(ssb(stk_mcmc)[, yr])                                                               
    ), main = "SSB Terminal Year"                                                                     
) 


save(fit_mle, fit_mcmc, stk_mle, stk_mcmc, stk_mle_sim, summary_df,
     file = "Robj/mle_vs_mcmc_comparison.RData")
