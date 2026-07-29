library(FLCore)
library(FLa4a)
library(FLBRP)
library(ggplot2);theme_set(theme_bw())
rm(list = ls())
source("analysis/plot_functions.R")
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

### ===================================================================== ###
### Final run ----------------------------------------------------------- ###
### ===================================================================== ###

## Proposal of the group: Exclude Sicilian LL - use only four indices
## Run 1: 1985 Constant M
stk85_c <- replaceZeros(window(swo.stk, start = 1985))
range(stk85_c)[c('minfbar', 'maxfbar')] <- c(1, 4)
fmod_85 <- ~ s(year, k = 20) + s(age, k = 8) + ti(year, age, k = c(6, 5))
srmod   <- ~ s(year, k = 20)
qmod    <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3), ~1)

fit85_c <- sca(stk85_c, idx_use, fmodel = fmod_85, srmodel = srmod, qmodel = qmod)
a4a85_c <- stk85_c + simulate(fit85_c, nsim = 1000)
plot(a4a85_c)

res85_c <- residuals(fit85_c, stk85_c, idx_gp)

res_df <- as.data.frame(residuals(fit85_c, stk85_c, idx_use))
  
ggplot(subset(res_df, qname == "catch.n"), aes(x = year, y = data)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_segment(aes(xend = year, yend = 0, color = data > 0), linewidth = 0.5) +
      geom_point(aes(color = data > 0), size = 1.5) +
      scale_color_manual(values = c("TRUE" = "#2c7bb6", "FALSE" = "#d7191c"), 
                         labels = c("TRUE" = "Positive", "FALSE" = "Negative")) +
      facet_wrap(~ age, ncol = 3) +
      labs(x = "Year", y = "Standardized Residuals", color = "Residual Sign") +
      theme_bw() +
      theme(legend.position = "bottom") + ylim(-3,3)

ggplot(subset(res_df, !qname %in% c("catch.n", "catch")), aes(x = year, y = data)) +              
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_segment(aes(xend = year, yend = 0, color = data > 0), linewidth = 0.5) +
      geom_point(aes(color = data > 0), size = 1.5) +
      scale_color_manual(values = c("TRUE" = "#2c7bb6", "FALSE" = "#d7191c"), 
                         labels = c("TRUE" = "Positive", "FALSE" = "Negative")) +
      facet_wrap(~ qname, ncol = 3) +
      labs(x = "Year", y = "Standardized Residuals", color = "Residual Sign") +
      theme_bw() +
      theme(legend.position = "bottom")

cd <- computeCatchDiagnostics(fit85_c, stk85_c)
obs_df <- as.data.frame(cd[["obs"]]); obs_df <- obs_df[, c("year", "data")]
est_df <- as.data.frame(cd[["est"]]); est_df <- est_df[, c("year", "data")]
    ## 1. Top-Left: Observation Error (oe)
    p1 <- plot_catch_err(cd[["oe"]], "Observation Error (oe)")
  
    ## 2. Top-Right: Estimation Error (ee)
    p2 <- plot_catch_err(cd[["ee"]], "Estimation Error (ee)")
  
    ## 3. Bottom-Left: Combined Error (oee)
    p3 <- plot_catch_err(cd[["oee"]], "Prediction Error (oee)")
  
    ## 4. Bottom-Right: Standardized Catch Residuals (resstd)
    res_df <- as.data.frame(cd[["resstd"]])
    p4 <- ggplot(res_df, aes(x = year, y = data)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_segment(aes(xend = year, yend = 0, color = data > 0), linewidth = 0.6) +
      geom_point(aes(color = data > 0), size = 1.8) +
      geom_smooth(se = FALSE, color = "black", linewidth = 0.6) +
      scale_color_manual(values = c("TRUE" = "#2b83ba", "FALSE" = "#d7191c")) +
      labs(title = "Standardized Catch Residuals (resstd)", x = "Year", y = "Residuals") +            
      theme_bw(base_size = 9) +
      theme(legend.position = "none")
  
    ## Combine into a 2x2 layout
    (p1 | p2) / (p3 | p4)

retro_runs <- function(stk, idxs, retro = 5, kfrac = 0.3, k_f = 18, k_sr = 20, qmod) {
  lst0 <- split(0:retro, 0:retro)
  lst0 <- lapply(lst0, function(x) {
    yr <- range(stk)["maxyear"] - x
    stk_sub <- window(stk, end = yr)
    idx_sub <- FLIndices(window(idxs, end = yr))
    KY <- unname(k_f - floor(x * kfrac))
    KZ <- unname(k_sr - floor(x * kfrac))
    fmod_sub <- as.formula(substitute(~ s(year, k = KY) + s(age, k = 8) + ti(year, age, k = c(6, 5)), list(KY = KY)))
    srmod_sub <- as.formula(substitute(~ s(year, k = KZ), list(KZ = KZ)))
    fit_sub <- sca(stk_sub, idx_sub, fmodel = fmod_sub, srmodel = srmod_sub, qmodel = qmod)
    stk_sub + fit_sub
  })
  FLStocks(lst0)
}

retro85_c <- retro_runs(stk85_c, idx_use, retro = 5, k_f = 18, k_sr = 20, qmod = qmod)
mohn_ssb <- icesAdvice::mohn(mohnMatrix(retro85_c, ssb))
mohn_f   <- icesAdvice::mohn(mohnMatrix(retro85_c, fbar))
plot(retro_fn) + theme(legend.position = "none")+
      labs(
        subtitle = sprintf("Mohn's rho: SSB = %.3f  |  Fbar = %.3f", mohn_ssb, mohn_f)
      )
plot_status(stk85_c, fit85_c,a4a85_c)

df <- as.data.frame(idx_use)
df <- df[slot == 'index',]
p1 <- ggplot(df, 
             aes(x = year, 
                 y = data, group = cname, color = cname)) +
  geom_line(alpha = 0.85) +
  facet_wrap(~cname)+
  scale_color_viridis_d(option = "turbo") +
  labs(title = "CPUE indices", x = "Age", y = "Relative biomass") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = c(0:9))
plot(idx_use) + geom_line(linesize = 1.2)
