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
for(i in 1:length(swo.idx)) {
  range(swo.idx[[i]],c("min","max")) = c(1,4)
  units(index(swo.idx[[i]])) <- "t"
}


### ------------------------------------------------------ ###
# Candidate runs:
### ------------------------------------------------------ ###

### ====================================================== ###
# 1. From 1985 - 5 indices
# Constant M
stk <- window(swo.stk, start = 1985)
# stk <- window(swo.stk, start = 1978)
stk <- replaceZeros(stk)
range(stk)[c('minfbar','maxfbar')] <- c(1,4) 
idx <- swo.idx[-2]

srmod <- ~s(year, k =20)
fmod3 <- ~s(year, k = 17) + s(age, k = 8) + ti(year, age, k = c(6,5))
# fmod3 <- ~s(year, k = 20) + factor(age) + ti(year, age, k = c(6,5))
qmod <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3),~1)
fit <- sca(stk, idx, fmodel = fmod3, qmodel = qmod, srmodel = srmod) 
# fit <- sca(stk, idx, fmodel = fmod3, qmodel = qmod) 

a4a.stk <- stk + simulate(fit, nsim = 1000)
res <- residuals(fit, stk, idx)
plot(res)
bubbles(res)
plot(a4a.stk, stk)

plot_status(stk, fit, a4a.stk)

fmod1 <- ~s(year, k = 20) + s(age, k = 8) + ti(year, age, k = c(6,5))
fmod2 <- ~s(year, k = 21) + s(age, k = 8) + ti(year, age, k = c(6,5))
fmod3 <- ~s(year, k = 22) + s(age, k = 8) + ti(year, age, k = c(6,5))
fmod4 <- ~s(year, k = 23) + s(age, k = 8) + ti(year, age, k = c(6,5))
fmod5 <- ~s(year, k = 24) + s(age, k = 8) + ti(year, age, k = c(6,5))
fmod6 <- ~s(year, k = 19) + s(age, k = 8) + ti(year, age, k = c(6,5))
fmod7 <- ~s(year, k = 18) + s(age, k = 8) + ti(year, age, k = c(6,5))
fmod8 <- ~s(year, k = 17) + s(age, k = 8) + ti(year, age, k = c(6,5))
fmod9 <- ~s(year, k = 16) + s(age, k = 8) + ti(year, age, k = c(6,5))
fmod10 <- ~s(year, k = 15) + s(age, k = 8) + ti(year, age, k = c(6,5))

fmods <- list(fmod1, fmod2, fmod3, fmod4, fmod5,fmod6,
fmod7,fmod8,fmod9,fmod10)


fits <-scas(FLStocks(stk), list(idx),fmodel = fmods,srmodel =  list(srmod),qmodel =  list(qmod))
plot(fits)

retro <- function(stk, idxs, retro=5, kfrac="missing", k, k2, ...){
  args <- list(...)
  if(missing(kfrac)) kfrac <- 0.3
  lst0 <- split(0:retro, 0:retro)
  lst0 <- lapply(lst0, function(x){
    yr <- range(stk)["maxyear"] - x
    args$stock <- window(stk, end=yr)
    args$indices <- FLIndices(window(idxs, end=yr))
    KY <- unname(k - floor(x*kfrac))
    KZ <- unname(k2 - floor(x*kfrac))
    # NOTE THIS HAS TO BE ADAPTED FOR THE MODEL USED  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # fmod <- substitute(~s(year, k = 20) + s(age, k = 6) + ti(year, age, k = c(6,5)))
    # fmod <- substitute(~s(year, k = KY) + factor(age) + ti(year, age, k = c(6,5), list(KY = KY)))
    srmod <- substitute(~s(year, k = KZ), list(KZ=KZ))
    fmod <- substitute(~s(year, k = KY) + s(age, k = 8) + ti(year, age, k = c(6,5)),list(KY=KY))
    args$fmodel <- as.formula(fmod)
    # args$srmodel <- as.formula(srmod)
    fit <- do.call("sca", args)
    args$stock + fit
  })
  FLStocks(lst0)
}

retro.lst <- retro(stk, idx, k=20, k2 = 20, qmodel = qmod)
p<-plot(retro.lst)+theme(legend.position="none")
p


### ====================================================== ###
# 1. From 1985 - 5 indices
# Lorenzen
stk.l <- window(swo.stk.l, start = 1985)
# stk <- window(swo.stk, start = 1978)
stk.l <- replaceZeros(stk.l)
range(stk.l)[c('minfbar','maxfbar')] <- c(1,4) 
idx <- swo.idx[-2]

srmod <- ~s(year, k =20)
fmod3 <- ~s(year, k = 20) + s(age, k = 8) + ti(year, age, k = c(6,5))
# fmod3 <- ~s(year, k = 20) + factor(age) + ti(year, age, k = c(6,5))
qmod <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3),~1)
fit.l <- sca(stk.l, idx, fmodel = fmod3, qmodel = qmod, srmodel = srmod) 
# fit <- sca(stk, idx, fmodel = fmod3, qmodel = qmod) 

a4a.stk.l <- stk.l + simulate(fit.l, nsim = 1000)
res <- residuals(fit.l, stk.l, idx)
plot(res)
bubbles(res)
plot(a4a.stk.l, stk.l)

plot_status(stk.l, fit.l)

plot(FLStocks(loren = a4a.stk.l, const = a4a.stk))

retro <- function(stk, idxs, retro=5, kfrac="missing", k, k2, ...){
  args <- list(...)
  if(missing(kfrac)) kfrac <- 0.5
  lst0 <- split(0:retro, 0:retro)
  lst0 <- lapply(lst0, function(x){
    yr <- range(stk)["maxyear"] - x
    args$stock <- window(stk, end=yr)
    args$indices <- FLIndices(window(idxs, end=yr))
    KY <- unname(k - floor(x*kfrac))
    KZ <- unname(k2 - floor(x*kfrac))
    # NOTE THIS HAS TO BE ADAPTED FOR THE MODEL USED  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # fmod <- substitute(~s(year, k = 20) + s(age, k = 6) + ti(year, age, k = c(6,5)))
    # fmod <- substitute(~s(year, k = KY) + factor(age) + ti(year, age, k = c(6,5), list(KY = KY)))
    srmod <- substitute(~s(year, k = KZ), list(KZ=KZ))
    fmod <- substitute(~s(year, k = KY) + factor(age) + ti(year, age, k = c(6,5)),list(KY=KY))
    args$fmodel <- as.formula(fmod)
    # args$srmodel <- as.formula(srmod)
    fit <- do.call("sca", args)
    args$stock + fit
  })
  FLStocks(lst0)
}

retro.lst <- retro(stk.l, idx, k=25, k2 = 20, qmodel = qmod)
p<-plot(retro.lst)+theme(legend.position="none")
p


plot(FLStocks(constntM = a4a.stk, lorenzen = a4a.stk.l))

### ====================================================== ###
# 1. From 1978 - 5 indices
# Constant M
# stk <- window(swo.stk, start = 1985)
stk <- window(swo.stk, start = 1978)
stk <- replaceZeros(stk)
range(stk)[c('minfbar','maxfbar')] <- c(1,4) 
idx <- swo.idx[-2]


# srmod2 <- ~I(as.numeric(year<=2024)) + 1 + s(year, k = 15,by = as.numeric(year>=1987))
srmod <- ~s(replace(year,year<=1985,1985), k = 20)
srmod <- ~s(year, k =20)+I(as.numeric(year<=2024)) + 1 
fmod3 <- ~s(year, k = 25) + s(age, k = 8) + ti(year, age, k = c(6,5))
qmod <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3),~1)
fit <- sca(stk, idx, fmodel = fmod3, qmodel = qmod, srmodel = srmod) 
# fit <- sca(stk, idx, fmodel = fmod3, qmodel = qmod) 

a4a.stk <- stk + simulate(fit, nsim = 1000)
res <- residuals(fit, stk, idx)
plot(res)
bubbles(res)
plot(a4a.stk, stk)

plot_status(stk, fit)


retro <- function(stk, idxs, retro=5, kfrac="missing", k, k2, ...){
  args <- list(...)
  if(missing(kfrac)) kfrac <- 0.5
  lst0 <- split(0:retro, 0:retro)
  lst0 <- lapply(lst0, function(x){
    yr <- range(stk)["maxyear"] - x
    args$stock <- window(stk, end=yr)
    args$indices <- FLIndices(window(idxs, end=yr))
    KY <- unname(k - floor(x*kfrac))
    KZ <- unname(k2 - floor(x*kfrac))
    # NOTE THIS HAS TO BE ADAPTED FOR THE MODEL USED  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # fmod <- substitute(~s(year, k = 20) + s(age, k = 6) + ti(year, age, k = c(6,5)))
    # fmod <- substitute(~s(year, k = KY) + factor(age) + ti(year, age, k = c(6,5), list(KY = KY)))
    srmod <- substitute(~s(year, k = KZ), list(KZ=KZ))
    fmod <- substitute(~s(year, k = KY) + factor(age) + ti(year, age, k = c(6,5)),list(KY=KY))
    args$fmodel <- as.formula(fmod)
    # args$srmodel <- as.formula(srmod)
    fit <- do.call("sca", args)
    args$stock + fit
  })
  FLStocks(lst0)
}

retro.lst <- retro(stk, idx, k=25, k2 = 20, qmodel = qmod)
p<-plot(retro.lst)+theme(legend.position="none")
p


### ====================================================== ###
# 1. From 1978 - 5 indices
# Lorenzen
# stk.l <- window(swo.stk.l, start = 1985)
stk.l <- window(swo.stk, start = 1978)
stk.l <- replaceZeros(stk.l)
range(stk.l)[c('minfbar','maxfbar')] <- c(1,4) 
idx <- swo.idx[-2]

srmod <- ~s(year, k =20)
fmod3 <- ~s(year, k = 25) + s(age, k = 8) + ti(year, age, k = c(6,5))
# fmod3 <- ~s(year, k = 25) + factor(age) + ti(year, age, k = c(6,5))
qmod <- list(~1, ~s(year, k = 3), ~1, ~s(year, k = 3),~1)
fit.l <- sca(stk.l, idx, fmodel = fmod3, qmodel = qmod, srmodel = srmod) 
# fit <- sca(stk, idx, fmodel = fmod3, qmodel = qmod) 

a4a.stk.l <- stk.l + simulate(fit.l, nsim = 1000)
res <- residuals(fit.l, stk.l, idx)
plot(res)
bubbles(res)
plot(a4a.stk.l, stk.l)

plot_status(stk.l, fit.l)


retro <- function(stk, idxs, retro=5, kfrac="missing", k, k2, ...){
  args <- list(...)
  if(missing(kfrac)) kfrac <- 0.5
  lst0 <- split(0:retro, 0:retro)
  lst0 <- lapply(lst0, function(x){
    yr <- range(stk)["maxyear"] - x
    args$stock <- window(stk, end=yr)
    args$indices <- FLIndices(window(idxs, end=yr))
    KY <- unname(k - floor(x*kfrac))
    KZ <- unname(k2 - floor(x*kfrac))
    # NOTE THIS HAS TO BE ADAPTED FOR THE MODEL USED  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # fmod <- substitute(~s(year, k = 20) + s(age, k = 6) + ti(year, age, k = c(6,5)))
    # fmod <- substitute(~s(year, k = KY) + factor(age) + ti(year, age, k = c(6,5), list(KY = KY)))
    srmod <- substitute(~s(year, k = KZ), list(KZ=KZ))
    fmod <- substitute(~s(year, k = KY) + factor(age) + ti(year, age, k = c(6,5)),list(KY=KY))
    args$fmodel <- as.formula(fmod)
    # args$srmodel <- as.formula(srmod)
    fit <- do.call("sca", args)
    args$stock + fit
  })
  FLStocks(lst0)
}

retro.lst <- retro(stk.l, idx, k=25, k2 = 20, qmodel = qmod)
p<-plot(retro.lst)+theme(legend.position="none")
p


ls()

cd <- computeCatchDiagnostics(fit85_c, stk85_c)
    library(ggplot2)                                                                                  
    library(gridExtra) # or library(gridExtra)                                                        
                                                                                                      
    ## Assuming 'cd' is the output from computeCatchDiagnostics(fit, stk)                             
                                                                                                      
    ## Helper function to summarize simulated FLQuant (median + 90% CI)                               
    sum_sim <- function(flq) {
      df <- as.data.frame(flq)
      agg <- aggregate(data ~ year, data = df, FUN = function(x) c(
        med = median(x, na.rm = TRUE),
        lwr = unname(quantile(x, 0.05, na.rm = TRUE)),
        upr = unname(quantile(x, 0.95, na.rm = TRUE))
      ))
      mat <- as.data.frame(agg$data)
      names(mat) <- c("med", "lwr", "upr")
      data.frame(year = agg$year, mat)
    }

    ## Extract baseline deterministic observed and estimated catch series
    obs_df <- as.data.frame(cd[["obs"]]); obs_df <- obs_df[, c("year", "data")]
    est_df <- as.data.frame(cd[["est"]]); est_df <- est_df[, c("year", "data")]
  
    ## Helper to build catch error comparison plots
    plot_catch_err <- function(flq_sim, title_text) {
      sim_data <- sum_sim(flq_sim)
      ggplot() +
        geom_ribbon(data = sim_data, aes(x = year, ymin = lwr, ymax = upr), 
                    fill = "#2b83ba", alpha = 0.25) +
        geom_line(data = sim_data, aes(x = year, y = med), 
                  color = "#2b83ba", linetype = "dashed", linewidth = 0.7) +
        geom_line(data = est_df, aes(x = year, y = data), 
                  color = "black", linewidth = 0.8) +
        geom_point(data = obs_df, aes(x = year, y = data), 
                   color = "#d7191c", size = 1.8) +
        labs(title = title_text, x = "Year", y = "Catch (t)") +
        theme_bw(base_size = 9)
    }
  
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


nyears <- 5
nsq    <- 3
ars.stk <- stk85_c + fit85_c
hc <- a4adiags::a4ahcxval(ars.stk, FLIndices(idx[[1]]), nyears = nyears, nsq = nsq, 
                srmodel = srmodel(fit85_c), 
                fmodel = fmodel(fit85_c), 
                qmodel = formula(qmodel(fit85_c)))

mohn_ssb <- icesAdvice::mohn(mohnMatrix(retro_85c, ssb))
mohn_f   <- icesAdvice::mohn(mohnMatrix(retro_85c, fbar))
mohn_ssb
