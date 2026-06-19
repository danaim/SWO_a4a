library(FLCore)
library(FLa4a)
library(FLBRP)
library(ggplot2);theme_set(theme_bw())
rm(list = ls())

indices <- read.csv("data/indices_all.csv")
indices

idx_bio <- FLIndices()

# GRC
dms <- list(age="all", year=1987:2024)
idx_bio[['GRC_LL']] <- FLIndexBiomass(index = FLQuant(indices$GRC_LL_std[indices$Year %in% 1987:2024], 
  dimnames=dms))

# LIG early
dms <- list(age="all", year=1991:2009)
idx_bio[['LIG_LL_early']] <- FLIndexBiomass(index = FLQuant(indices$w_LIG_LL_early_std[indices$Year %in% 1991:2009], 
  dimnames=dms))

# LIG
dms <- list(age="all", year=2009:2022)
idx_bio[['LIG_LL']] <- FLIndexBiomass(index = FLQuant(indices$w_LIG_LL_std[indices$Year %in% 2009:2022], 
  dimnames=dms))

# SPN
dms <- list(age="all", year=1988:2023)
idx_bio[['SPN_LL']] <- FLIndexBiomass(index = FLQuant(indices$w_SPN_LL_std[indices$Year %in% 1988:2023], 
  dimnames=dms))

# MOR
dms <- list(age="all", year=2012:2024)
idx_bio[['MOR_LL']] <- FLIndexBiomass(index = FLQuant(indices$MOR_LL_std[indices$Year %in% 2012:2024], 
  dimnames=dms))

for(i in 1:length(idx_bio)) {
  range(idx_bio[[i]],c("startf","endf")) = c(0,1)
  range(idx_bio[[i]],c("min","max")) = c(1,3)
  units(index(idx_bio[[i]])) <- "t"
}

saveRDS(idx_bio, file = 'Robj/swo_bio_idx.rds')

plot(idx_bio)
index(idx_bio[[1]])
