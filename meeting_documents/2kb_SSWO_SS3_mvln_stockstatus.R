#South Atlantic SWO 2026
#R-4.3.1

#SS3 results
####Kobe plot########

library(r4ss)
library(ss3diags)
library(reshape)
library(reshape2)
library(kobe) 
library(grid)
library(plyr)
library(ggplot2)
library(png)
library(dplyr)

#ss3
dir0 <- "D:/0_ICCAT_Assessment_Files/4_ALB_ATL/2026/SALB/"
#dir0 <- "F:/0_ICCAT_Assessment_Files/4_ALB_ATL/2026/SALB/"
dirRec <- paste0(dir0,"Recommendation/")

StartYr <- 1956 #start year of stock assessment 
EndYr <- 2024 #end year of stock assessment 
PrjYr <- 2040 #end year of projection 
outprob =c(0.5,0.025,0.975) #CI

limB10 <- 0.1 #to provide probability: stock (Bratio) below limB
limB20 <- 0.2 #to provide probability: stock (Bratio) below limB

#only Scenario 1 be used for the stock status
###read, mvln, summary##############
dirs <- paste0(dir0,"Recommendation/Stock_Synthesis/SALB_newF/")
#run <- list.files(dirs)
nrun <- c("S05_50th-VBF_50th-M")

dirns <- paste0(dirs,nrun)
#rename
nrun2 <- c("S05_50th-VBF_50th-M")


#1 read SS3 results #############

ss3detall <- NULL
mvln <- NULL

i <- 1
#ss3 reference case
  ss3rep1 <- SS_output(dir = dirns[i],verbose = TRUE, printstats = TRUE,covar=TRUE)
  
  #output deterministic, already adjusted SSB at the end of the year.
  ###!!!! ADJUSTED 
  #result scenario1
  
  ss3rep <- ss3rep1
  ss3det <- data.frame(year=c(StartYr:EndYr))
  #ss3rep$derived_quants$Label[c(285:353)+1] #Bratio
  ss3det$stock <- ss3rep$derived_quants$Value[c(285:353)+1] #Label, Value,StdDev, SSB/ssbmsy
  ss3det$stocksd <- ss3rep$derived_quants$StdDev[c(285:353)+1] #Value,StdDev
  #ss3rep$derived_quants$Label[c(215:283)] #F
  ss3det$harvest <- ss3rep$derived_quants$Value[c(215:283)] #Label, Value,StdDev, F/Fmsy
  ss3det$harvestsd <- ss3rep$derived_quants$StdDev[c(215:283)] #Value,StdDev
  #ss3rep$derived_quants$Label[c(370)] # MSY 
  ss3det$msy <- ss3rep$derived_quants$Value[c(370)] # 24333 (23190- 25476)
  ss3det$msysd <- ss3rep$derived_quants$StdDev[c(370)] 
  
  ss3det$stockLCI <- ss3det$stock-1.96*ss3det$stocksd
  ss3det$stockUCI <- ss3det$stock+1.96*ss3det$stocksd
  ss3det$harvestLCI <- ss3det$harvest-1.96*ss3det$harvestsd
  ss3det$harvestUCI <- ss3det$harvest+1.96*ss3det$harvestsd
  ss3det$msyLCI <- ss3det$msy-1.96*ss3det$msysd
  ss3det$msyUCI <- ss3det$msy+1.96*ss3det$msysd
  
    
  ss3det$model <- paste0("SS3_",nrun2[i])
  
  ss3detall <- rbind(ss3detall,ss3det)

  #mvln
  sspar(mfrow=c(1,1),plot.cex = 0.9)
  mvln1 = SSdeltaMVLN(ss3rep1, plot = T, mc = 10000, Fref="MSY",
                      years = c(1956:2025)) 
  mvnlkb <- mvln1$kb
  mvnlkb$run <- paste0("SS3_",nrun2[i])
  mvln <- rbind(mvln,mvnlkb)


#save all output kbs, and assessment report
write.csv(ss3det,file=paste0(dirRec,"kobe/SALB_all_SS3_deterministic_adj.csv"),row.names=F)
save(mvln,file=paste0(dirRec,"kobe/SALB_all_MVLN.Rdata"))

#ssb year adjustment;'SSB at the end of the year

Skbs <- mvln
tm <- subset(Skbs,,c(year,run,iter,stock))
tm$year <- tm$year-1 #ssb adjustment

Skbs <- Skbs[,-5] #remove "stock"
Skbs <- merge(Skbs,tm,all=T)

Skbs <- Skbs[Skbs$year>=1957 & Skbs$year<= 2024, ]

#SS3 MVLN all runs# Skbs
Skbs$bcrit10 <- Skbs$stock/limB10
Skbs$bcrit20 <- Skbs$stock/limB20

###save all iterations
save(Skbs,file=paste0(dirRec,"kobe/SALB_MVLN_adj_1957_2024.Rdata"))

