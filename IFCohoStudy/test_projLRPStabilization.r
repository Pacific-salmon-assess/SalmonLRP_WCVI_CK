


library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(tmbstan)
library(here)
library(zoo)

setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
cohoDir<-paste(rootDir,"/IFCohoStudy",sep="")

setwd(codeDir)

sourceAll <- function(){
  source("ProjLRP_Functions.r")
  source("plotFunctions.r")
  source("helperFunctions.r")
}
sourceAll()



# ======================================================================
# Read-in Coho data:  
# =====================================================================
setwd(cohoDir)

CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"

CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")


# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)
CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)

# Roll up escapements, and get Gen Mean of that
genYrs <- 3
AggEscp <- CoEscpDat %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Escp, genYrs, gm_mean, fill = NA, align="right"))


# ======================================================================
# Specify data sets used for projections  
# =====================================================================

# Additional projection parameters
BroodYrLag <- 2
pList <- seq(0.2,1,by=0.2)
year <- 2018


# Only use SR data for brood years that have recruited by specified year
# (note: most recent brood year is calculated by subtracting BroodYearLag (e.g. 2 years) from current year)
SRDat <- CoSRDat %>%  filter(BroodYear <= year-BroodYrLag)
EscpDat.yy <- CoEscpDat %>% filter(yr <= year)

SRDat$yr_num <- group_by(SRDat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
SRDat$CU_ID <- group_by(SRDat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
EscDat <- EscpDat.yy %>%  right_join(unique(SRDat[,c("CU_ID", "CU_Name")]))



# ===================================================================
# Estimate LRPs
# ==================================================================

# Read in projection outputs to create input lists for logistic regression


OMsToInclude<-c("IM.base.nyrs10000_nsims10", "IM.base.nyrs7000_nsims100", "IM.base.nyrs30_nsims5000",
                "IM.base.nyrs100_nsims2000")

for (i in 1:length(OMsToInclude)) {
  filename<-paste("projLRPDat_",OMsToInclude[i],".csv",sep="")
  dat.i<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
  dat.i<-dat.i %>% filter(year > max(SRDat$yr_num)+4)
  dat.i$OM.Name<-OMsToInclude[i]
  if (i == 1) projLRPDat<-dat.i
  if (i > 1) projLRPDat<-rbind(projLRPDat,dat.i)
  
  filename<-paste( "projSpwnDat_",OMsToInclude[i],".csv",sep="")
  spDat.i<-read.csv(here(cohoDir,"SamSimOutputs", "simData",filename))
  spDat.i$OM.Name<-OMsToInclude[i]
  if (i == 1) projCUSpDat<-spDat.i
  if (i > 1) projCUSpDat<-rbind(projCUSpDat,spDat.i)
}



# Calculate LRP and associated probability interval based on distribution of sAg 
# -- Calculate and save LRPs by OM 
LRPs_byOM<-projLRPDat %>% group_by(OM.Name,ppnCUsLowerBM) %>% 
  summarise(LRP.50=median(sAg), LRP.95=quantile(sAg,0.95),LRP.05=quantile(sAg,0.05))


LRPs<-LRPs_byOM



# ===================================================================
# Estimate Alternate LRPs
# ==================================================================

# Read in projection outputs to create input lists for logistic regression

OMsToInclude<-c("IM.base.nyrs30_nsims5000")
probThresh<-0.99

for (i in 1:length(OMsToInclude)) {
  filename<-paste("projLRPDat_",OMsToInclude[i],".csv",sep="")
  dat.i<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
  dat.i<-dat.i %>% filter(year > max(SRDat$yr_num)+4)
  dat.i$OM.Name<-OMsToInclude[i]
  if (i == 1) projLRPDat<-dat.i
  if (i > 1) projLRPDat<-rbind(projLRPDat,dat.i)
  
  # filename<-paste( "projSpwnDat_",OMsToInclude[i],".csv",sep="")
  # spDat.i<-read.csv(here(cohoDir,"SamSimOutputs", "simData",filename))
  # spDat.i$OM.Name<-OMsToInclude[i]
  # if (i == 1) projCUSpDat<-spDat.i
  # if (i > 1) projCUSpDat<-rbind(projCUSpDat,spDat.i)
}


# All replicates combined:


projLRPDat<-projLRPDat

minBreak<-0#round(min(projLRPDat$sAg),digits=-2)
maxBreak<-round(max(projLRPDat$sAg),digits=-2)

breaks<-seq(minBreak, maxBreak,by=200)

projLRPDat$bins<-cut(projLRPDat$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)))

tmp<-projLRPDat %>% group_by(bins) %>% summarise(nSims=(length(ppnCUsLowerBM)))

# Note: temporary step to filter out bins with < 100 nSims in below call
tmp2<-projLRPDat %>% group_by(bins) %>% summarise(nSimsProp1=(length(ppnCUsLowerBM[ppnCUsLowerBM == 1.0]))) %>%
  add_column(nSims=tmp$nSims) %>% filter(nSims>=100)




projLRPDat<-tmp2 %>% add_column(prob=tmp2$nSimsProp1/tmp2$nSims)
projLRPDat$diff<-abs(probThresh-projLRPDat$prob)

LRP<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))


plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19, xlim=c(0,85000), cex=0.2,
     xlab="Aggregate Abundance", ylab="Pr (All CUs > Lower Benchmark)")
abline(h=probThresh, lty=2)
#abline(v=LRP, col="orange", lwd=2)


barplot(height = projLRPDat$nSims,names.arg = projLRPDat$bins)




# 
# 
# 
# # If looking at individual replicates:
# 
# nProj<-max(projLRPDat$iteration)
# 
# LRP<-rep(NA,nProj)
# 
# for (i in 1:nProj) {
#   
#   projLRPDat.i<-projLRPDat %>% filter(iteration==i)
#   
#   minBreak<-round(min(projLRPDat.i$sAg),digits=-2)
#   maxBreak<-round(max(projLRPDat.i$sAg),digits=-2)
#   
#   breaks<-seq(minBreak, maxBreak,by=250)
#   
#   projLRPDat.i$bins<-cut(projLRPDat.i$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)))
#   
#   tmp<-projLRPDat.i %>% group_by(bins) %>% summarise(nSims=(length(ppnCUsLowerBM)))
#   
#   # Note: temporary step to filter out bins with < 100 nSims in below call
#   tmp2<-projLRPDat.i %>% group_by(bins) %>% summarise(nSimsProp1=(length(ppnCUsLowerBM[ppnCUsLowerBM == 1.0]))) %>%
#     add_column(nSims=tmp$nSims) %>% filter(nSims>=100)
#   
#   projLRPDat.i<-tmp2 %>% add_column(prob=tmp2$nSimsProp1/tmp2$nSims)
#   projLRPDat.i$diff<-abs(probThresh-projLRPDat.i$prob)
#   
#   LRP[i]<-as.numeric(as.character(projLRPDat.i$bins[projLRPDat.i$diff == min(projLRPDat.i$diff)]))
#   
#   LRP.i<-LRP[i]
#   
#   plot(projLRPDat.i$bins,projLRPDat.i$prob, pch=19, 
#        xlab="Aggregate Abundance", ylab="Pr (All CUs > Lower Benchmark)")
#   abline(h=probThresh, lty=2)
#   abline(v=LRP.i, col="orange", lwd=2)
#  
#  
# }
# 
# hist(LRP)
# abline(v=median(LRP), col="red")
# abline(v=quantile(LRP,0.95)[[1]], col="red", lty=2) 
# abline(v=quantile(LRP,0.05)[[1]], col="red", lty=2) 
# 
# print(quantile(LRP,c(0.05,0.5,0.95)))
