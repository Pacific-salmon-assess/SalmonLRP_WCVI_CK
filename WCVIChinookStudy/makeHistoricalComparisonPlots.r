# This file runs function calls that compare LRP estimates for IF Coho that have been
# calculated from a range of different estimation methods.
# e.g., logistic-regression LRPs vs. projected LRPs vs. proportion-based LRPs

# Prior to running this code, the following files must be sourced to produce output csvs:
#   - runFraserCoho.r (to estimate logistic LRPs)
# (will eventually add runFraserCoho_projLRP.r, but have not added this in yet, so not yet plotting projected LRPs)

#library(MASS) # dose.p function to get SE around P95
library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
# library(tmbstan)
library(here)
library(zoo)

setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
wcviCKDir<-paste(rootDir,"/WCVIChinookStudy",sep="")

setwd(codeDir)

sourceAll <- function(){
  source("plotFunctions.r")
  source("helperFunctions.r")
}

sourceAll()

# ======================================================================
# Read-in Chinook data:
# =====================================================================
setwd(wcviCKDir)

EscpDat <- read.csv("DataIn/Inlet_Sum.csv")
# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
colnames(EscpDat)[colnames(EscpDat)=="Inlet_ID"] <- "CU"
colnames(EscpDat)[colnames(EscpDat)=="MU_Name"] <- "WCVI"
colnames(EscpDat)[colnames(EscpDat)=="BroodYear"] <- "yr"
colnames(EscpDat)[colnames(EscpDat)=="Spawners"] <- "Escp"

SRDat <- read.csv("DataIn/WCVI_SRbyCU.csv")



# Restrict data set to years 1990, year where all inlets start to have data
EscpDat <- EscpDat %>% filter(yr >= 1990)
SRDat <- SRDat %>% filter(BroodYear >= 1990)
genYrs <- 4

# Roll up escpaments, and get Gen Mean of htat
AggEscp <- EscpDat %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Escp, genYrs, gm_mean, fill = NA, align="right"))

# ===================================================================================
# Plot annual status with bars to show years in which LRP was breached
# ===================================================================================




# Specify which thresholds should be tested in the data-based proportional thresholds
ps_Prop<-c(1.0)

# Provide WSP assessment results
WSP_estYr<-2006
WSP_AboveLRP<-FALSE

LRP_estYr<-2020
retroYears<-1990:LRP_estYr
plotStatusBarsChinook_byYear(LRP_estYr, retroYears,  genYrs, AggEscp,
                             EscpDat=EscpDat, pLRP=0.5, ps_Prop=ps_Prop,
                             WSP_estYr=2016, WSP_AboveLRP=WSP_AboveLRP,
                             outDir = wcviCKDir,
                             fName = paste("statusPlot-withBars",
                                           LRP_estYr,sep=""))


