# Code by L. Warkentin & K. Holt
# Started on: October 16, 2020

# This file contains the code required to explore South Coast Chum escapement datasets and run retrospective analyses of 
#   data-based LRPs that use logistic regressions

library(rsample)
library(tidyverse)
#library(ggplot2)
#library(tidyr)
#library(dplyr)
options(scipen=1000000)
library(gridExtra)
library(reshape2)
library(TMB)
# Note: this package list has been copied from runFraserCoho.r; may not all be required for chum ... 

setwd('..') # This works if working directory is in ~/SalmonLRP_RetroEval/SCChumStudy folder
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
chumDir<-paste(rootDir,"/SCChumStudy",sep="")

setwd(codeDir)

sourceAll <- function(){
  source("benchmarkFunctions.r")
  source("LRPFunctions.r")
  source("plotFunctions.r")
  source("retroFunctions.r")
  source("helperFunctions.r")
}
sourceAll()

# Load TMB models
compile("TMB_Files/SR_IndivRicker_NoSurv.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_NoSurv"))

# Switch to chum directory
setwd(chumDir)

# Run to re-do infilling and brood table
# source("prepare_data.r") 

# ====================================================================
# Read in data and format for using in retrospective analysis
# ====================================================================

# Read in chum wild escapement data (infilled) by CU
ChumEscpDat <- read.csv("DataOut/WildEscape_w.Infill_ByCU.csv")
ChumEscpDat$MU <- "SC Chum" # Add MU column - FLAG this may not be necessary - can't find MU in any of the function scripts. Check with Kendra
ChumEscpDat$CU <- substr(ChumEscpDat$CU_Name, 1,1) # pull out CU ID from raw CU column name
ChumEscpDat$CU_Name <- substr(ChumEscpDat$CU_Name, 5, 100) # pull out just the name of the CU, replace CU_Name with that (remove the CU number)
# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
  colnames(ChumEscpDat)[colnames(ChumEscpDat)=="Year"] <- "yr"
  colnames(ChumEscpDat)[colnames(ChumEscpDat)=="SiteEsc"] <- "Escp" # FLAG: check that this is the right column to use (infilled escapement)

# Read in chum stock-recruit data
ChumSRDat <- read.csv("DataOut/SRdatWild.csv")
# Rename columns to be used in analysis, follow format of IFCohoStudy/DataIn/IFCoho_SRbyCU.csv
ChumSRDat$CU_ID <- substr(ChumSRDat$CU, 1,1) # make new column with just CU id number
ChumSRDat$CU_Name <- substr(ChumSRDat$CU, 5,100) # make new column with just CU name
ChumSRDat <- ChumSRDat[ , !names(ChumSRDat) %in% "CU"] # remove raw CU column
names(ChumSRDat)[names(ChumSRDat) =="Year"] <- "BroodYear" # FLAG: Check that Year is actually brood year in the infilling code
names(ChumSRDat)[names(ChumSRDat) =="Escape"] <- "Spawners" # FLAG: Check that Escape column is spawners (as opposed to the Return column)
names(ChumSRDat)[names(ChumSRDat) == "Recruit"] <- "Recruits" 
  
# FLAG: Should probably limit stock-recruit data to year > 1959/1960 to allow for full brood year returns up to age 6. 
# This may be done automatically, see retroFunctions.r line 20 
# This is done automatically using the BroodYrLag variable

# ==================================================================================
# Call functions to plot data availability:
# ====================================================================================
plot_CU_DataObs_Over_Time(ChumEscpDat, chumDir, plotName="SC_Chum_DataByCU")
plot_Num_CUs_Over_Time(ChumEscpDat, chumDir, plotName="SC_Chum_N_CUs")

# Note: these next 2 two escpt plots need to have formatting fixed
plot_CU_Escp_Over_Time(ChumEscpDat, chumDir, plotName="SCChum Esc", samePlot = T)
plot_CU_Escp_Over_Time(ChumEscpDat, chumDir, plotName="SCChum Esc Separate", samePlot = F)

# ==================================================================================================================
# Run retrospective analyses:
# =====================================================================================================================

TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)

# Run annual restrospective analyses using CUs ===========================

# Loop over p values and run annual retrospective analyses for each level of p
ps <- 0.90 # for now just use one p value
#ps <- c(seq(0.6, 0.95,.05), 0.99) 

# LW Notes
# Choosing values for runAnnualRetro function:
# - genYrs: average age of return.  runFraserCoho.r uses 3. I use 4 for chum (age 4 returns made up 69% of returns based on 2013 data). Update for 2018 age composition data.
# - BroodYrLag: How many years will it take between a fish returning and when you can construct the full recruitment.
#     The BroodYrLag parameter used when prepping for the LRP calculations represents the number of years 
#     it will take to reconstruct recruitment from a single brood year  (i.e., how many years of recruitment 
#     do you need to observe to estimate recruitment from a brood year).  Another way of defining it is the 
#     number of age classes that recruitment (maturation) occurs over. 
#     For the coho example, fish mostly return at age 3 or 4, so BroodYrLag is 2.
#     For the chum example, it looks like fish return at ages 3,4,5,6, so BroodYrLag is 4.
#     This parameter gets used in the runAnnualRetro() function to remove the first few years of data for 
#     which recruitment is NA because the BY has not yet been fully observed.

# not sure if this step is needed
ChumEscpDat <- ChumEscpDat[ChumEscpDat$yr >= 1958 ,] # remove years without full recruitment
ChumSRDat <- ChumSRDat[ChumSRDat$BroodYear >= 1958 ,] # remove years without full recruitment

source(paste0(codeDir, "/LRPFunctions.r"))

undebug(runAnnualRetro)
for(pp in 1:length(ps)){
  # Run with Binomial LRP model with individual model Ricker
  runAnnualRetro(EscpDat=ChumEscpDat, SRDat=ChumSRDat, startYr=1970, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
                 BMmodel = "SR_IndivRicker_NoSurv", LRPmodel="BinLogistic", integratedModel=T,
                 useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_",ps[pp]*100, sep=""),
                 bootstrapMode = F, plotLRP=T)
  # 
  # # Run with Bernoulli LRP model with individual model Ricker 
  # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurv_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
}


# Debugging
# Look at outputs after phase 1 fitting
# Note these are scaled
# read in estimates
All_Ests <- read.csv("2020-12-18_estimates_phase1_debugging.csv", stringsAsFactors = FALSE)
Scale <- 1000
# Plot to see what is going on with ricker parameters
N_Stocks <- length(unique(ChumEscpDat$CU_Name))

All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests$CU_ID[!(All_Ests$Param %in% c("logMuA", "logSigmaA", "B_0", "B_1", "Agg_LRP", "gamma", "logsigma", "prod", "Rec_Preds", "Logit_Preds"))] <- rep(1:(N_Stocks)) 

All_Ests$Estimate[All_Ests$Param %in% c( "logB", "logSigma")] <- exp(All_Ests$Estimate[All_Ests$Param %in% c( "logB", "logSigma")] )
#All_Ests$Param[All_Ests$Param == "logA"] <- "A"
All_Ests$Param[All_Ests$Param == "logB"] <- "B"
All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
All_Ests$Param[All_Ests$Param == "logSigma"] <- "sigma"
All_Ests$Param[All_Ests$Param == "SMSY"] <- "Smsy"
All_Ests[All_Ests$Param %in% c("Sgen", "SMSY", "SRep", "cap", "Agg_LRP", "^A$"), ] <-  All_Ests %>% filter(Param %in% c("Sgen", "SMSY", "SRep","cap","Agg_LRP")) %>% 
  mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
Preds <- All_Ests %>% filter(Param == "Logit_Preds")
Preds_Rec <- All_Ests %>% filter(Param == "Rec_Preds")
All_Ests <- All_Ests %>% filter(!(Param %in% c( "logSgen", "Logit_Preds", "Rec_Preds"))) 

ests_l <- All_Ests %>% filter(!is.na(CU_ID)) %>% select(Estimate, Param, CU_ID) %>% pivot_wider(names_from=Param, values_from=Estimate)

t <- merge(ChumSRDat, ests_l, by="CU_ID")

png("Figures/fig_debug_integrated_model_ricker_phase1.png", width=8, height=8, res=300, units="in")
layout(mat=matrix(1:9, byrow = TRUE, ncol=3))
for(i in 1:N_Stocks) {
  dat <- t[t$CU_ID==i,]
  plot(dat$Recruits ~ dat$Spawners, type="p", main=i)
  curve(dat$A * x * exp(- dat$B * x), add=TRUE) 
}
dev.off()

