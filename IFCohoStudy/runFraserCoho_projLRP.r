# ============================================================================
# Calculation of Projected Limit Reference Points for Interior Fraser Coho
# Kendra Holt, Last update: Apr 27, 2021
# ================================================================================

# Projected LRPs represent the aggregate abundance that is associated with a specified probability
# that a required proportion of CUs will be above their lower benchmarks in projections
#--  E.g., the LRP may represent the aggregate abundance that is projected to have a 50% probability
# that all CUs (100% of CUs) will be above their lower benchmarks

# Projections are done using the LRP branch of the samSim model (https://github.com/Pacific-salmon-assess/samSim)
# The code in this file is divided into the following sections:
#     (1) Read-in Coho data
#     (2) Specify initial parameters & data sets for projections 
#     (3) Run base projections (Using 4 different OM models at present)
#     (4) Run sensitivity analysis projections
#     (5) Estimate and save LRPs, and associated plots
#     (6) Plot CU-level spawner abundance projections (Optional)
#     (7) Make comparison plots among scenarios (NOT CURRENTLY WORKING)

# ===============================================================================


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
  source("get.mv.logistic.tau.r")
}
sourceAll()

# Load TMB models for fitting Bayesian Ricker stock recruit models;
#   outputs from these model fits are used to parameterize samSim

compile("TMB_Files/SR_IndivRicker_Surv_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_noLRP"))

compile("TMB_Files/SR_HierRicker_Surv_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv_noLRP"))

compile("TMB_Files/SR_IndivRicker_SurvCap_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap_noLRP"))

compile("TMB_Files/SR_HierRicker_SurvCap_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap_noLRP"))


# ======================================================================
#(1)  Read-in Coho data:  
# =====================================================================
setwd(cohoDir)

  CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")
 # Restrict data set to years 1998+ based on recommendation from Michael Arbeider
 CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)


# ======================================================================
# (2) Specify initial parameters & datasets for projections  
# =====================================================================

# Subset data up to current year
 year <- 2018 # - last year of data for parameterization
 BroodYrLag <- 2 # - number of years between brood year and first age of return (2 years for coho)
 
 # Only use SR data for brood years that have recruited by specified year
 # (note: most recent brood year is calculated by subtracting BroodYearLag (e.g. 2 years) from current year)
 SRDat <- CoSRDat %>%  filter(BroodYear <= year-BroodYrLag)
 SRDat$yr_num <- group_by(SRDat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
 SRDat$CU_ID <- group_by(SRDat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
 
 # 
 # # Calculate tau for error in age-at-maturity and save to input data file
 # CUages<-SRDat %>% select(Year = BroodYear, Age_3_Recruits, Recruits) %>% mutate(age3 = Age_3_Recruits/Recruits, age4=1-(Age_3_Recruits/Recruits))
 # CUages<-CUages%>%select(-Age_3_Recruits, -Recruits) %>% add_column(CU=SRDat$CU_ID, CU_Names=SRDat$CU_Name)
 # CU.tau <- NA
 # for (i in 1: length(unique(SRDat$CU_Name))){
 #   CUages.byCU <- CUages %>% filter(CU_Names== unique(CUages$CU_Names)[i]) %>% select(-c(Year, CU, CU_Names))
 #   # Added by K.Holt
 #   CUages.byCU[CUages.byCU$age3 == 1,1]<-0.99
 #   CUages.byCU[CUages.byCU$age4 == 0,2]<-0.01
 #   
 #   CU.tau[i] <- "get.mv.logistic.tau"(CUages.byCU)$best.tau
 # }
 # 
 
 
# TMB input parameters:
TMB_Inputs_HM <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                      logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5)


TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5)




# Prior means come from running "compareRickerModelTypes.r"
#cap_priorMean_HM<-c(10.957092, 5.565526, 11.467815, 21.104274, 14.803877)



# --- without bias correction:
cap_priorMean_HM_noBiasCor<-c(11.018110,  4.420246, 10.890888, 18.513363, 14.887224)
# --- with bias correction:
cap_priorMean_HM<-c(11.522086,  4.786252, 11.840564, 19.035260, 16.215605)



TMB_Inputs_HM_priorCap <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                               logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5,
                               cap_mean=cap_priorMean_HM, cap_sig=sqrt(2))

TMB_Inputs_HM_priorCap_noBiasCor <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                               logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5,
                               cap_mean=cap_priorMean_HM, cap_sig=sqrt(2))

# Prior means come from running "compareRickerModelTypes.r"
#cap_priorMean_IM<-c(11.153583,  5.714955, 11.535779, 21.379558, 14.889006)


# without bias correction
cap_priorMean_IM_noBiasCor<-c( 11.151905, 5.713514, 11.534271, 21.377327, 14.886351)
# --- with bias correction:
cap_priorMean_IM<-c(12.986952,  6.601032, 14.012006, 23.395788, 18.368697)


TMB_Inputs_IM_priorCap <- list(Scale = 1000, logA_Start = 1, Tau_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5,
                               cap_mean=cap_priorMean_IM, cap_sig=sqrt(2))

TMB_Inputs_IM_priorCap_noBiasCor <- list(Scale = 1000, logA_Start = 1, Tau_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5,
                               cap_mean=cap_priorMean_IM_noBiasCor, cap_sig=sqrt(2))

# Create output directories for Projected LRP outputs
figDir <- here(cohoDir, "Figures")
if (file.exists(figDir) == FALSE){
  dir.create(figDir)
}
figDir2 <- here(figDir, "ProjectedLRPs")
if (file.exists(figDir2) == FALSE){
  dir.create(figDir2)
}

projOutDir <- here(cohoDir, "DataOut")
if (file.exists(projOutDir) == FALSE){
  dir.create(projOutDir)
}
projOutDir2 <- here(projOutDir, "ProjectedLRPs")
if (file.exists(projOutDir2) == FALSE){
  dir.create(projOutDir2)
}


# ===================================================================
# (3) Run Base Projections
# ==================================================================

setwd(codeDir)
devtools::install_github("Pacific-salmon-assess/samSim", ref="LRP")




set.seed(100)



scenarioName <- "IM_mSurv_biasCorr"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=5000, nProj=1000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)




set.seed(100)
# 
# # Create samSim input files for current scenario
# scenarioName <- "Test_randomTrue"
# BMmodel <- "SR_IndivRicker_Surv"
# TMB_Inputs <- TMB_Inputs_IM
# projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
#                                 useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
#                                 nMCMC=10000, nProj=5000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE, 
#                                 biasCorrectEst=T, biasCorrectProj=T)


set.seed(100)


scenarioName <- "HM.Base"
BMmodel <- "SR_HierRicker_Surv"
TMB_Inputs <- TMB_Inputs_HM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=10000, nProj=5000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)




set.seed(100)

scenarioName <- "IMCap.Base"
BMmodel <- "SR_IndivRicker_SurvCap"
TMB_Inputs <- TMB_Inputs_IM_priorCap

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=10000, nProj=5000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)




set.seed(100)

scenarioName <- "HMCap.Base"

BMmodel <- "SR_HierRicker_SurvCap"
TMB_Inputs <- TMB_Inputs_HM_priorCap

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=1000, nProj=5000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)




# Create samSim input files for current scenario
scenarioName <- "IM.1.0cvER"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE)



# Create samSim input files for current scenario
scenarioName <- "IM.0.75cvER"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456*0.75, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE)





# Create samSim input files for current scenario
scenarioName <- "IM.0.25cvER"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456*0.25, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE)



# Create samSim input files for current scenario
scenarioName <- "IM.0.001cvER"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456*0.001, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE)




# Create samSim input files for current scenario
scenarioName <- "IM.1.0GammaSig"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=1,agePpnConst=TRUE)


# Create samSim input files for current scenario
scenarioName <- "IM.0.75GammaSig"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=0.75,agePpnConst=TRUE)

# Create samSim input files for current scenario
scenarioName <- "IM.0.50GammaSig"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=0.50,agePpnConst=TRUE)


# Create samSim input files for current scenario
scenarioName <- "IM.0.25GammaSig"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=0.25,agePpnConst=TRUE)




scenarioName <- "HM.Base"
BMmodel <- "SR_HierRicker_Surv"
TMB_Inputs <- TMB_Inputs_HM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=5000, nProj=2000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE)



scenarioName <- "IMCap.Base"
BMmodel <- "SR_IndivRicker_SurvCap"
TMB_Inputs <- TMB_Inputs_IM_priorCap

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE)


scenarioName <- "HMCap.Base"

BMmodel <- "SR_HierRicker_SurvCap"
TMB_Inputs <- TMB_Inputs_HM_priorCap

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=5000, nProj=2000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE)






# ==================================================================
# (4) Run Sensitivity Analyses 
# ====================================================================


# ===================================================================
#  (5) Estimate and Save LRPs (and associated plots)
# ==================================================================

# Specify threshold to use when calculating LRP
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThreshList<-c(0.50, 0.66, 0.90, 0.95) # probability theshhold; the LRP is set as the aggregate abundance that has this 
                    # probability that the propCUThreshold is met

OMsToInclude<-c("IM","IM_mSurv_biasCorr")

# Specify scenarios to calculate LRPs and make plots for.
# These scenarios will be looped over below with a LRP (and LRP plot) saved for each scenario

OMsToInclude<-c("IM.Base", "IM.0.25GammaSig","IM.0.50GammaSig","IM.0.75GammaSig","IM.1.0GammaSig","HM.Base",
                "IMCap.Base", "HMCap.Base","IM.0.001cvER","IM.0.25cvER","IM.0.75cvER","IM.1.0cvER")

OMsToInclude<-c("IM.Base","IM.Base_noBiasCorr","HM.Base","HM.Base_noBiasCorr",
                "IMCap.Base","IMCap.Base_noBiasCorr","HMCap.Base","HMCap.Base_noBiasCorr")

OMsToInclude<-c("Test_randomTrue", "HM.Base", "IMCap.Base", "HMCap.Base")

OMsToInclude<-c("Test_30Yrs_Thin", "Test_30Yrs_Thin_highER")

# Save LRPs for all OM scenarios
write.csv(LRP_Ests, paste(projOutDir2, "ProjectedLRPs.csv", sep="/"), row.names=F)
# Save LRP projection summaries used for calculating and plotting LRP (Optional)
write.csv(projLRPDat.plot, paste(projOutDir2, "ProjectedLRP_data.csv", sep="/"), row.names=F)


# ===================================================================================
# (5b) Estimate and Save Model-Averaged LRPs, and associated plots (Optional)
# =============================================================================

# Specify threshold to use when calculating LRP
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThreshList<-c(0.50, 0.66, 0.90, 0.95) # probability theshhold; the LRP is set as the aggregate abundance that has this 
# probability that the propCUThreshold is met

OMsToCombine<-c("IM.Base", "HM.Base", "IMCap.Base", "HMCap.Base")



for (p in 1:length(probThreshList)) {
  
  probThresh<-probThreshList[p]
  
  # Loop over OM Scenarios 
  for (i in 1:length(OMsToCombine)) {
    # Read in samSim outputs for OM
    filename<-paste("projLRPDat_",OMsToCombine[i],".csv",sep="")
    tmp<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
    tmp<-tmp %>% filter(year > max(SRDat$yr_num)+4)
    if (i == 1) {
      projLRPDat <- tmp
    } else {
      projLRPDat <- rbind(projLRPDat, tmp)
    }
  }
    # Create bins for projected spawner abundances
    minBreak<-0
    maxBreak<-round(max(projLRPDat$sAg),digits=-2)
    binSize<-200 # Note: bin size is currently set here
    breaks<-seq(minBreak, maxBreak,by=binSize)  
    
    # Set bin labels as the mid-point
    projLRPDat$bins<-cut(projLRPDat$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)))
    
    # Summarize nSims in each bin
    tmp<-projLRPDat %>% group_by(bins) %>% summarise(nSims=(length(ppnCUsLowerBM)))
    
    # Filter out bins with < 100 nSims
    tmp2<-projLRPDat %>% group_by(bins) %>% summarise(nSimsProp1=(length(ppnCUsLowerBM[ppnCUsLowerBM == propCUThresh]))) %>%
      add_column(nSims=tmp$nSims) %>% filter(nSims>=100)
    
    # For each bin, calculate probability that required proportion of CUs above benchmark
    projLRPDat<-tmp2 %>% add_column(prob=tmp2$nSimsProp1/tmp2$nSims)
    # For each bin, calculate the difference between the threshold probability and the calculated probability 
    projLRPDat$diff<-abs(probThresh-projLRPDat$prob)
    
    # Save projection summaries used to create plots
    projLRPDat$OM.Name<-"Combined_Base"
   # if (i == 1) projLRPDat.plot<-projLRPDat
  #  if (i > 1) projLRPDat.plot<-rbind(projLRPDat.plot,projLRPDat)
   
    projLRPDat.plot<-projLRPDat
     
    # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
    LRP<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))
    
    # Create a table of LRP estimates to be saved for each OM model
    if (p == 1) {
      LRP_Ests<-data.frame("Combined_Base", probThresh, propCUThresh, LRP, binSize)
      names(LRP_Ests)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
    } else {
      tmp.df<-data.frame("Combined_Base", probThresh, propCUThresh, LRP, binSize)
      names(tmp.df)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
      LRP_Ests<-rbind(LRP_Ests,tmp.df)
    }
    
    # Plot projected LRP abundance relationship =============================================================== 
    pdf(paste(cohoDir,"/Figures/ProjectedLRPs/", "Combined_Base", "_ProjLRPCurve_prob",probThresh,".pdf", sep=""), 
        width=6, height=6) 
    
    plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19, xlim=c(0,85000), cex=0.2,
         xlab="Aggregate Abundance", ylab="Pr (All CUs > Lower Benchmark)")
    abline(h=probThresh, lty=2)
    abline(v=LRP, col="orange", lwd=2)
    
    dev.off()
    
    # Option to plot histogram of nSims in each Agg Abundance Bin
    #barplot(height = projLRPDat$nSims,names.arg = projLRPDat$bins)
 
}


LRP_Ests_SingleModels <- read.csv(paste(projOutDir2, "ProjectedLRPs.csv", sep="/"))
LRP_Ests<-rbind(LRP_Ests_SingleModels,LRP_Ests)

# Save LRPs for all OM scenarios
write.csv(LRP_Ests, paste(projOutDir2, "ProjectedLRPs.csv", sep="/"), row.names=F)



# ===================================================================
# (6) Plot CU-level Spawner Abundance Projections (Optional)
# ==================================================================


OMsToInclude<-c("Test", "Test_BiasCorr")

for (i in 1:length(OMsToInclude)) {
 
  
  filename<-paste( "projSpwnDat_",OMsToInclude[i],".csv",sep="")
  spDat.i<-read.csv(here(cohoDir,"SamSimOutputs", "simData",filename))
  spDat.i$OM.Name<-OMsToInclude[i]
  # if (i == 1) projCUSpDat<-spDat.i
  # if (i > 1) projCUSpDat<-rbind(projCUSpDat,spDat.i)
  
  projCUSpDat.i<-spDat.i
  
  # Function to plot spawner abundance projections, by CU
  makeCUSpawnerProjPlot<- function(i, projSpwnDat, CUNames) {
    plotDat<-projSpwnDat %>% filter(CU==i) %>% group_by(year, expRate) %>% 
      summarise(medSpawn=median(spawners), lwr=quantile(spawners,0.10),upr=quantile(spawners,0.90))
    p <- ggplot(data=plotDat, mapping=aes(x=year,y=medSpawn, colour=factor(expRate))) +
      geom_ribbon(data=plotDat, aes(ymin = lwr, ymax = upr, x=year, fill=factor(expRate)), alpha=0.2) +
      geom_line(mapping=aes(x=year, y=medSpawn)) +
      geom_line(data=plotDat %>% filter(year < 18), col="black", size=1) +
      ggtitle(CUNames[i]) +
      xlab("Year") + ylab("Spawners") +
      theme_classic()  
  }
  
  ps<-lapply(1:length(unique(SRDat$CU_Name)), makeCUSpawnerProjPlot, projSpwnDat = projCUSpDat.i,CUNames=unique(SRDat$CU_Name))
  
  pdf(paste(cohoDir,"/Figures/ProjectedLRPs/", OMsToInclude[i], "_CUSpawnerProj.pdf", sep=""), 
      width=9, height=6)
  do.call(grid.arrange,  ps)
  dev.off()
  
}





# ===================================================================
# (7) Plots of LRP stabilization with number of trials
# ==================================================================

# Specify threshold to use when calculating LRP
# # Note: may want to loop over probThresholds as well; still needs to be added
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThresh<-0.5 # probability theshhold; the LRP is set as the aggregate abundance that has this
# probability that the propCUThreshold is met

OMsToInclude<-c("Test3")



# Read in samSim outputs for OM
filename<-paste("projLRPDat_",OMsToInclude,".csv",sep="")
projLRPDat<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
CUpars <- read.csv(paste(cohoDir, "SamSimInputs/CUPars.csv",sep="/"))

nTrials <- 2000

# Loop over nTrials
for (i in 1:500) {
  projLRPDat<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
  projLRPDat<-projLRPDat %>% filter(year > CUpars$ageMaxRec[1]*10) %>% filter(iteration < (i*10))
  
  # Create bins for projected spawner abundances
  minBreak<-0
  maxBreak<-round(max(projLRPDat$sAg),digits=-2)
  binSize<-200 # Note: bin size is currently set here
  breaks<-seq(minBreak, maxBreak,by=binSize)
  
  # Set bin labels as the mid-point
  projLRPDat$bins<-cut(projLRPDat$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)))
  
  # Summarize nSims in each bin
  tmp<-projLRPDat %>% group_by(bins) %>% summarise(nSims=(length(ppnCUsLowerBM)))
  
  # Filter out bins with < 100 nSims
  tmp2<-projLRPDat %>% group_by(bins) %>% summarise(nSimsProp1=(length(ppnCUsLowerBM[ppnCUsLowerBM == propCUThresh]))) %>%
    add_column(nSims=tmp$nSims) %>% filter(nSims>=10)
  
  # For each bin, calculate probability that required proportion of CUs above benchmark
  projLRPDat<-tmp2 %>% add_column(prob=tmp2$nSimsProp1/tmp2$nSims)
  # For each bin, calculate the difference between the threshold probability and the calculated probability
  projLRPDat$diff<-abs(probThresh-projLRPDat$prob)
  
  # Save projection summaries used to create plots
  projLRPDat$OM.Name<-OMsToInclude
  if (i == 1) projLRPDat.plot<-projLRPDat
  if (i > 1) projLRPDat.plot<-rbind(projLRPDat.plot,projLRPDat)
  
  # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
  LRP<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))
  
  # Create a table of LRP estimates to be saved for each OM model
  if (i ==1) {
    LRP_Ests_nTrials<-data.frame(OMsToInclude, i*10, probThresh, propCUThresh, LRP, binSize)
    names(LRP_Ests_nTrials)<-c("OM", "nTrials","ProbThresh", "PropCURequired", "LRP", "binSize")
  } else {
    tmp.df<-data.frame(OMsToInclude, i*10, probThresh, propCUThresh, LRP, binSize)
    names(tmp.df)<-c("OM", "nTrials", "ProbThresh", "PropCURequired", "LRP", "binSize")
    LRP_Ests_nTrials<-rbind(LRP_Ests_nTrials,tmp.df)
  }
  
}

# Save LRPs for all OM scenarios
write.csv(LRP_Ests_nTrials, paste(projOutDir2, "ProjectedLRPs_nTrials_p0.5.csv", sep="/"), row.names=F)
# Save LRP projection summaries used for calculating and plotting LRP (Optional)
LRP_nTrials <- as.data.frame(read.csv(paste(projOutDir2, "ProjectedLRPs_nTrials_p0.5.csv", sep="/")))

LRP_nTrials_plot <- ggplot(LRP_nTrials, aes(nTrials,LRP))+geom_line() #+
#ylim(14000,19000) #+
#geom_vline(xintercept=100, linetype="dashed")


ggsave(paste(cohoDir,"/Figures/LRP_ntrials_p0.5.png",sep=""), plot = LRP_nTrials_plot,
       width = 6, height = 4, units = "in")
#
# # Save LRPs for all OM scenarios
# write.csv(LRP_Ests, paste(projOutDir2, "ProjectedLRPs.csv", sep="/"), row.names=F)
# # Save LRP projection summaries used for calculating and plotting LRP (Optional)
# write.csv(projLRPDat.plot, paste(projOutDir2, "ProjectedLRP_data.csv", sep="/"), row.names=F)




# ===================================================================
# (8) Compare among-CU correlation in observed and projected escapements
# ==================================================================

# # ---- Code to compare observed and projected among-CU correlations in spawner abundance

OMsToTest<-c("IM.Base", "IM.0.25GammaSig","IM.0.50GammaSig","IM.0.75GammaSig","IM.1.0GammaSig")

OMsToTest<-c("IM.Base", "IM.0.25GammaSig","IM.0.50GammaSig","IM.0.75GammaSig","IM.1.0GammaSig", "HM.Base",
             "IMCap.Base", "HMCap.Base")


OMsToTest<-c("IM.0.001cvER","IM.0.25cvER","IM.Base","IM.0.75cvER","IM.1.0cvER")


OMsToTest<-c("Test", "Test_BiasCorr")

propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThresh<-0.5

for (j in 1:length(OMsToTest)) {

  filename<-paste( "projSpwnDat_",OMsToTest[j],".csv",sep="")
  spDat<-read.csv(here(cohoDir,"SamSimOutputs", "simData",filename))

  spDat<-as_tibble(spDat)
  spDat<-spDat%>%select(-X)

  RecCorMat<-array(NA,dim=c(5,5,max(spDat$iteration)))
  SpwnCorMat<-array(NA,dim=c(5,5,max(spDat$iteration)))

  for (i in 1:max(spDat$iteration)) {

    recruits.i<-spDat %>% filter(iteration==i & expRate==0.325) %>%  select(-spawners)
    recruits.i<-recruits.i %>% select(-expRate, -iteration)
    cor_mat<-recruits.i %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=recruits) %>% select(-year) %>% cor()
    RecCorMat[,,i]<-cor_mat

    spawners.i<-spDat %>% filter(iteration==i & expRate==0.125) %>%  select(-recruits)
    spawners.i<-spawners.i %>% select(-expRate, -iteration)
    cor_mat<-spawners.i %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=spawners) %>% select(-year) %>% cor()
    SpwnCorMat[,,i]<-cor_mat

    if (i ==1) SpwnCorrValues<-SpwnCorMat[,,i][lower.tri(SpwnCorMat[,,i])==TRUE]
    if (i > 1) SpwnCorrValues<-c(SpwnCorrValues,SpwnCorMat[,,i][lower.tri(SpwnCorMat[,,i])==TRUE])
  }
  
  OM_Name<-rep(OMsToTest[j],max(spDat$iteration))
  
  if (j ==1) SpwnCorr.df<-data.frame(OM_Name,SpwnCorrValues)
  if (j > 1) {
    tmp<-data.frame(OM_Name,SpwnCorrValues)
    SpwnCorr.df<-rbind(SpwnCorr.df,tmp)
  }
  
}


# Add observed escapement correlations to correlation data frame
spawners.obs<-data.frame(SRDat$BroodYear, SRDat$CU_ID, SRDat$Spawners)
names(spawners.obs)<-c("year", "CU", "spawners")
cor_mat<-spawners.obs %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=spawners) %>% select(-year) %>% cor()
SpwnCorrValues.Obs<-cor_mat[lower.tri(cor_mat)==TRUE]
tmp<-data.frame(OM_Name = "Observed",SpwnCorrValues = SpwnCorrValues.Obs)
SpwnCorr.df<-rbind(SpwnCorr.df,tmp)
  
# Save LRPs for all OM scenarios
write.csv(SpwnCorr.df, paste(projOutDir2, "SpwnCorr.df.csv", sep="/"), row.names=F)




dat<-as_tibble(SpwnCorr.df) %>% filter(OM_Name %in% c("Observed","IM.Base", "IM.0.25GammaSig","IM.0.50GammaSig","IM.0.75GammaSig",
                                                      "IM.1.0GammaSig"))


g <- ggplot(dat,aes(y=SpwnCorrValues,x=as.factor(OM_Name))) + geom_boxplot(width=0.5) +
  scale_x_discrete(limits=c("Observed","IM.Base", "IM.0.25GammaSig","IM.0.50GammaSig","IM.0.75GammaSig",
                            "IM.1.0GammaSig"),
                   labels=c("Obs","gammaSig=0","gammaSig=0.03" ,"gammaSig=0.06","gammaSig=0.09","gammaSig=0.12")) +
                    xlab("Sensitivity Analysis Scenario") + ylab("Between-CU Correlation")


dat<-as_tibble(SpwnCorr.df) %>% filter(OM_Name %in% c("Observed","IM.0.001cvER","IM.0.25cvER","IM.Base","IM.0.75cvER","IM.1.0cvER"))



g2 <- ggplot(dat,aes(y=SpwnCorrValues,x=as.factor(OM_Name))) + geom_boxplot(width=0.5) +
  scale_x_discrete(limits=c("Observed","IM.0.001cvER","IM.0.25cvER","IM.Base","IM.0.75cvER","IM.1.0cvER"),
                   labels=c("Obs", "~0cvER", "0.11cvER", "0.23cvER", "0.34cvER", "0.46cvER")) +
  xlab("Sensitivity Analysis Scenario") + ylab("Between-CU Correlation")

 ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareEscpCorrelation_cvER.png",sep=""), plot = g2,
        width = 4, height = 3, units = "in") 

 
 
 
 
 
 
dat<-as_tibble(SpwnCorr.df) %>% filter(OM_Name %in% c("Observed", "Test_BiasCorr"))
 
 
g3 <- ggplot(dat,aes(y=SpwnCorrValues,x=as.factor(OM_Name))) + geom_boxplot(width=0.5) +
   scale_x_discrete(limits=c("Observed","Test_BiasCorr"),
                    labels=c("Obs","Test_BiasCorr")) +
   xlab("Sensitivity Analysis Scenario") + ylab("Between-CU Correlation")
 
 
 
 
# =============================================================================================
# (9) Look at posterior samples from MCMC model parameterization
# =============================================================================================
 
 Mod <- "SR_IndivRicker_Surv_noLRP"
 OM<-"Test3"
 
 muLSurv<-SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
 muLSurv <- muLSurv$muLSurv
 
 post<-read.csv(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/SamSimInputs/",OM,"/","SR_IndivRicker_Surv","_mcmc.csv", sep=""))
 adjProd<-post$alpha + post$gamma * rep(muLSurv, length=nrow(post))
 post<- post %>% add_column(adjProd = adjProd)

 CU_list<-unique(SRDat[, "CU_Name"])
 CUID_list<-unique(SRDat[, "CU_ID"])
 nCUs<-length(CU_list)
 
 
 library(gridExtra)
 
 # ps<-list()
 # ps<-lapply(1:nCUs, plotPostHist, post=post, parName="Sgen", CUNames=CU_list)
 # png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",OM,"_SgenPost.png", sep=""))
 # do.call(grid.arrange, ps)
 # dev.off()
 
 ps<-list()
 ps<-lapply(1:nCUs, plotPostHist, post=post, parName="alpha",  CUNames=CU_list)
 png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",OM,"_AlphaPost.png", sep=""))
 do.call(grid.arrange, ps)
 dev.off()
 
 ps<-list()
 ps<-lapply(1:nCUs, plotPostHist, post=post, parName="sigma",  CUNames=CU_list)
 png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",OM,"_SigmaPost.png", sep=""))
 do.call(grid.arrange, ps)
 dev.off()
 
 ps<-list()
 ps<-lapply(1:nCUs, plotPostHist, post=post, parName="adjProd", CUNames=CU_list)
 png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",OM,"AdjProdPost.png", sep=""))
 do.call(grid.arrange, ps)
 dev.off()
 
 if(Mod %in% c("SR_HierRicker_SurvCap_noLRP", "SR_IndivRicker_SurvCap_noLRP","SR_HierRicker_Surv_noLRP", "SR_IndivRicker_Surv_noLRP")) {
   png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",OM,"GammaPost.png", sep=""))
   margPost<-post %>% filter(stk==1) %>% select("gamma")
   margPost<-as.data.frame(margPost$gamma)
   ggplot(margPost, aes(x=margPost)) + geom_histogram() + labs(x="gamma", y = "Count")
   dev.off()
 }
 
 if(Mod %in% c("SR_HierRicker_SurvCap_noLRP", "SR_IndivRicker_SurvCap_noLRP")) {
   ps<-list()
   ps<-lapply(1:nCUs, plotPostHist, post=post, parName="cap", CUNames=CU_list)
   png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",OM,"CapPost.png", sep=""))
   do.call(grid.arrange, ps)
   dev.off()
   
 } else {
   
   ps<-list()
   ps<-lapply(1:nCUs, plotPostHist, post=post, parName="beta", CUNames=CU_list)
   png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",OM,"BetaPost.png", sep=""))
   do.call(grid.arrange, ps)
   dev.off()
   
 }
 
 
 
 # =============================================================================================
 # (9b) Look at posterior samples from MCMC model parameterization, with a focus on those with low Sgen
 # =============================================================================================
 
 Mod <- "SR_HierRicker_Surv_noLRP"
 OM<-"Test"
 
 CU_list<-unique(SRDat[, "CU_Name"])
 x<-5
 
 muLSurv<-SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
 muLSurv <- muLSurv$muLSurv
 
 post<-read.csv(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/SamSimInputs/",OM,"/","SR_IndivRicker_Surv","_mcmc.csv", sep=""))
 #adjProd<-exp(post$alpha + post$gamma * rep(muLSurv, length=nrow(post)))
 adjProd<-post$alpha + post$gamma * rep(muLSurv, length=nrow(post))
 post<- post %>% add_column(adjProd = adjProd)
 
 post.all<-post %>% filter(stk == x)
 post.sub<-post %>% filter(stk == x, Sgen <= 340)
 #post.sub<-post %>% filter(stk == x, Sgen <= 560)
 
 
 # Plot Alphas ==================================
 margPost.all<-as.data.frame(post.all$alpha)
 margPost.sub<-as.data.frame(post.sub$alpha)
 
 h<-ggplot(margPost.all, aes(x=post.all$alpha)) + geom_histogram() + 
   labs(title=CU_list[x], x="alpha", y = "Count") +
   geom_histogram(data = margPost.sub, aes(x=post.sub$alpha),fill="red") 
 
 # Plot Betas ==================================
 margPost.all<-as.data.frame(post.all$beta)
 margPost.sub<-as.data.frame(post.sub$beta)
 
 h2<-ggplot(margPost.all, aes(x=post.all$beta)) + geom_histogram() + 
   labs(title=CU_list[x], x="beta", y = "Count") +
   geom_histogram(data = margPost.sub, aes(x=post.sub$beta),fill="red") 

 
 # Plot Sgens ==================================
 margPost.all<-as.data.frame(post.all$Sgen)
 margPost.sub<-as.data.frame(post.sub$Sgen)
 
 h3<-ggplot(margPost.all, aes(x=post.all$Sgen)) + geom_histogram() + 
   labs(title=CU_list[x], x="Sgen", y = "Count") +
   geom_histogram(data = margPost.sub, aes(x=post.sub$Sgen),fill="red") 
 
 
 # Plot adjProd ==================================
 margPost.all<-as.data.frame(post.all$adjProd)
 margPost.sub<-as.data.frame(post.sub$adjProd)
 
 h4<-ggplot(margPost.all, aes(x=post.all$adjProd)) + geom_histogram() + 
   labs(title=CU_list[x], x="adjProd", y = "Count") +
   geom_histogram(data = margPost.sub, aes(x=post.sub$adjProd),fill="red") 
 
 
 
 # Plot gamma ==================================
 margPost.all<-as.data.frame(post.all$gamma)
 margPost.sub<-as.data.frame(post.sub$gamma)
 
 h5<-ggplot(margPost.all, aes(x=post.all$gamma)) + geom_histogram() + 
   labs(title=CU_list[x], x="gamma", y = "Count") +
   geom_histogram(data = margPost.sub, aes(x=post.sub$gamma),fill="red") 
 
 
 
 # Plot Alpha vs Beta ==================================
 
 v<-ggplot(as.data.frame(post.all), aes(x=alpha, y=beta)) + geom_point() +
   geom_point(data=as.data.frame(post.sub), aes(x=alpha, y=beta), col="red")
 
 # Plot Adj Prod vs Beta ==================================
 
 v2<-ggplot(as.data.frame(post.all), aes(x=adjProd, y=beta)) + geom_point() +
   geom_point(data=as.data.frame(post.sub), aes(x=adjProd, y=beta), col="red")
 
 
 # Plot Alpha vs Gamma ==================================
 
 v3<-ggplot(as.data.frame(post.all), aes(x=alpha, y=gamma)) + geom_point() +
   geom_point(data=as.data.frame(post.sub), aes(x=alpha, y=gamma), col="red")

 
 
 # ===================================================================
 # (11) Make Comparison Plots Among Scenarios (NOT CURRENTLY WORKING)
 # ==================================================================
 
 
 # Create plots to compare LRP estimates 
 
 
 projLRPDat<-read.csv(file=paste(projOutDir2, "ProjectedLRPs.csv", sep="/"))
 
 
 plotDat<-projLRPDat %>% filter(OM %in% c("IM.Base", "HM.Base", "IMCap.Base", "HMCap.Base"))
 
 g<-ggplot(data=plotDat, aes(x=as.factor(ProbThresh), y=LRP, col=OM)) + geom_point(size=2) +
   xlab("Probability Threshold") + ylab("Estimated LRP") +
   labs(col = "SR Model")
 
 
 g2<-ggplot(data=plotDat, aes(x=OM, y=LRP, col=as.factor(ProbThresh))) + geom_point() +
   scale_x_discrete(limits=c("IM.Base","HM.Base","IMCap.Base","HMCap.Base"),
                    labels=c("IM", "HM", "IMcap", "HMcap")) +
   xlab("Stock Recruit Model") + ylab("Estimated LRP")
 
 
 
 # Note: The below code needs to be updated for new projected LRP method (Apr 26, 2021)
 
 
 # # Plot to compare LRP among different SRR structures =====================
 # 
 # OMsToPlot<-c("IM.Base", "HM.Base", "IMCap.base", "HMCap.base")
 # 
 # p<-1.0
 # plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
 # g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
 #              geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
 #              xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
 #               scale_x_discrete(limits=c("IM.Base","IMCap.base","HM.Base","HMCap.base"),
 #                                labels=c("IM", "IM.Cap", "HM", "HM.Cap"))
 # ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareLRP_byOM_p=",p,".pdf",sep=""), plot = g,
 #        width = 4, height = 3, units = "in")    
 # 
 # 
 # p<-0.8
 # plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
 # g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
 #   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
 #   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
 #   scale_x_discrete(limits=c("IM.Base","IMCap.base","HM.Base","HMCap.base"),
 #                    labels=c("IM", "IM.Cap", "HM", "HM.Cap"))
 # ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareLRP_byOM_p=",p,".pdf",sep=""), plot = g,
 #        width = 4, height = 3, units = "in")
 # 
 # 
 # # Plot to show combined OM scenarios =====================================
 # 
 # OMsToPlot<-c("IM.Base","IMCap.base", "CombinedIM")
 # 
 # p<-1.0
 # 
 # plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
 # 
 # g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
 #   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
 #   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
 #   scale_x_discrete(limits=c("IM.Base","IMCap.base", "CombinedIM"),
 #                    labels=c("IM", "IM.Cap", "IM.Composite"))
 # 
 # ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareLRP_compositeIM_p=",p,".pdf",sep=""), plot = g,
 #        width = 4, height = 3, units = "in")    
 # 
 # 
 # Plot to show sensitivity analysis to variability in gamma
 
 # OMsToPlot<-c("IM.Base","IM.medGammaSig","IM.medGammaSig")
 
 # p<-1.0
 # plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
 # g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
 #   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
 #   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
 #   scale_x_discrete(limits=c("IM.Base","IM.60RecCor","IM.40RecCor","IM.10RecCor", "IM-.20RecCor","IM-.40RecCor"),
 #                    labels=c("Base(MPD)", "60%Corr", "40%Corr","20%Corr","-20%Corr","-40%Corr"))
 # ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareLRP_SAnalRecCorr_p=",p,".pdf",sep=""), plot = g,
 #        width = 4, height = 3, units = "in")  
 # 
 # p<-0.8
 # plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
 # g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
 #   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
 #   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
 #   scale_x_discrete(limits=c("IM.Base","IM.80RecCor","IM.60RecCor","IM.40RecCor","IM.10RecCor"),
 #                    labels=c("Base(MPD)", "80%Corr", "60%Corr","40%Corr","10%Corr"))
 # ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareLRP_SAnalRecCorr_p=",p,".pdf",sep=""), plot = g,
 #        width = 4, height = 3, units = "in")  
 # 
 # 
 # 
 
 
 
 
  
# # Plot to show sensitivity analysis to ER variability
# 
# OMsToPlot<-c("IM.Base","IM.cvER1.5","IM.cvER2.0","IM.cvER2.5","IM.cvER3.0")
# 
# p<-1.0
# plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
#   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   scale_x_discrete(limits=c("IM.Base", "IM.cvER1.5","IM.cvER2.0", "IM.cvER2.5","IM.cvER3.0"),
#                    labels=c("Base", "1.5CV","2.0CV", "2.5CV","3.0CV"))
# ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareLRP_SAnalCvER_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in") 
# 
# p<-0.8
# plotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50)) +
#   geom_point() + geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   scale_x_discrete(limits=c("IM.Base", "IM.cvER1.5","IM.cvER2.0", "IM.cvER2.5","IM.cvER3.0"),
#                    labels=c("Base", "1.5CV","2.0CV", "2.5CV","3.0CV"))
# ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareLRP_SAnalCvER_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in") 
# 
# 
# 
# # Plot to compare projected and empirical LRPs
# 
# OMsToPlot<-c("IM.Base", "HM.Base", "IMCap.base", "HMCap.base")
# 
# p<-1.0
# projPlotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# projPlotDat <- projPlotDat %>% add_column(Method=rep("Projected", 4))
# 
# # Empirical LRPs with p = 0.99 and likelihood penalty
# HM.Cap<-c(32962.96235,	28565.95578,	37359.96893)
# IM<-c(24358.94812,	20298.83689,	28419.05935)
# IM.Cap<-c(33787.86878,	29323.28802,	38252.44954)
# HM<-c(19055.00439,	15281.82266,	22828.18612)
# 
# empiricalPlotDat<-data.frame(OM.Name=c("HM.Base", "HMCap.base","IM.Base", "IMCap.base"),
#                              ppnCUsLowerBM = rep(1,4), 
#                              LRP.50=rep(NA,4),LRP.95=rep(NA,4),LRP.05=rep(NA,4),
#                              Method = "Empirical")
# 
# empiricalPlotDat[1,3:5] <- HM
# empiricalPlotDat[2,3:5] <- HM.Cap
# empiricalPlotDat[3,3:5] <- IM
# empiricalPlotDat[4,3:5] <- IM.Cap
# 
# empiricalPlotDat<-as_tibble(empiricalPlotDat)
# 
# plotDat<-bind_rows(projPlotDat, empiricalPlotDat)
#   
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50, colour=as.factor(Method))) +
#   geom_point(position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0,position=position_dodge(width=0.5)) +
#   scale_x_discrete(limits=c("IM.Base","IMCap.base","HM.Base","HMCap.base"),
#                    labels=c("IM", "IM.Cap", "HM", "HM.Cap")) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   labs(colour = "Method") 
#   
# 
# 
# ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareProjVsEmp_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")  
#    
# # ---- with p = 0.80
# p<-0.80
# projPlotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# projPlotDat <- projPlotDat %>% add_column(Method=rep("Proj", 4))
# 
# 
# # Empirical LRPs with p = 0.80 and likelihood penalty
# HM.Cap<-c(22246.52014, 19686.91885, 24806.12142)
# IM<-c(16702.85366, 14317.14271, 19088.56461)
# IM.Cap<-c(22765.51884, 20171.41927,25359.61841)
# HM<-c(13624.39163,	11523.63963,	15725.14363)
# 
# empiricalPlotDat<-data.frame(OM.Name=c("HM.Base", "HMCap.base","IM.Base", "IMCap.base"),
#                              ppnCUsLowerBM = rep(1,4), 
#                              LRP.50=rep(NA,4),LRP.95=rep(NA,4),LRP.05=rep(NA,4),
#                              Method = "Empirical")
# 
# empiricalPlotDat[1,3:5] <- HM
# empiricalPlotDat[2,3:5] <- HM.Cap
# empiricalPlotDat[3,3:5] <- IM
# empiricalPlotDat[4,3:5] <- IM.Cap
# 
# empiricalPlotDat<-as_tibble(empiricalPlotDat)
# 
# plotDat<-bind_rows(projPlotDat, empiricalPlotDat)
# 
# g<-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50, colour=as.factor(Method))) +
#   geom_point(position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0,position=position_dodge(width=0.5)) +
#   scale_x_discrete(limits=c("IM.Base","IMCap.base","HM.Base","HMCap.base"),
#                    labels=c("IM", "IM.Cap", "HM", "HM.Cap")) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   labs(colour = "Method") 
# 
# ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareProjVsEmp_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in") 
# 
# 
# # Empirical regression using p = 0.999
# p<-1.0
# projPlotDat<- LRPs %>% filter(ppnCUsLowerBM ==p & OM.Name %in% OMsToPlot)
# projPlotDat <- projPlotDat %>% add_column(Method=rep("Projected", 4))
# 
# # Empirical LRPs with p = 0.99 and likelihood penalty
# HM.Cap<-c(32962.96235,	28565.95578,	37359.96893)
# IM<-c(24358.94812,	20298.83689,	28419.05935)
# IM.Cap<-c(33787.86878,	29323.28802,	38252.44954)
# HM<-c(19055.00439,	15281.82266,	22828.18612)
# 
# empiricalPlotDat<-data.frame(OM.Name=c("HM.Base", "HMCap.base","IM.Base", "IMCap.base"),
#                              ppnCUsLowerBM = rep(1,4), 
#                              LRP.50=rep(NA,4),LRP.95=rep(NA,4),LRP.05=rep(NA,4),
#                              Method = "Empirical")
# 
# empiricalPlotDat[1,3:5] <- HM
# empiricalPlotDat[2,3:5] <- HM.Cap
# empiricalPlotDat[3,3:5] <- IM
# empiricalPlotDat[4,3:5] <- IM.Cap
# 
# empiricalPlotDat<-as_tibble(empiricalPlotDat)
# 
# plotDat<-bind_rows(projPlotDat, empiricalPlotDat)
# 
# g <-ggplot(data=plotDat, mapping=aes(x=OM.Name, y=LRP.50, colour=as.factor(Method))) +
#   geom_point(position=position_dodge(width=0.5)) + 
#   geom_errorbar(aes(ymin=LRP.05, ymax=LRP.95), width=0,position=position_dodge(width=0.5)) +
#   scale_x_discrete(limits=c("IM.Base","IMCap.base","HM.Base","HMCap.base"),
#                    labels=c("IM", "IM.Cap", "HM", "HM.Cap")) +
#   xlab("Operating Model") + ylab(paste("LRP (p=", p,")", sep="")) +
#   labs(colour = "Method") 
# 
# 
# 
# ggsave(paste(cohoDir,"/Figures/ProjectedLRPs/compareProjVsEmp_p=",p,".pdf",sep=""), plot = g,
#        width = 4, height = 3, units = "in")  
# 
# 
# 
# # Plot correlation matrix to show scenarios
# library(ggcorrplot)
# cormat<-read.csv(paste(cohoDir,"/SamSimInputs/IM.base/cohoCorrMat.csv", sep=""), header=F)
# rownames(cormat) <- unique(SRDat$CU_Name)
# colnames(cormat) <- unique(SRDat$CU_Name)
# ggcorrplot(cormat,
#            hc.order = TRUE,
#            type = "lower",
#            outline.color = "white", 
#            lab=T)
# 
# 
# 
# 
# # Plot ER distribution
# 
# canER<-0.125
# cvER<-0.456
# sigCanER<-cvER*canER
# 
# shape1<- canER^2 * (((1-canER)/sigCanER^2)-(1/canER))
# shape2<-shape1 * (1/canER-1)
# 
# sampBeta<-rbeta(1000,shape1,shape2)
