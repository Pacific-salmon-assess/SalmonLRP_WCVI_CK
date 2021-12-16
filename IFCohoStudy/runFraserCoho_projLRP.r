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
#        (5.1) Estimate and Save Model-Averaged LRPs, and associated plots (Optional)
#        (5.2) Plot LRP estimates relative to aggregate abundance time series
#     (6) Plot CU-level spawner abundance projections (Optional)
#     (7) Plots of LRP stabilization with number of trials
#     (8) Compare among-CU correlation in observed and projected escapements
#     (9) Look at posterior samples from MCMC model parameterization
#     (10) Make example plot to show calculation of projected LRPs
#     (11) Plot Correlation Matrices for Recruitment Residuals
#     (12) Make Projected Curve Comparison Plots Among Scenarios
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

 
 CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
 # Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
 colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
 colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
 colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
 colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"
 
 AggEscp <- CoEscpDat %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>%
   mutate(Gen_Mean = rollapply(Agg_Escp, 3, gm_mean, fill = NA, align="right"))
 
 
 # Summary of CU-level escapements based on Natural.Spawners
 # -- Read-in Natural spawners at the stream level
 CoEscpDat_bySubpop<-read.csv("DataIn/IFCoho_escpBySubpop.csv")
 # --Change column names to yr, CU_Name, Escp, Subpop_Name, CU_ID
 CoEscpDat_bySubpop<-CoEscpDat_bySubpop %>% select(MU_Name=MU_Name, yr=Return.Year, CU_Name=Conservation.Unit, Escp=Natural.Returns, Subpop_Name=Sub.Population)
 tmp.df<-data.frame(CU_Name=unique(CoEscpDat_bySubpop$CU_Name), CU_ID=seq(1,length(unique(CoEscpDat_bySubpop$CU_Name)),by=1))
 CoEscpDat_bySubpop <- left_join(CoEscpDat_bySubpop,tmp.df)
 # -- Calculate CU-level Natural Spawners
 CoEscpDat.Natural <- as_tibble(CoEscpDat_bySubpop) %>% group_by(MU_Name, CU_Name, CU_ID, yr) %>% select(MU_Name, CU_Name, CU_ID, yr, Escp) %>%summarize(Escp=sum(Escp))
 CoEscpDat.Natural <- as.data.frame(CoEscpDat.Natural)
 # -- Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
 colnames(CoEscpDat.Natural)[colnames(CoEscpDat.Natural)=="CU_ID"] <- "CU"
 colnames(CoEscpDat.Natural)[colnames(CoEscpDat.Natural)=="MU_Name"] <- "MU"
 
 AggEscp.Natural <- CoEscpDat.Natural %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>%
   mutate(Gen_Mean = rollapply(Agg_Escp, 3, gm_mean, fill = NA, align="right"))
 
 
 
# ======================================================================
# (2) Specify initial parameters & datasets for projections  
# =====================================================================

# Subset data up to current year ====================================================
 year <- 2020 # - last year of data for parameterization
 BroodYrLag <- 4 # - number of years between brood year and first age of return (2 years for coho)
 
 # Only use SR data for brood years that have recruited by specified year
 # (note: most recent brood year is calculated by subtracting BroodYearLag (e.g. 2 years) from current year)
 SRDat <- CoSRDat %>%  filter(BroodYear <= year-BroodYrLag)
 SRDat$yr_num <- group_by(SRDat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
 SRDat$CU_ID <- group_by(SRDat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
 # 

 # Caluclate proportion-at-age parameters directly from data and input to CUpars.csv ============
 CUages<-SRDat %>% select(Year = BroodYear, Age_3_Recruits, Recruits) %>% mutate(age3 = Age_3_Recruits/Recruits, age4=1-(Age_3_Recruits/Recruits))
 CUages<-CUages%>%select(-Age_3_Recruits, -Recruits) %>% add_column(CU=SRDat$CU_ID, CU_Names=SRDat$CU_Name)
 # CU.meanAge3 <- NA
 # CU.meanAge4 <- NA
 # CU.tau <- NA
 # for (i in 1: length(unique(SRDat$CU_Name))){
 #   CUages.byCU <- CUages %>% filter(CU_Names== unique(CUages$CU_Names)[i]) %>% select(-c(Year, CU, CU_Names))
 #   # Added by K.Holt
 #   CUages.byCU[CUages.byCU$age3 == 1,1]<-0.99
 #   CUages.byCU[CUages.byCU$age4 == 0,2]<-0.01
 #   CU.meanAge3[i] <- mean(CUages.byCU$age3)
 #   CU.meanAge4[i] <- mean(CUages.byCU$age4)
 #   CU.tau[i] <- "get.mv.logistic.tau"(CUages.byCU)$best.tau
 # }
 # 
 # CUPars<-read.csv("SamSimInputs/CUPars.csv")
 # CUPars$tauCycAge<-CU.tau
 # CUPars$meanRec3<-CU.meanAge3
 # CUPars$meanRec4<-CU.meanAge4
 # write.csv(CUPars, "SamSimInputs/CUPars.csv", row.names=F)
 # 
 # Save plot of annual age proportions by CU
 plotAgeProp_byCU(CUages, outDir = paste(cohoDir,"Figures",sep="/"), plotName="coho-ObsAgeProp-byCU")

 
 
# Specify TMB input parameters for SR model fitting ==========================
# TMB input parameters:
TMB_Inputs_HM <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                      logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5)


TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5)


# Prior means come from running "compareRickerModelTypes_onlySR.r", with bias correction
## Note: using expansion by 1.40 for prior cap (not 1.5, like Arbeider et al., Korman et al)
cap_priorMean_IM<-c(11.661606, 4.220915, 13.584558, 20.156285, 17.128963)

TMB_Inputs_IM_priorCap <- list(Scale = 1000, logA_Start = 1, Tau_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5,
                               cap_mean=cap_priorMean_IM, cap_sig=sqrt(2))


# Prior means come from running "compareRickerModelTypes_onlySR.r"
## Note: using expansion by 1.40 for prior cap (not 1.5, like Arbeider et al., Korman et al)
cap_priorMean_HM<-  c(10.953847, 4.644284, 12.856508, 19.590172, 16.087002)

TMB_Inputs_HM_priorCap <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                               logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.5,
                               cap_mean=cap_priorMean_HM, cap_sig=sqrt(2))



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


scenarioName <- "Ricker"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)


scenarioName <- "Ricker_priorCap"
BMmodel <- "SR_IndivRicker_SurvCap"
TMB_Inputs <- TMB_Inputs_IM_priorCap

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1,gammaSigScalar=NULL, agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)


scenarioName <- "RickerHM"
BMmodel <- "SR_HierRicker_Surv"
TMB_Inputs <- TMB_Inputs_HM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1,gammaSigScalar=NULL, agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)





scenarioName <- "RickerHM_priorCap"
BMmodel <- "SR_HierRicker_SurvCap"
TMB_Inputs <- TMB_Inputs_HM_priorCap

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1,gammaSigScalar=NULL, agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)




# ==================================================================
# (4) Run Sensitivity Analyses 
# ====================================================================


scenarioName <- "Ricker_0.25cvER"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.25, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)



scenarioName <- "Ricker_0.75cvER"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=10000, nProj=20000, cvER=0.442*0.75, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)


scenarioName <- "Ricker_1.0cvER"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=10000, nProj=20000, cvER=0.442, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)


scenarioName <- "Ricker_1.0sigGamma"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1, gammaSigScalar=1.0, agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)



scenarioName <- "Ricker_0.75sigGamma"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1, gammaSigScalar=0.75, agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)


scenarioName <- "Ricker_0.5sigGamma"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1, gammaSigScalar=0.5, agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)

scenarioName <- "Ricker_0.25sigGamma"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE, 
                                recCorScalar=1, gammaSigScalar=0.25, agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T)


scenarioName <- "Ricker_ER2.5"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE,
                                recCorScalar=1, gammaSigScalar=0.25, agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T, ER=0.025)

scenarioName <- "Ricker_ER22.5"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE,
                                recCorScalar=1, gammaSigScalar=0.25, agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T, ER=0.225)


scenarioName <- "Ricker_ER32.5"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=14000, nProj=20000, cvER=0.442*0.5, cvERSMU = 0.442, annualcvERCU=FALSE,
                                recCorScalar=1, gammaSigScalar=0.25, agePpnConst=TRUE,
                                biasCorrectEst=T, biasCorrectProj=T, ER=0.325)




# ===================================================================
#  (5) Estimate and Save LRPs (and associated plots)
# ==================================================================

# Specify threshold to use when calculating LRP
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThreshList<-c(0.50, 0.66, 0.90, 0.99) # probability theshhold; the LRP is set as the aggregate abundance that has this 
                    # probability that the propCUThreshold is met

OMsToInclude<-c("Ricker", "Ricker_priorCap", "RickerHM", "RickerHM_priorCap",
                "Ricker_0.25cvER", "Ricker_0.75cvER","Ricker_1.0cvER",
                "Ricker_0.25sigGamma","Ricker_0.5sigGamma","Ricker_0.75sigGamma","Ricker_1.0sigGamma",
                "Ricker_ER2.5","Ricker_ER22.5","Ricker_ER32.5")


# Specify scenarios to calculate LRPs and make plots for.
# These scenarios will be looped over below with a LRP (and LRP plot) saved for each scenario


for (p in 1:length(probThreshList)) {
  
  probThresh<-probThreshList[p]
  
  # Loop over OM Scenarios 
  for (i in 1:length(OMsToInclude)) {
    
    # Read in samSim outputs for OM
    filename<-paste("projLRPDat_",OMsToInclude[i],".csv",sep="")
    projLRPDat<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
    #projLRPDat<-projLRPDat %>% filter(year > max(SRDat$yr_num)+4)
    
    projLRPDat<-projLRPDat %>% filter(year > max(SRDat$yr_num) + 4)
    
    # Create bins for projected spawner abundances
    binSize<-200 # Note: bin size is currently set here
    minBreak<-0
    maxBreak<-round(max(projLRPDat$sAg),digits=-2) + binSize
    breaks<-seq(minBreak, maxBreak,by=binSize)  
    
    # Set bin labels as the mid-point
    projLRPDat$bins<-cut(projLRPDat$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)),include.lowest = T)
    
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
    projLRPDat$OM.Name<-OMsToInclude[i]
    if (i == 1) projLRPDat.plot<-projLRPDat
    if (i > 1) projLRPDat.plot<-rbind(projLRPDat.plot,projLRPDat)
    
    # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
    LRP<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))
    
    # Create a table of LRP estimates to be saved for each OM model
    if (i ==1 & p == 1) {
      LRP_Ests<-data.frame(OMsToInclude[i], probThresh, propCUThresh, LRP, binSize)
      names(LRP_Ests)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
    } else {
      tmp.df<-data.frame(OMsToInclude[i], probThresh, propCUThresh, LRP, binSize)
      names(tmp.df)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
      LRP_Ests<-rbind(LRP_Ests,tmp.df)
    }
    
    # Plot projected LRP abundance relationship =============================================================== 
    pdf(paste(cohoDir,"/Figures/ProjectedLRPs/", OMsToInclude[i], "_ProjLRPCurve_prob",probThresh,".pdf", sep=""), 
        width=6, height=6) 
    
    # plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19, xlim=c(0,85000), ylim=c(0,1.0),cex=0.2,
    #      xlab="Aggregate Abundance", ylab="Pr (All CUs > Lower Benchmark)")
    plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19, ylim=c(0,1.0),cex=0.2,
         xlab="Aggregate Abundance", ylab="Pr (All CUs > Lower Benchmark)")
    abline(h=probThresh, lty=2)
    abline(v=LRP, col="orange", lwd=2)
    
    dev.off()
    
    # Option to plot histogram of nSims in each Agg Abundance Bin
    #barplot(height = projLRPDat$nSims,names.arg = projLRPDat$bins)
    
  }
  
}


# Save LRPs for all OM scenarios
write.csv(LRP_Ests, paste(projOutDir2, "ProjectedLRPs.csv", sep="/"), row.names=F)
# Save LRP projection summaries used for calculating and plotting LRP (Optional)
write.csv(projLRPDat.plot, paste(projOutDir2, "ProjectedLRP_data.csv", sep="/"), row.names=F)



# ===================================================================================
# (5.1) Estimate and Save Model-Averaged LRPs, and associated plots (Optional)
# =============================================================================

# Specify threshold to use when calculating LRP
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThreshList<-c(0.50, 0.66, 0.90, 0.95) # probability theshhold; the LRP is set as the aggregate abundance that has this 
# probability that the propCUThreshold is met

OMsToCombine<-c("Ricker", "Ricker_priorCap")



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




# =============================================================================
# (5.2) Plot LRP estimates relative to observed aggregate abundance time series
# =============================================================================


#AggEscp<-AggEscp %>% filter(yr >= 2000)

AggEscp.Natural<-AggEscp.Natural %>% filter(yr >= 2000)

LRP_Ests <- read.csv(paste(projOutDir2, "ProjectedLRPs.csv", sep="/"))

colList<-c("#E69F00", "#56B4E9", "#009E73", "#D55E00")

LRP_Ricker0.5<- LRP_Ests[LRP_Ests$OM=="Ricker" & LRP_Ests$ProbThresh==0.5,]$LRP
LRP_RickerCap0.5<- LRP_Ests[LRP_Ests$OM=="Ricker_priorCap" & LRP_Ests$ProbThresh==0.5,]$LRP

LRP_Ricker0.66<- LRP_Ests[LRP_Ests$OM=="Ricker" & LRP_Ests$ProbThresh==0.66,]$LRP
LRP_RickerCap0.66<- LRP_Ests[LRP_Ests$OM=="Ricker_priorCap" & LRP_Ests$ProbThresh==0.66,]$LRP

png(paste(cohoDir,"/Figures/coho-EscpSeries-wProjLRP.png", sep=""), width=450, height=350)

par(mar=c(4,4,1,1))

plot(AggEscp.Natural$yr, AggEscp.Natural$Gen_Mean, typ="l", ylim=c(0,max(AggEscp.Natural$Gen_Mean)), bty="l",
     xlab="Year", ylab="Aggregate Spawner Abundance", cex.lab=1.3, lwd=2, cex.axis=1.1)
points(AggEscp.Natural$yr, AggEscp.Natural$Gen_Mean,pch=19,col="grey30", cex=1)


abline(h=LRP_Ricker0.5, col=colList[1],lty=1, lwd=2)
#abline(h=LRP_Ricker0.66, col=colList[2],lty=1)
abline(h=LRP_RickerCap0.5, col=colList[1],lty=2,lwd=2)
#abline(h=LRP_RickerCap0.66, col=colList[2],lty=2)

# legend(x=2000, y=15000, cex=0.8,
#        legend=c("Ricker_p0.5", "Ricker_p0.66", "RickerCap_p0.5", "RickerCap_p0.66"),
#        col = c(colList[1], colList[2], colList[1], colList[2]), 
#        lty=c(1,1,2,2), bty="n")

legend(x=2000, y=10000,cex=1.3,
       legend=c("Ricker", "Ricker_priorCap"),
       col = c(colList[1], colList[1]),
       lty=c(1,2), bty="n", lwd=c(2,2))

dev.off()


# ===================================================================
# (6) Plot CU-level Spawner Abundance Projections (Optional)
# ==================================================================


OMsToInclude<-c("Ricker", "Ricker_priorCap")

for (i in 1:length(OMsToInclude)) {
 
  
  filename<-paste( "projSpwnDat_",OMsToInclude[i],".csv",sep="")
  spDat.i<-read.csv(here(cohoDir,"SamSimOutputs", "simData",filename))
  spDat.i$OM.Name<-OMsToInclude[i]
  # if (i == 1) projCUSpDat<-spDat.i
  # if (i > 1) projCUSpDat<-rbind(projCUSpDat,spDat.i)
  
  projCUSpDat.i<-spDat.i
  
  # Function to plot spawner abundance projections, by CU
  makeCUSpawnerProjPlot<- function(i, projSpwnDat, CUNames) {
  
    CU.name<-NA
     
    if(as.character(CUNames[i]) == "Middle_Fraser") CU.name<-"Middle Fraser"
    if(as.character(CUNames[i]) == "South_Thompson") CU.name<-"South Thompson"
    if(as.character(CUNames[i]) == "North_Thompson") CU.name<-"North Thompson"
    if(as.character(CUNames[i]) == "Lower_Thompson") CU.name<-"Lower Thompson"
    if(as.character(CUNames[i]) == "Fraser_Canyon") CU.name<-"Fraser Canyon"
    
    plotDat<-projSpwnDat %>% filter(CU==i) %>% group_by(year) %>% 
      summarise(medSpawn=median(spawners), lwr=quantile(spawners,0.10),upr=quantile(spawners,0.90))
    # plotDat<-projSpwnDat %>% filter(CU==i) %>% group_by(year, expRate) %>% 
    #   summarise(medSpawn=median(spawners), lwr=quantile(spawners,0.10),upr=quantile(spawners,0.90))
    #p <- ggplot(data=plotDat, mapping=aes(x=year,y=medSpawn, colour=factor(expRate))) +
     p <- ggplot(data=plotDat, mapping=aes(x=year,y=medSpawn)) +
      geom_ribbon(data=plotDat, aes(ymin = lwr, ymax = upr, x=year), alpha=0.2) +
      geom_line(mapping=aes(x=year, y=medSpawn)) +
      geom_line(data=plotDat %>% filter(year < 18), col="black", size=1) +
      ggtitle(CU.name) +
      xlab("Year") + ylab("Spawners") +
      theme_classic() +
       theme(axis.text=element_text(size=20), #change font size of axis text
           axis.title=element_text(size=20), #change font size of axis titles
           plot.title=element_text(size=20)) #change font size of plot title

     #text=element_text(size=20)
          
  }
  
  ps<-lapply(1:length(unique(SRDat$CU_Name)), makeCUSpawnerProjPlot, projSpwnDat = projCUSpDat.i,CUNames=unique(SRDat$CU_Name))
  
  png(paste(cohoDir,"/Figures/ProjectedLRPs/coho-CUSpawnerProj-",OMsToInclude[i],".png",sep=""),res=80, 
      width=850, height=750)
  do.call(grid.arrange,  ps)
  dev.off()
  
}





# ===================================================================
# (7) Plots of LRP stabilization with number of trials
# ==================================================================

# Specify threshold to use when calculating LRP
# # Note: may want to loop over probThresholds as well; still needs to be added
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThresh<-0.99 # probability theshhold; the LRP is set as the aggregate abundance that has this
# probability that the propCUThreshold is met

OMsToInclude<-c("test")



# Read in samSim outputs for OM
filename<-paste("projLRPDat_",OMsToInclude,".csv",sep="")
projLRPDat<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
CUpars <- read.csv(paste(cohoDir, "SamSimInputs/CUPars.csv",sep="/"))

nTrials <- 50000
nTrials.binSize<-5000

minTrials<-5000

nTrialsSeq<-seq(minTrials,nTrials,by=nTrials.binSize)

# Loop over nTrials
for (i in 1:length(nTrialsSeq)) {
  projLRPDat<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
  projLRPDat<-projLRPDat %>% filter(year > (max(SRDat$yr_num)+CUpars$ageMaxRec[1])) %>% filter(iteration < nTrialsSeq[i])
  
  # Create bins for projected spawner abundances
  minBreak<-0
  maxBreak<-round(max(projLRPDat$sAg),digits=-2)
  binSize<-200 # Note: bin size is currently set here
  breaks<-seq(minBreak, maxBreak+binSize,by=binSize)
  
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
  projLRPDat$OM.Name<-OMsToInclude
  if (i == 1) projLRPDat.plot<-projLRPDat
  if (i > 1) projLRPDat.plot<-rbind(projLRPDat.plot,projLRPDat)
  
  # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
  LRP<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))
  
  # Create a table of LRP estimates to be saved for each OM model
  if (i ==1) {
    LRP_Ests_nTrials<-data.frame(OMsToInclude, nTrialsSeq[i], probThresh, propCUThresh, LRP, binSize)
    names(LRP_Ests_nTrials)<-c("OM", "nTrials","ProbThresh", "PropCURequired", "LRP", "binSize")
  } else {
    tmp.df<-data.frame(OMsToInclude, nTrialsSeq[i], probThresh, propCUThresh, LRP, binSize)
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


ggsave(paste(cohoDir,"/Figures/LRP_ntrials_p",probThresh,".png",sep=""), plot = LRP_nTrials_plot,
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

OMsToTest<-c("Ricker", "Ricker_0.5sigGamma","Ricker_0.75sigGamma",
             "Ricker_1.0sigGamma")

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

    recruits.i<-spDat %>% filter(iteration==i) %>%  select(-spawners)
    recruits.i<-recruits.i %>% select(-expRate, -iteration)
    cor_mat<-recruits.i %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=recruits) %>% select(-year) %>% cor()
    RecCorMat[,,i]<-cor_mat

    spawners.i<-spDat %>% filter(iteration==i) %>%  select(-recruits)
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
 

# Calculate observed escapement correlations for natural spawners (but, don't add to data frame)
natSpawners.obs<-data.frame(CoEscpDat.Natural$yr, CoEscpDat.Natural$CU, CoEscpDat.Natural$Escp)
names(natSpawners.obs)<-c("year", "CU", "spawners")
cor_mat_nat<-natSpawners.obs %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=spawners) %>% select(-year) %>% cor()
natSpwnCorrValues.Obs<-cor_mat_nat[lower.tri(cor_mat_nat)==TRUE]



tmp<-data.frame(OM_Name = "Observed",SpwnCorrValues = SpwnCorrValues.Obs)



 
# Save LRPs for all OM scenarios
write.csv(SpwnCorr.df, paste(projOutDir2, "SpwnCorr.df.csv", sep="/"), row.names=F)
# 
# 
# dat<-as_tibble(SpwnCorr.df) %>% filter(OM_Name %in% c("Observed","Ricker_0.25cvER", "Ricker","Ricker_0.75cvER","Ricker_1.0cvER"))
# 
# 
# g <- ggplot(dat,aes(y=SpwnCorrValues,x=as.factor(OM_Name))) + geom_boxplot(width=0.5) +
#   scale_x_discrete(limits=c("Observed","Ricker_0.25cvER", "Ricker","Ricker_0.75cvER","Ricker_1.0cvER"),
#                    labels=c("Obs","0.25cv", "0.50cv(base)", "0.75cv", "1.0cv")) +
#                     xlab("Sensitivity Analysis Scenario") + ylab("Between-CU Correlation")
# 
# 
# ggsave(paste(cohoDir,"/Figures/coho-corrEffect_cvER",probThresh,".png",sep=""), plot = g,
#        width = 4.5, height = 3.5, units = "in")

 


dat<-as_tibble(SpwnCorr.df) %>% filter(OM_Name %in% c("Observed","Ricker","Ricker_0.25sigGamma", "Ricker_0.5sigGamma","Ricker_0.75sigGamma",
                                                      "Ricker_1.0sigGamma"))



g <- ggplot(dat,aes(y=SpwnCorrValues,x=as.factor(OM_Name))) + geom_boxplot(width=0.5) +
  scale_x_discrete(limits=c("Observed","Ricker", "Ricker_0.25sigGamma", "Ricker_0.5sigGamma","Ricker_0.75sigGamma",
                            "Ricker_1.0sigGamma"),
                   labels=c("Obs","0(base)", "0.0225sig", "0.045sig", "0.0625sig", "0.09sig")) +
                    xlab("Sensitivity Analysis Scenario") + ylab("Between-CU Correlation")


ggsave(paste(cohoDir,"/Figures/coho-corrEffect_sigGamma",probThresh,".png",sep=""), plot = g,
       width = 4.5, height = 3.5, units = "in")


 
# =============================================================================================
# (9) Look at posterior samples from MCMC model parameterization
# =============================================================================================
 
 Mod <- "SR_IndivRicker_Surv"
 OM<-"Ricker"
 
 #Mod <- "SR_IndivRicker_SurvCap"
 #OM<-"Ricker_priorCap"
 
 muLSurv<-SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
 muLSurv <- muLSurv$muLSurv
 
 post<-read.csv(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/SamSimInputs/",OM,"/",Mod,"_mcmc.csv", sep=""))
 adjProd<-exp(post$alpha + post$gamma * rep(muLSurv, length=nrow(post)))
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
 png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/Figures/ModelFits/",OM,"_AlphaPost.png", sep=""))
 do.call(grid.arrange, ps)
 dev.off()
 
 ps<-list()
 ps<-lapply(1:nCUs, plotPostHist, post=post, parName="sigma",  CUNames=CU_list)
 png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/Figures/ModelFits/",OM,"_SigmaPost.png", sep=""))
 do.call(grid.arrange, ps)
 dev.off()
 
 ps<-list()
 ps<-lapply(1:nCUs, plotPostHist, post=post, parName="adjProd", CUNames=CU_list)
 png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/Figures/ModelFits/",OM,"AdjProdPost.png", sep=""))
 do.call(grid.arrange, ps)
 dev.off()
 
 ps<-list()
 ps<-lapply(1:nCUs, plotPostHist, post=post, parName="beta", CUNames=CU_list)
 png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/Figures/ModelFits/",OM,"BetaPost.png", sep=""))
 do.call(grid.arrange, ps)
 dev.off()
 
 # Plotting code for gamma included directly here bc it differs from other pars bc only one value for all CUs
 # -- (could eventually be moved to plotFunctions.r)
 if(Mod %in% c("SR_HierRicker_SurvCap", "SR_IndivRicker_SurvCap","SR_HierRicker_Surv", "SR_IndivRicker_Surv")) {
   png(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/Figures/ModelFits/",OM,"GammaPost.png", sep=""))
   margPost<-post %>% filter(stk==1) %>% select("gamma")
   margPost<-data.frame(gamma = margPost$gamma)
   round.int<-2
   g<-ggplot(margPost, aes(x=gamma)) + geom_histogram() + labs(x="gamma", y = "Count") +
     geom_vline(aes(xintercept = mean(gamma)),col='red',size=1) +
     geom_text(aes(label=round(mean(gamma),round.int),y=0,x=mean(gamma)*1.25),
               vjust=-1,col='red',size=5) +
     geom_vline(aes(xintercept = quantile(gamma,0.95)[[1]]),col='blue',size=1) +
     geom_text(aes(label=round(quantile(gamma,0.95)[[1]],round.int),y=200,x=quantile(gamma,0.95)[[1]]*1.25),
               vjust=-1,col='blue',size=5) +
     geom_vline(aes(xintercept = quantile(gamma,0.05)[[1]]),col='blue',size=1) +
     geom_text(aes(label=round(quantile(gamma,0.05)[[1]],round.int),y=200,
                   x=quantile(gamma,0.05)[[1]]*1.25), vjust=-1,col='blue',size=5)
   print(g)
   dev.off()
 }
 

 
 # =============================================================================================
 # (9.1) Save posterior quantiles
 # =============================================================================================

 
 Mod <- "SR_IndivRicker_SurvCap"
 OM<-"Ricker_priorCap"
 
 # Mod <- "SR_IndivRicker_Surv"
 # OM<-"Ricker"
 
 muLSurv<-SRDat  %>% group_by(CU_ID) %>% summarise(muLSurv=mean(log(STAS_Age_3)))
 muLSurv <- muLSurv$muLSurv
 
 post<-read.csv(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/SamSimInputs/",OM,"/",Mod,"_mcmc.csv", sep=""))
 adjProd<-exp(post$alpha + post$gamma * rep(muLSurv, length=nrow(post)))
 post<- post %>% add_column(adjProd = adjProd)
 
 CU_list<-unique(SRDat[, "CU_Name"])
 CUID_list<-unique(SRDat[, "CU_ID"])
 nCUs<-length(CU_list)
 
 for (i in 1:5) {
 
 tmp<-post %>% filter(stk==i)
 Name<-unique(SRDat$CU_Name)[i]
 tmp2<-tmp %>% gather("variable", "value") %>% 
   group_by(variable) %>% 
   summarize(Mean = mean(value), 
             P05 = quantile(value, probs = .05),
             P50 = quantile(value, probs = .5),
             P75 = quantile(value, probs = .95)) %>% add_column(CU_Name = Name, .before="variable")
  if (i == 1) postSummary<-as.data.frame(tmp2)
  if ( i > 1) {
    postSummary<-rbind(postSummary,tmp2)
  }
 }
  
  write.csv(postSummary,file=paste(projOutDir,"/postSummary_",Mod,".csv",sep=""))
 
  
 
 
 
 
 # =============================================================================================
 # (9.2) Look at posterior samples from MCMC model parameterization, with a focus on those with low Sgen
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

 
 
 
 # ================================================================================================
 # (10) Make example plot to show calculation of projected LRPs
 # ====================================================================================================
 
 # Run example scenario:
 scenarioName <- "examplePlot"
 BMmodel <- "SR_IndivRicker_Surv"
 TMB_Inputs <- TMB_Inputs_IM
 
 run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
          useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
          nMCMC=30000, nProj=30000, cvER = 0.456*0.5, recCorScalar=1,gammaSigScalar=NULL,agePpnConst=TRUE,
          biasCorrectEst=T, biasCorrectProj=T)
 
 
 # Specify threshold to use when calculating LRP
 propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
 probThreshList<-c(0.50, 0.66, 0.90, 0.99) # probability threshhold; the LRP is set as the aggregate abundance that has this 
                                            # -- probability that the propCUThreshold is met
 
 
 for (p in 1:length(probThreshList)) {
   
   probThresh<-probThreshList[p]
     
     # Read in samSim outputs for OM
     filename<-paste("projLRPDat_","examplePlot",".csv",sep="")
     projLRPDat<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
     projLRPDat<-projLRPDat %>% filter(year > max(SRDat$yr_num)+4)
     
     # Create bins for projected spawner abundances
     binSize<-200 # Note: bin size is currently set here
     minBreak<-0
     maxBreak<-round(max(projLRPDat$sAg),digits=-2) + binSize
     breaks<-seq(minBreak, maxBreak,by=binSize)  
     
     # Set bin labels as the mid-point
     projLRPDat$bins<-cut(projLRPDat$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)),include.lowest = T)
     
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
     projLRPDat$OM.Name<-"examplePlot"
     projLRPDat.plot<-projLRPDat
     
     # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
     LRP<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))
     
     # Create a table of LRP estimates to be saved for each OM model
     if (p == 1) {
       LRP_Ests<-data.frame("examplePlot", probThresh, propCUThresh, LRP, binSize)
       names(LRP_Ests)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
     } else {
       tmp.df<-data.frame("examplePlot", probThresh, propCUThresh, LRP, binSize)
       names(tmp.df)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
       LRP_Ests<-rbind(LRP_Ests,tmp.df)
     }
  } # end of p loop     

      colList<-c("#E69F00", "#56B4E9", "#009E73", "#D55E00")
 
     # Plot projected LRP abundance relationship =============================================================== 
     png(paste(cohoDir,"/Figures/", "methods-Example-ProjectedLRP.png", sep=""), width=480, height=400)
     
     plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19, xlim=c(0,85000), ylim=c(0,1.0),cex=0.2,
          xlab="Aggregate Spawner Abundance", ylab="Pr (All CUs > Lower Benchmark)", cex.lab=1.5,
          cex.axis=1.2)
     
     abline(h=probThreshList[1],lty=3)
     lines(rep(LRP_Ests$LRP[LRP_Ests$ProbThresh==probThreshList[1]],2),c(0,probThreshList[1]), lty=1, lwd=2, col=colList[1])
     
     abline(h=probThreshList[2],lty=3)
     lines(rep(LRP_Ests$LRP[LRP_Ests$ProbThresh==probThreshList[2]],2),c(0,probThreshList[2]), lty=1, lwd=2, col=colList[2])
     
     abline(h=probThreshList[3],lty=3)
     lines(rep(LRP_Ests$LRP[LRP_Ests$ProbThresh==probThreshList[3]],2),c(0,probThreshList[3]), lty=1, lwd=2, col=colList[3])
     
     abline(h=probThreshList[4],lty=3)
     lines(rep(LRP_Ests$LRP[LRP_Ests$ProbThresh==probThreshList[4]],2),c(0,probThreshList[4]), lty=1, lwd=2, col=colList[4])
     
     dev.off()
     
     
     # Option to plot histogram of nSims in each Agg Abundance Bin
     #barplot(height = projLRPDat$nSims,names.arg = projLRPDat$bins)

  
# ===================================================================
# (11) Plot Correlation Matrices for Recruitment Residuals
# ==================================================================
     
 Mod<-"Ricker_priorCap"
 cormat<-read.csv(paste(cohoDir,"/SamSimInputs/",Mod,"/corrMat.csv", sep=""), header=F)
     
 CUNames <- unique(SRDat$CU_Name)
 CUNames[CUNames=="Middle_Fraser"]<-"Middle Fraser"
 CUNames[CUNames=="Fraser_Canyon"]<-"Fraser Canyon"
 CUNames[CUNames=="Lower_Thompson"]<-"Lower Thompson"
 CUNames[CUNames=="South_Thompson"]<-"South Thompson"
 CUNames[CUNames=="North_Thompson"]<-"North Thompson"
 
     
 # Plot correlation matrix to show scenarios
  # -- Method 1: ggplot
 library(ggcorrplot)
 rownames(cormat) <- CUNames
 colnames(cormat) <- CUNames
 g <-ggcorrplot(cormat,
        hc.order = TRUE,
                type = "lower",
                outline.color = "white",
                lab=T)
 
 # -- Method 2: corrplot
 library(corrplot)
 cormat<-as.matrix(cormat)
 rownames(cormat) <- CUNames
 colnames(cormat) <- CUNames

 png(filename=paste(cohoDir, "/Figures/coho-RecuitResidCorrelation-", Mod, ".png", sep=""), width=4, height=4.5, units="in", res=500)
  corrplot(cormat, method="circle", p.mat=cormat, insig="p-value", type="lower",addCoef.col="black", tl.col = "black")
 dev.off()


 
 
 
 
 
 
 # ===================================================================
 # (12) Make Projected Curve Comparison Plots Among Scenarios
 # ==================================================================
 
 # # Comparison 1: Among SR models ========================
 # OMsToPlot<-c("Ricker", "Ricker_priorCap", "Combined")
 # png(paste(cohoDir,"/Figures/", "coho-projLRPCurve-byOM.png", sep=""),
 #     width=450, height=450)
 # par(mfrow=c(2,2), mar=c(2,2,2,1), oma=c(2,2,0,0))
 # LabelsToPlot<-OMsToPlot
 # headerText<-"none"

 # #Comparison 2: Among sigGamma scenarios ========================
 # OMsToPlot<-c("Ricker", "Ricker_0.5sigGamma",
 #              "Ricker_0.75sigGamma", "Ricker_1.0sigGamma")
 # png(paste(cohoDir,"/Figures/", "coho-projLRPCurve-bySigGamma.png", sep=""),
 #     width=450, height=450)
 # par(mfrow=c(2,2), mar=c(2,2,2,1), oma=c(2,2,0,0))
 # LabelsToPlot<-c(0, 0.045, 0.0675, 0.09)
 # headerText<-"sigma"
 
 
 #Comparison 3: Among ER scenarios ========================
 OMsToPlot<-c("Ricker_ER2.5", "Ricker",
              "Ricker_ER22.5", "Ricker_ER32.5")
 png(paste(cohoDir,"/Figures/", "coho-projLRPCurve-byER.png", sep=""),
     width=450, height=450)
 par(mfrow=c(2,2), mar=c(2,2,2,1), oma=c(2,2,0,0))
 LabelsToPlot<-c(2.5, 12.5, 22.5, 32.5)
 headerText<-"ER"
 
 
 
 # Specify threshold to use when calculating LRP
 propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
 probThreshList<-c(0.50, 0.66, 0.90, 0.99) # probability threshhold; the LRP is set as the aggregate abundance that has this 
 # -- probability that the propCUThreshold is met

 for (i in 1:length(OMsToPlot)) {
 
   
   if (OMsToPlot[i]!="Combined") {
     # Read in samSim outputs for OM
     filename<-paste("projLRPDat_",OMsToPlot[i],".csv",sep="")
     projLRPDat<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
     projLRPDat<-projLRPDat %>% filter(year > max(SRDat$yr_num)+4)
   } else {
     filename1<-"projLRPDat_Ricker.csv"
     filename2<-"projLRPDat_Ricker_priorCap.csv"
     projLRPDat1<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename1))
     projLRPDat2<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename2))
     projLRPDat<-rbind(projLRPDat1,projLRPDat2)
     projLRPDat<-projLRPDat %>% filter(year > max(SRDat$yr_num)+4)
   }
   
   # Create bins for projected spawner abundances
   binSize<-200 # Note: bin size is currently set here
   minBreak<-0
   maxBreak<-round(max(projLRPDat$sAg),digits=-2) + binSize
   breaks<-seq(minBreak, maxBreak,by=binSize)  
   
   # Set bin labels as the mid-point
   projLRPDat$bins<-cut(projLRPDat$sAg,breaks=breaks,labels=as.character(rollmean(breaks,k=2)),include.lowest = T)
   
   # Summarize nSims in each bin
   tmp<-projLRPDat %>% group_by(bins) %>% summarise(nSims=(length(ppnCUsLowerBM)))
   
   # Filter out bins with < 100 nSims
   tmp2<-projLRPDat %>% group_by(bins) %>% summarise(nSimsProp1=(length(ppnCUsLowerBM[ppnCUsLowerBM == propCUThresh]))) %>%
     add_column(nSims=tmp$nSims) %>% filter(nSims>=100)
   
 
 for (p in 1:length(probThreshList)) {
   
   probThresh<-probThreshList[p]
   
   # For each bin, calculate probability that required proportion of CUs above benchmark
   projLRPDat<-tmp2 %>% add_column(prob=tmp2$nSimsProp1/tmp2$nSims)
   # For each bin, calculate the difference between the threshold probability and the calculated probability 
   projLRPDat$diff<-abs(probThresh-projLRPDat$prob)
   
   # Save projection summaries used to create plots
   projLRPDat$OM.Name<-OMsToPlot[i]
   projLRPDat.plot<-projLRPDat
   
   # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
   if (min(abs(projLRPDat$diff)) < 0.02) {
      LRP<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))
   } else {
     LRP<-NA
    }
   
   # Create a table of LRP estimates to be saved for each OM model
   if (p == 1) {
     LRP_Ests<-data.frame(OMsToPlot[i], probThresh, propCUThresh, LRP, binSize)
     names(LRP_Ests)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
   } else {
     tmp.df<-data.frame(OMsToPlot[i], probThresh, propCUThresh, LRP, binSize)
     names(tmp.df)<-c("OM", "ProbThresh", "PropCURequired", "LRP", "binSize")
     LRP_Ests<-rbind(LRP_Ests,tmp.df)
   }
 } # end of p loop     
 
 colList<-c("#E69F00", "#56B4E9", "#009E73", "#D55E00")
 
 plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19, xlim=c(0,90000), ylim=c(0,1.0),cex=0.2,
      xlab="", ylab="", cex.axis=1.1)

 abline(h=probThreshList[1],lty=3)
 lines(rep(LRP_Ests$LRP[LRP_Ests$ProbThresh==probThreshList[1]],2),c(0,probThreshList[1]), lty=1, lwd=2, col=colList[1])
 
 abline(h=probThreshList[2],lty=3)
 lines(rep(LRP_Ests$LRP[LRP_Ests$ProbThresh==probThreshList[2]],2),c(0,probThreshList[2]), lty=1, lwd=2, col=colList[2])
 
 abline(h=probThreshList[3],lty=3)
 lines(rep(min(LRP_Ests$LRP[LRP_Ests$ProbThresh==probThreshList[3]]),2),c(0,probThreshList[3]), lty=1, lwd=2, col=colList[3])
 
 abline(h=probThreshList[4],lty=3)
 lines(rep(min(LRP_Ests$LRP[LRP_Ests$ProbThresh==probThreshList[4]]),2),c(0,probThreshList[4]), lty=1, lwd=2, col=colList[4])
 
 if (headerText == "none") mtext(LabelsToPlot[i], side = 3, outer = FALSE, line=0.5)
  
 if (headerText == "sigma") {
   sigVal<-LabelsToPlot[i]
   mtext(paste ("sigGamma =",sigVal), side = 3, outer = FALSE, line=0.5)
 }
 
 if (headerText == "ER") {
   mtext(paste ("ER =",LabelsToPlot[i]), side = 3, outer = FALSE, line=0.5)
 }
 
 }  # end of OM loop
 
 mtext("Aggregate Spawner Abundances", side = 1, outer = TRUE, line=0.5, cex=1.2, font=1)
 mtext("Pr (All CUs > Lower Benchmark)", side = 2, outer = TRUE, line=0.5, cex=1.2, font=1)
 
 
 dev.off()
 
 
 
 
 
 
 
 # ===================================================================
 # (13) Make Comparison Plots Among Scenarios (NOT CURRENTLY WORKING)
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
