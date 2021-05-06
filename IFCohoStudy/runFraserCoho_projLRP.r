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
 
 
# TMB input parameters:
TMB_Inputs_HM <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                      logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)


TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)


# Prior means come from running "compareRickerModelTypes.r"
cap_priorMean_HM<-c(10.957092, 5.565526, 11.467815, 21.104274, 14.803877)

TMB_Inputs_HM_priorCap <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                               logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                               cap_mean=cap_priorMean_HM, cap_sig=sqrt(2))

# Prior means come from running "compareRickerModelTypes.r"
cap_priorMean_IM<-c(11.153583,  5.714955, 11.535779, 21.379558, 14.889006)

TMB_Inputs_IM_priorCap <- list(Scale = 1000, logA_Start = 1, Tau_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                               cap_mean=cap_priorMean_IM, cap_sig=sqrt(2))



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


# Create samSim input files for current scenario
scenarioName <- "IM.base"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456, recCorScalar=1,gammaSigScalar=NULL)


# Create samSim input files for current scenario
scenarioName <- "IM.highGammaSig"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456, recCorScalar=1,gammaSigScalar=1.5)

# Create samSim input files for current scenario
scenarioName <- "IM.medGammaSig"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=5000, nProj=2000, cvER = 0.456, recCorScalar=1,gammaSigScalar=0.5)



# Create samSim input files for current scenario
scenarioName <- "IM.base.test"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM
projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=500, nProj=10, cvER = 0.456, recCorScalar=1,gammaSigScalar=NULL)


scenarioName <- "HM.base"
BMmodel <- "SR_HierRicker_Surv"
TMB_Inputs <- TMB_Inputs_HM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=5000, nProj=2000, cvER = 0.456, recCorScalar=1)


scenarioName <- "IMCap.base"
BMmodel <- "SR_IndivRicker_SurvCap"
TMB_Inputs <- TMB_Inputs_IM_priorCap

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=5000, nProj=2000, cvER = 0.456, recCorScalar=1)


scenarioName <- "HMCap.base"

BMmodel <- "SR_HierRicker_SurvCap"
TMB_Inputs <- TMB_Inputs_HM_priorCap

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=T,
                                nMCMC=5000, nProj=2000, cvER = 0.456, recCorScalar=1)




# ==================================================================
# (4) Run Sensitivity Analyses 
# ====================================================================


# Create samSim input files for current scenario
scenarioName <- "IM.cvER1.5"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=NA, nProj=5000, cvER = 0.456*1.5, recCorScalar=1)


# Create samSim input files for current scenario
scenarioName <- "IM.cvER2.0"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=NA, nProj=5000, cvER = 0.456*2, recCorScalar=1)


scenarioName <- "IM.cvER3.0"
BMmodel <- "SR_IndivRicker_Surv"
TMB_Inputs <- TMB_Inputs_IM

projSpawners <-run_ScenarioProj(SRDat = SRDat, BMmodel = BMmodel, scenarioName=scenarioName,
                                useGenMean = F, genYrs = genYrs,  TMB_Inputs, outDir=cohoDir, runMCMC=F,
                                nMCMC=NA, nProj=5000, cvER = 0.456*3.0, recCorScalar=1)




# ===================================================================
#  (5) Estimate and Save LRPs (and associated plots)
# ==================================================================

# Specify threshold to use when calculating LRP
  # # Note: may want to loop over probThresholds as well; still needs to be added
propCUThresh <- 1.0 # required proportion of CUs above lower benchmark
probThresh<-0.50 # probability theshhold; the LRP is set as the aggregate abundance that has this 
                    # probability that the propCUThreshold is met

# Specify scenarios to calculate LRPs and make plots for.
# These scenarios will be looped over below with a LRP (and LRP plot) saved for each scenario

OMsToInclude<-c("IM.base","IM.highGammaSig", "IM.modGammaSig") 


# Loop over OM Scenarios 
for (i in 1:length(OMsToInclude)) {
 
  # Read in samSim outputs for OM
  filename<-paste("projLRPDat_",OMsToInclude[i],".csv",sep="")
  projLRPDat<-read.csv(here(cohoDir, "SamSimOutputs", "simData",filename))
  projLRPDat<-projLRPDat %>% filter(year > max(SRDat$yr_num)+4)

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
  projLRPDat$OM.Name<-OMsToInclude[i]
  if (i == 1) projLRPDat.plot<-projLRPDat
  if (i > 1) projLRPDat.plot<-rbind(projLRPDat.plot,projLRPDat)
  
  # Calculate the LRP as aggregate abundance bin with the minimum difference from threshold
  LRP<-as.numeric(as.character(projLRPDat$bins[projLRPDat$diff == min(projLRPDat$diff)]))
  
  # Create a table of LRP estimates to be saved for each OM model
  if (i ==1) {
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
  
  plot(as.numeric(as.character(projLRPDat$bins)),projLRPDat$prob, pch=19, xlim=c(0,85000), cex=0.2,
       xlab="Aggregate Abundance", ylab="Pr (All CUs > Lower Benchmark)")
  abline(h=probThresh, lty=2)
  abline(v=LRP, col="orange", lwd=2)
  
  dev.off()
  
  # Option to plot histogram of nSims in each Agg Abundance Bin
  #barplot(height = projLRPDat$nSims,names.arg = projLRPDat$bins)
  
}


# Save LRPs for all OM scenarios
write.csv(LRP_Ests, paste(projOutDir2, "ProjectedLRPs.csv", sep="/"), row.names=F)
# Save LRP projection summaries used for calculating and plotting LRP (Optional)
write.csv(projLRPDat.plot, paste(projOutDir2, "ProjectedLRP_data.csv", sep="/"), row.names=F)






# ===================================================================
# (6) Plot CU-level Spawner Abundance Projections (Optional)
# ==================================================================


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
# (7) Make Comparison Plots Among Scenarios (NOT CURRENTLY WORKING)
# ==================================================================

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
# # Plot to show sensitivity analysis to recruitment correlation
# 
# #OMsToPlot<-c("IM.Base","IM.80RecCor","IM.60RecCor","IM.40RecCor","IM.10RecCor")
# 
# OMsToPlot<-c("IM.Base","IM.60RecCor","IM.40RecCor","IM.10RecCor", "IM-.20RecCor","IM-.40RecCor")
# 
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
# # ---- Temporary code to test effect of correlation scalar on projected correlations
# 
# MPD_recCormat<-read.csv(paste(cohoDir,"/SamSimInputs/IM.base/cohoCorrMat.csv", sep=""), header=F)
# 
# OMsToTest<-c("IM.Base","IM.80RecCor","IM.60RecCor","IM.40RecCor","IM.40RecCor","IM-.40RecCor")
# 
# #OMsToTest<-c("IM.cvER1.25", "IM.cvER1.5","IM.cvER2.0","IM.cvER2.5","IM.cvER3.0")
# 
# filename<-paste( "projSpwnDat_",OMsToTest[6],".csv",sep="")
# spDat<-read.csv(here(cohoDir,"SamSimOutputs", "simData",filename))
# 
# # Calculate correlation matrix in MPD recruitment residuals ========================
# #resids<-as_tibble(data.frame(stock=data$stk,year=data$yr, resid=obj$report()$R_Resid))
# spDat<-as_tibble(spDat)
# spDat<-spDat%>%select(-X)
# 
# RecCorMat<-array(NA,dim=c(5,5,max(spDat$iteration)))
# SpwnCorMat<-array(NA,dim=c(5,5,max(spDat$iteration)))
# 
# for (i in 1:max(spDat$iteration)) {
#   
#   recruits.i<-spDat %>% filter(iteration==i & expRate==0.125) %>%  select(-spawners) 
#   recruits.i<-recruits.i %>% select(-expRate, -iteration) 
#   cor_mat<-recruits.i %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=recruits) %>% select(-year) %>% cor()
#   RecCorMat[,,i]<-cor_mat
# 
#   spawners.i<-spDat %>% filter(iteration==i & expRate==0.125) %>%  select(-recruits) 
#   spawners.i<-spawners.i %>% select(-expRate, -iteration) 
#   cor_mat<-spawners.i %>% pivot_wider(names_from = CU, names_prefix="CU", values_from=spawners) %>% select(-year) %>% cor()
#   SpwnCorMat[,,i]<-cor_mat
# 
# }
# 
# 
# RecCorMat_Ave<-apply(RecCorMat, c(1,2), mean)
# SpwnCorMat_Ave<-apply(SpwnCorMat, c(1,2), mean)
# 
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
