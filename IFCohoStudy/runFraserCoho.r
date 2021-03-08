# Code by K.Holt & B.Davis
# Started on: February 7, 2020

# This file contains the code required to explore Interior Fraser Coho escapement datasets and run retrospective analyses of 
#   data-based LRPs that use logistic regressions

#library(MASS) # dose.p function to get SE around P95
library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(tmbstan)
library(here)

setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
cohoDir<-paste(rootDir,"/IFCohoStudy",sep="")

setwd(codeDir)

sourceAll <- function(){
  source("benchmarkFunctions.r")
  source("LRPFunctions.r")
  source("plotFunctions.r")
  source("retroFunctions.r")
  source("helperFunctions.r")
  source("LRdiagnostics.r")
}

sourceAll()

# Load TMB models

compile("TMB_Files/LRP_Logistic_Only.cpp")
dyn.load(dynlib("TMB_Files/LRP_Logistic_Only"))

compile("TMB_Files/SR_HierRicker_Surv_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv_LowAggPrior"))

compile("TMB_Files/SR_IndivRicker_Surv_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_LowAggPrior"))

compile("TMB_Files/SR_HierRicker_SurvCap_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap_LowAggPrior"))

compile("TMB_Files/SR_IndivRicker_SurvCap_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap_LowAggPrior"))

compile("TMB_Files/SR_IndivRicker_Surv_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_noLRP"))

#compile("TMB_Files/Proj_LRP.cpp")
#dyn.load(dynlib("TMB_Files/Proj_LRP"))

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

CoEscpDat_bySubpop<-read.csv("DataIn/IFCoho_escpBySubpop.csv")

# Change column names to yr, CU_Name, Escp, Subpop_Name, CU_ID
  CoEscpDat_bySubpop<-CoEscpDat_bySubpop %>% select(yr=Return.Year, CU_Name=Conservation.Unit, Escp=Natural.Returns, Subpop_Name=Sub.Population)
  tmp.df<-data.frame(CU_Name=unique(CoEscpDat_bySubpop$CU_Name), CU_ID=seq(1,length(unique(CoEscpDat_bySubpop$CU_Name)),by=1))
  CoEscpDat_bySubpop <- left_join(CoEscpDat_bySubpop,tmp.df)

setwd(codeDir)


# Run annual restrospective analyses over various levels of p, restrict escapement dataset to 1998+ ============================================================
# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)
CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)
CoEscpDat_bySubpop <- CoEscpDat_bySubpop %>% filter(yr >= 1998)

# Roll up escpaments, and get Gen Mean of htat
AggEscp <- CoEscpDat %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Escp, 3, gm_mean, fill = NA, align="right"))


# ==================================================================================
# Call functions to plot data availability:
# ====================================================================================
# plot_CU_DataObs_Over_Time(CoEscpDat, cohoDir, plotName="Fr_Co_DataByCU")
# plot_Num_CUs_Over_Time(CoEscpDat, cohoDir, plotName="Fr_Co_N_CUs")
# 
# # Note: these next 2 two escpt plots need to have formatting fixed
# plot_CU_Escp_Over_Time(CoEscpDat, cohoDir, plotName="IFC Esc", samePlot = T)
# plot_CU_Escp_Over_Time(CoEscpDat, cohoDir, plotName="IFC Esc Separate", samePlot = F)
# plot_Subpop_Escp_Over_Time(CoEscpDat_bySubpop, cohoDir, plotName="IFC Esc Separate - by Subpop", samePlot = F)
# plot_Subpop_Escp_Over_Time(CoEscpDat_bySubpop, cohoDir, plotName="IFC Esc - by Subpop", samePlot = T)


# ==================================================================================================================
# Run retrospective analyses:
# =====================================================================================================================

B_penalty_mu<-4489
B_penalty_sigma<-2048

# TMB input parameters:
TMB_Inputs_HM <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1, 
                B_penalty_mu=B_penalty_mu, B_penalty_sigma=B_penalty_sigma)


TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                      B_penalty_mu=B_penalty_mu, B_penalty_sigma=B_penalty_sigma)


# Prior means come from running "compareRickerModelTypes.r"
cap_priorMean_HM<-c(10.957092, 5.565526, 11.467815, 21.104274, 14.803877)

TMB_Inputs_HM_priorCap <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                   logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                   gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                   cap_mean=cap_priorMean_HM, cap_sig=sqrt(2),B_penalty_mu=B_penalty_mu, B_penalty_sigma=B_penalty_sigma)

# Prior means come from running "compareRickerModelTypes.r"
cap_priorMean_IM<-c(11.153583,  5.714955, 11.535779, 21.379558, 14.889006)

TMB_Inputs_IM_priorCap <- list(Scale = 1000, logA_Start = 1, Tau_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                               cap_mean=cap_priorMean_IM, cap_sig=sqrt(2),B_penalty_mu=B_penalty_mu, B_penalty_sigma=B_penalty_sigma)


TMB_Inputs_Subpop <- list(Scale = 1000)

# TMB_Inputs_Subpop <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
#                           logMuA_sig = 2, Tau_dist = 0.1, Tau_A_dist = 0.1, 
#                           gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)


# # Test projected LRP
# devtools::install_github("Pacific-salmon-assess/samSim", ref="LRP")
# ps <- c(0.6, 0.8, 1.0)
# 
# for(pp in 1:length(ps)){
# 
#   runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
#                BMmodel = "SR_IndivRicker_Surv", LRPmodel="Proj", integratedModel=F,
#                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Proj.IndivRickerSurv_",ps[pp]*100, sep=""),
#                bootstrapMode = F, plotLRP=T)
# }



# Run annual restrospective analyses using CUs ===========================

# Note: if .cpp files are already compiled (.dll and .o files are present in SalmonLRP_RetroEval/Code/TMB_Files ), 
# code below will run and produce figures, but will give this error message: 
# "Error in .Call("FreeADFunObject", ptr, PACKAGE = DLL) : "FreeADFunObject" not available for .Call() for package "<TMB model name>""
# Deleting the .dll and .o files in /TMB_Files folder and recompiling .cpp files gets rid of this error. 

# Note: BroodYrLag = the number of years to subract from the current year to get to the most recent brood year

# Loop over p values and run annual retrospective analyses for each level of p
  ps <- c(0.8, 0.99)
  for(pp in 1:length(ps)){
    # Run with Binomial LRP model with hierarchical Ricker 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_HierRicker_Surv_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bin.HierRickerSurv_",ps[pp]*100, sep=""),
                  bootstrapMode = F, plotLRP=T, runLogisticDiag=T)
    
    # # Run with Bernoulli LRP model with hierarchical Ricker 
    # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
    #                BMmodel = "SR_HierRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
    #                useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bern.HierRickerSurv_",ps[pp]*100, sep=""),
    #                bootstrapMode = F, plotLRP=T)
    
    # Run with Binomial LRP model with individual model Ricker 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_IndivRicker_Surv_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bin.IndivRickerSurv_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T,runLogisticDiag=T)
    
    # # Run with Bernoulli LRP model with individual model Ricker 
    # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
    #                BMmodel = "SR_IndivRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
    #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurv_",ps[pp]*100, sep=""),
    #                bootstrapMode = F, plotLRP=T)
    
    # Run with Binomial LRP model with hierarchical Ricker, with prior on capacity 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_HierRicker_SurvCap_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_HM_priorCap, outDir=cohoDir, RunName = paste("Bin.HierRickerSurvCap_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T,runLogisticDiag=T)
  
    
    # # Run with Bernoulli LRP model with hierarchical Ricker, with prior on capacity 
    # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
    #                BMmodel = "SR_HierRicker_SurvCap", LRPmodel="BernLogistic", integratedModel=T,
    #                useGenMean=F, TMB_Inputs=TMB_Inputs_HM_priorCap, outDir=cohoDir, RunName = paste("Bern.HierRickerSurvCap_",ps[pp]*100, sep=""),
    #                bootstrapMode = F, plotLRP=T)
    
    # Run with Binomial LRP model with individual model Ricker, with prior on capacity 
     runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                    BMmodel = "SR_IndivRicker_SurvCap_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
                    useGenMean=F, TMB_Inputs=TMB_Inputs_IM_priorCap, outDir=cohoDir, RunName = paste("Bin.IndivRickerSurvCap_",ps[pp]*100, sep=""),
                    bootstrapMode = F, plotLRP=T,runLogisticDiag=T)
    
  
    # # Run with Bernoulli LRP model with individual model Ricker, with prior on capacity 
    # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
    #                BMmodel = "SR_IndivRicker_SurvCap", LRPmodel="BernLogistic", integratedModel=T,
    #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM_priorCap, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurvCap_",ps[pp]*100, sep=""),
    #                bootstrapMode = F, plotLRP=T)
    
    
  }



# Run annual restrospective analyses using subpopulations ===========================

ps <- c(seq(0.6, 0.95,.05), 0.99)
for(pp in 1:length(ps)){

runAnnualRetro(EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
               BMmodel = "ThreshAbund_Subpop1000_ST", LRPmodel="BernLogistic", integratedModel=F,
               useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bern.SPopAbundThreshST_",ps[pp]*100, sep=""),
               bootstrapMode = F, plotLRP=T)

runAnnualRetro(EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
               BMmodel = "ThreshAbund_Subpop1000_LT", LRPmodel="BernLogistic", integratedModel=F,
               useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bern.SPopAbundThreshLT_",ps[pp]*100, sep=""),
               bootstrapMode = F, plotLRP=T)

runAnnualRetro(EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
               BMmodel = "ThreshAbund_Subpop1000_ST", LRPmodel="BinLogistic", integratedModel=F,
               useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bin.SPopAbundThreshST_",ps[pp]*100, sep=""),
               bootstrapMode = F, plotLRP=T)

runAnnualRetro(EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
               BMmodel = "ThreshAbund_Subpop1000_LT", LRPmodel="BinLogistic", integratedModel=F,
               useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bin.SPopAbundThreshLT_",ps[pp]*100, sep=""),
               bootstrapMode = F, plotLRP=T)

}

# Run all available combinations of number of CUs ============================================================

# For CU approach ===========================================================

nCUs<-c(5,3,4)
ps<-c(0.60, 0.80,0.99)


for (pp in 1:length(ps)){

  runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2,
               genYrs=3, p = ps[pp], BMmodel="SR_IndivRicker_Surv", LRPmodel="BinLogistic", integratedModel=T,
               useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bin.IndivRickerSurv_",ps[pp]*100, sep=""))

  runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, 
               genYrs=3, p = ps[pp], BMmodel="SR_IndivRicker_SurvCap", LRPmodel="BinLogistic", integratedModel=T,
               useGenMean=F, TMB_Inputs=TMB_Inputs_IM_priorCap, outDir=cohoDir, RunName = paste("Bin.IndivRickerSurvCap_",ps[pp]*100, sep="")) 
  
  runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2,
             genYrs=3, p = ps[pp], BMmodel="SR_HierRicker_Surv", LRPmodel="BinLogistic", integratedModel=T,
             useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bin.HierRickerSurv_",ps[pp]*100, sep=""))

  runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, 
               genYrs=3, p = ps[pp], BMmodel="SR_HierRicker_SurvCap", LRPmodel="BinLogistic", integratedModel=T,
               useGenMean=F, TMB_Inputs=TMB_Inputs_HM_priorCap, outDir=cohoDir, RunName = paste("Bin.HierRickerSurvCap_",ps[pp]*100, sep="")) 
  
  runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2015, endYr=2018, BroodYrLag=2,
               genYrs=3, p = ps[pp], BMmodel="ThreshAbund_Subpop1000_ST", LRPmodel="BinLogistic", integratedModel=F,
               useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bin.SPopAbundThreshST_",ps[pp]*100, sep=""))

  
}




# =======================================================================
# Make Plots for Retrospective Analysis: 
# =======================================================================


# Plot aggregate status by method and number of CUs

yearList<-2015:2018

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "SR_HierRicker_Surv",p=0.80, Dir=cohoDir,
                     inputPrefix="Bin.HierRickerSurv_80",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "SR_IndivRicker_Surv",p=0.80, Dir=cohoDir,
                     inputPrefix="Bin.IndivRickerSurv_80",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "SR_HierRicker_SurvCap",p=0.80, Dir=cohoDir,
                      inputPrefix="Bin.HierRickerSurvCap_80",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "SR_IndivRicker_SurvCap",p=0.80, Dir=cohoDir,
                      inputPrefix="Bin.IndivRickerSurvCap_80",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "ThreshAbund_Subpop1000_ST", p=0.95, Dir=cohoDir,
                     inputPrefix="Bin.SPopAbundThreshST_80", plotAveLine=TRUE)






# 
# 
# # Plot CVs on LRP estimates by method and number of CUs
# yearList<-2015:2018
# 
# plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_HierRicker_Surv", p=0.80, Dir=cohoDir,
#                    inputPrefix = "retroStatus_CU.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))
# 
# plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_HierRicker_Surv", p=0.90, Dir=cohoDir,
#                    inputPrefix = "retroStatus_CU.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))
# 
# plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_HierRicker_Surv", p=0.95, Dir=cohoDir,
#                    inputPrefix = "retroStatus_CU.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))
# 
# plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "SR_HierRicker_Surv", p=0.80, Dir=cohoDir,
#                    inputPrefix = "retroStatus_CU.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))
# 
# plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "SR_HierRicker_Surv", p=0.90, Dir=cohoDir,
#                    inputPrefix = "retroStatus_CU.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))
# 
# plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "ThreshAbund_Subpop1000_ST", p=0.95, Dir=cohoDir,
#                    inputPrefix = "retroStatus_SPop.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))


# Plot annual retrospective status with bars for various options for LRP metric 

# Make a list of all modelled combinations to be shown in plots

modelFitList<-c("Bin.IndivRickerSurv_60",
                "Bin.HierRickerSurv_60",
                "Bin.IndivRickerSurvCap_60",
                "Bin.HierRickerSurvCap_60",
                "Bin.SPopAbundThreshST_60",
                "Bin.IndivRickerSurv_80",
                "Bin.HierRickerSurv_80",
                "Bin.IndivRickerSurvCap_80",
                "Bin.HierRickerSurvCap_80",
                "Bin.SPopAbundThreshST_80",
                "Bin.IndivRickerSurv_99",
                "Bin.HierRickerSurv_99",
                "Bin.IndivRickerSurvCap_99",
                "Bin.HierRickerSurvCap_99",
                "Bin.SPopAbundThreshST_99")

# Specify which thresholds should be tested in the data-based proportional thresholds
ps_Prop<-c(0.6, 0.8, 1.0)

# Provide WSP assessment results
WSP_estYr<-2014
WSP_AboveLRP<-TRUE

LRP_estYr<-2015
retroYears<-2000:LRP_estYr
plotStatus_byYear(LRP_estYr, retroYears, Dir=cohoDir, genYrs=3, AggEscp, EscpDat=CoEscpDat, 
                  modelFitList=modelFitList, ps_Prop=ps_Prop, WSP_estYr, WSP_AboveLRP,
                  outDir = cohoDir, fName = paste("statusPlot_withBars",LRP_estYr,sep=""))

LRP_estYr<-2016
retroYears<-2000:LRP_estYr
plotStatus_byYear(LRP_estYr, retroYears, Dir=cohoDir, genYrs=3, AggEscp, EscpDat=CoEscpDat, 
                  modelFitList=modelFitList, ps_Prop=ps_Prop, WSP_estYr, WSP_AboveLRP,
                  outDir = cohoDir, fName = paste("statusPlot_withBars",LRP_estYr,sep=""))

LRP_estYr<-2017
retroYears<-2000:LRP_estYr
plotStatus_byYear(LRP_estYr, retroYears, Dir=cohoDir, genYrs=3, AggEscp, EscpDat=CoEscpDat, 
                  modelFitList=modelFitList, ps_Prop=ps_Prop, WSP_estYr, WSP_AboveLRP,
                  outDir = cohoDir, fName = paste("statusPlot_withBars",LRP_estYr,sep=""))

LRP_estYr<-2018
retroYears<-2000:LRP_estYr
plotStatus_byYear(LRP_estYr, retroYears, Dir=cohoDir, genYrs=3, AggEscp, EscpDat=CoEscpDat, 
                  modelFitList=modelFitList, ps_Prop=ps_Prop, WSP_estYr, WSP_AboveLRP,
                  outDir = cohoDir, fName = paste("statusPlot_withBars",LRP_estYr,sep=""))


# ============================================================================================
# Create plot to compare LRP estimates among methods for a single year
# ===========================================================================================

modelFitList<-c("Bin.IndivRickerSurv",
                "Bin.HierRickerSurv",
                "Bin.IndivRickerSurvCap",
                "Bin.HierRickerSurvCap",
                "Bin.SPopAbundThreshST")

pList<-c(60, 80, 99)

LRP_estYr<-2018

plotLRPCompare(LRP_estYr, modelFitList, pList, outDir = cohoDir, 
               fName = paste("compareLRPs_withBars",LRP_estYr,sep=""))

# ============================================================================================
# Create plot to compare retrospective time series of LRP estimates
# ===========================================================================================

EscpDat.yy <- CoEscpDat %>% filter(yr <= 2018) 

modelFitList<-c("Bin.IndivRickerSurv",
                "Bin.HierRickerSurv",
                "Bin.IndivRickerSurvCap",
                "Bin.HierRickerSurvCap",
                "Bin.SPopAbundThreshST")

L_Names<-c("Sgen: IM", "Sgen: HM","Sgen: IM.HiSrep","Sgen: HM.HiSrep", "Distributional")  

pList<-c(60, 80, 99)


plotAnnualRetro_Compare(Dat=EscpDat.yy, Names = modelFitList, pList=pList, L_Names = L_Names, outDir=cohoDir, useGenMean = T, genYrs = 3)
  
# ============================================================================================
# Create plot to compare status as a function of nCUs for different LRP options
# ===========================================================================================
  

modelFitList<-c("Bin.IndivRickerSurv",
                "Bin.HierRickerSurv",
                "Bin.IndivRickerSurvCap",
                "Bin.HierRickerSurvCap",
                "Bin.SPopAbundThreshST")


plotAggStatus_byNCUs_Compare(estYear=2018, nCUList=c(5,4,3), Names=modelFitList, p=0.80, Dir=cohoDir,plotAveLine=TRUE)  


# Notes on comparing LW results run (2020-10-26) with results on dropbox https://www.dropbox.com/sh/otiku88jc2eu8cx/AACagLQd85blX7jc8-yBYs4Na?dl=0 
#   1. DataOut/ModelFits/AllEsts_Hier_Ricker_Surv.csv: 
#         sigmaA estimate is 0.321763 vs. -1.133939234 on dropbox
#         other values look the same
#   2. DataOut/AnnualRetrospective/annualRetro_LRPs.csv:
#         values are identical



# ========================================================================
# Run leave-one-out logistic model diagnostic
# ========================================================================


year <- 2018
p<-0.99
useGenMean <- FALSE
genYrs<-3
BroodYrLag<-2
useBern_Logistic <- FALSE


#CoEscpDat
#CoSRDat




#LOO_LRdiagnostics <- function(remove.EnhStocks=TRUE, n=18){
  
  
  # Step 1: estimate logistic regression iteratively, removing a single year 
  # each time
  
  # predPpnAboveBM <- NA
  # 
  # for (i in 1:n){
  #   # Estimate logistic regression using function Get.LRP
  #   zz <- Get.LRP(remove.EnhStocks = TRUE, LOO=i)
  #   All_Ests <- zz$out$All_Ests
  #   
  #   if(i==1){ # These remain constant over iterations
  #     # Step 2: Get observed time-series of aggregate raw abundances that includes all
  #     # data and then scale to units near 1-10 
  #     AggAbundRaw <- zz$out$Logistic_Data$xx
  #     digits <- count.dig(AggAbundRaw)
  #     ScaleSMU <- min(10^(digits -1 ), na.rm=T)
  #     AggAbund <- AggAbundRaw/ScaleSMU
  #     # Get time-series of observed ppns of CUs> benchamark, including all
  #     # data
  #     obsPpnAboveBM <- zz$out$Logistic_Data$yy
  #     # Get threshold p value (ppn of CUs>benchmark) used to estimate LRP
  #     p <- zz$LRPppn
  #     #dir <- "DataOut/"
  #   }
  #   
  #   # Step 3: Get predicted ppn of CUs above their lower benchmark for the year 
  #   # that was held out
  #   B_0 <- All_Ests %>% filter(Param=="B_0") %>% pull(Estimate)
  #   B_1 <- All_Ests %>% filter(Param=="B_1") %>% pull(Estimate)
  #   #predPpnAboveBM <- inv_logit(B_0 + B_1*AggAbund)
  #   predPpnAboveBM[i] <- inv_logit(B_0 + B_1*AggAbund[i])
  # } # End of for i in 1:18
  # 
  # # Step 4: Calculate Hit Ratio
  # 
  # # In which years did the model predict aggregate abundances >LRP?
  # yHat <- predPpnAboveBM > p
  # # In which years were observed aggregate abundances >LRP?
  # y <- obsPpnAboveBM > p
  # 
  # # Confusion Matrix
  # confMat <- table(y, yHat)
  # 
  # # What is the accuracy in classifying observed aggregate abundances?
  # # Hit ratio = ratio of correct classification
  # hitRatio <- sum(diag(confMat))/sum(confMat)
  # hitRatio <- round(hitRatio, digits=2)
  # 
  # 
  # 
  # #}# End of Function 2: LOO_LRdiagnostics()
  # 
