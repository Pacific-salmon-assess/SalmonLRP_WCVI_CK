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

compile("TMB_Files/SR_HierRicker_Surv_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv_LowAggPrior"))

compile("TMB_Files/SR_IndivRicker_Surv_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_LowAggPrior"))

compile("TMB_Files/SR_HierRicker_SurvCap_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap_LowAggPrior"))

compile("TMB_Files/SR_IndivRicker_SurvCap_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap_LowAggPrior"))

compile("TMB_Files/SR_HierRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv"))

compile("TMB_Files/SR_IndivRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv"))

compile("TMB_Files/SR_HierRicker_SurvCap.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap"))

compile("TMB_Files/SR_IndivRicker_SurvCap.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap"))

# compile("TMB_Files/LRP_Logistic_Only.cpp")
# dyn.load(dynlib("TMB_Files/LRP_Logistic_Only"))

# Only needed if running projection models =======
#compile("TMB_Files/SR_IndivRicker_Surv_noLRP.cpp")
#dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_noLRP"))

# Only needed for leave-one-out logistic regression diagnostics
compile("TMB_Files/LRP_BasicLogistic_Only_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/LRP_BasicLogistic_Only_LowAggPrior"))

# Only needed for leave-one-out logistic regression diagnostics
compile("TMB_Files/LRP_Logistic_Only_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/LRP_Logistic_Only_LowAggPrior"))

compile("TMB_Files/LRP_Logistic_Only.cpp")
dyn.load(dynlib("TMB_Files/LRP_Logistic_Only"))
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

# Calculate penalty for sub-population approach
subPop_B_penalty_lwr<-1000 # set at abundance below which no one CU could be above subpop have at least half of subpops above 1000 fish
subPop_B_penalty_upr<-8000 # set at minimum abundance at which all 5 CUs could have half of subpops above 1000 fish 
      # (i.e., MFr = 1 subpop>1000, FCany=1 subpop>1000, LThomp=1 subpop>1000, NThomp=2 subpop>3000, SThomp = 2 subpop> 1000)
B_penalty_mu<-mean(c(subPop_B_penalty_lwr,subPop_B_penalty_upr))
dum<-optim(par=200, fn = getSD, method="Brent",lower=1, upper=5000, low_lim=subPop_B_penalty_lwr, hi_lim=subPop_B_penalty_upr)
B_penalty_sigma<-dum$par

TMB_Inputs_Subpop <- list(Scale = 1000, B_penalty_mu=B_penalty_mu, B_penalty_sigma=B_penalty_sigma)



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
  ps <- c(0.6, 0.7, 0.8, 0.9, 0.99)
  for(pp in 1:length(ps)){
    
    
    # Run with Binomial LRP model with hierarchical Ricker 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_HierRicker_Surv_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bin.HierRickerSurv_",ps[pp]*100, sep=""),
                  bootstrapMode = F, plotLRP=T, runLogisticDiag=T)
    
    # Run with Binomial LRP model with individual model Ricker 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_IndivRicker_Surv_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bin.IndivRickerSurv_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T,runLogisticDiag=T)

    # Run with Binomial LRP model with hierarchical Ricker, with prior on capacity 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_HierRicker_SurvCap_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_HM_priorCap, outDir=cohoDir, RunName = paste("Bin.HierRickerSurvCap_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T,runLogisticDiag=T)
  
    # Run with Binomial LRP model with individual model Ricker, with prior on capacity 
     runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                    BMmodel = "SR_IndivRicker_SurvCap_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
                    useGenMean=F, TMB_Inputs=TMB_Inputs_IM_priorCap, outDir=cohoDir, RunName = paste("Bin.IndivRickerSurvCap_",ps[pp]*100, sep=""),
                    bootstrapMode = F, plotLRP=T,runLogisticDiag=T)
    
     # # Run with Bernoulli LRP model with individual model Ricker
     runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                    BMmodel = "SR_IndivRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
                    useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurv_",ps[pp]*100, sep=""),
                    bootstrapMode = F, plotLRP=T,runLogisticDiag=T)
     
     # Run with Bernoulli LRP model with hierarchical Ricker
     runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                    BMmodel = "SR_HierRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
                    useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bern.HierRickerSurv_",ps[pp]*100, sep=""),
                    bootstrapMode = F, plotLRP=T,runLogisticDiag=T)
     
  }



# Run annual restrospective analyses using subpopulations ===========================

ps <- c(0.6, 0.7, 0.8, 0.9, 0.99)
for(pp in 1:length(ps)){

runAnnualRetro(EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
               BMmodel = "ThreshAbund_Subpop1000_ST", LRPmodel="BernLogistic", LRPfile="LRP_Logistic_Only",integratedModel=F,
               useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bern.SPopAbundThreshST_",ps[pp]*100, sep=""),
               bootstrapMode = F, plotLRP=T)


runAnnualRetro(EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
               BMmodel = "ThreshAbund_Subpop1000_ST", LRPmodel="BinLogistic", LRPfile="LRP_Logistic_Only_LowAggPrior",integratedModel=F,
               useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir, RunName = paste("Bin.SPopAbundThreshST_",ps[pp]*100, sep=""),
               bootstrapMode = F, plotLRP=T,runLogisticDiag=T)

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





# ================================================================================
#Summarize and save logistic model fit diagnostics 
# ===============================================================================

# Create vectors of diagnostic stats to fill for each model
modelName<-NA
PearChiSq<-NA
DevChiSq<-NA
quasiR2<-NA
Wald<-NA
hitRatio<-NA


year <- 2015
p<-0.99
useBern_Logistic <- FALSE

# Model 1 =====================
mm<-1
modelName[mm]<-"IM"
# Specify model fit to evaluate 
BMmodel<-"SR_IndivRicker_Surv_LowAggPrior"
RunName <- paste("Bin.IndivRickerSurv_",p*100, sep="")
TMB_Inputs<-TMB_Inputs_IM
outputDir <- paste(cohoDir,"DataOut/AnnualRetrospective", RunName, sep="/")

# Extract logistic regression data from .rda file
load(paste(outputDir,"/logisticFitDiagStats_",year,".rda",sep=""))
PearChiSq[mm]<-round(logisticDiagStats$p.PearChiSq,2)
DevChiSq[mm]<-round(logisticDiagStats$p.DevChiSq,2)
quasiR2[mm]<-round(logisticDiagStats$quasiR2,2)
Wald[mm]<-round(logisticDiagStats$p.Wald,2)

# Run leave-one-out logistic model diagnostic:
hitRatio[mm]<-LOO_LRdiagnostics_cohoModel(year=year, p = p, useBern_Logistic=useBern_Logistic,
                    RunName=RunName, outputDir = outputDir, TMB_Inputs=TMB_Inputs)



# Model 2 =====================
mm<-2
modelName[mm]<-"HM"
# Specify model fit to evaluate 
BMmodel<-"SR_HierRicker_Surv_LowAggPrior"
RunName <- paste("Bin.HierRickerSurv_",p*100, sep="")
TMB_Inputs<-TMB_Inputs_HM
outputDir <- paste(cohoDir,"DataOut/AnnualRetrospective", RunName, sep="/")

# Extract logistic regression data from .rda file
load(paste(outputDir,"/logisticFitDiagStats_",year,".rda",sep=""))
PearChiSq[mm]<-round(logisticDiagStats$p.PearChiSq,2)
DevChiSq[mm]<-round(logisticDiagStats$p.DevChiSq,2)
quasiR2[mm]<-round(logisticDiagStats$quasiR2,2)
Wald[mm]<-round(logisticDiagStats$p.Wald,2)

# Run leave-one-out logistic model diagnostic:
hitRatio[mm]<-LOO_LRdiagnostics_cohoModel(year=year, p = p, useBern_Logistic=useBern_Logistic,
                                          RunName=RunName, outputDir = outputDir, TMB_Inputs=TMB_Inputs)


# Model 3 =====================
mm<-3
modelName[mm]<-"IM_Cap"
# Specify model fit to evaluate 
BMmodel<-"SR_IndivRicker_SurvCap_LowAggPrior"
RunName <- paste("Bin.IndivRickerSurvCap_",p*100, sep="")
TMB_Inputs<-TMB_Inputs_IM_priorCap
outputDir <- paste(cohoDir,"DataOut/AnnualRetrospective", RunName, sep="/")

# Extract logistic regression data from .rda file
load(paste(outputDir,"/logisticFitDiagStats_",year,".rda",sep=""))
PearChiSq[mm]<-round(logisticDiagStats$p.PearChiSq,2)
DevChiSq[mm]<-round(logisticDiagStats$p.DevChiSq,2)
quasiR2[mm]<-round(logisticDiagStats$quasiR2,2)
Wald[mm]<-round(logisticDiagStats$p.Wald,2)

# Run leave-one-out logistic model diagnostic:
hitRatio[mm]<-LOO_LRdiagnostics_cohoModel(year=year, p = p, useBern_Logistic=useBern_Logistic,
                                          RunName=RunName, outputDir = outputDir, TMB_Inputs=TMB_Inputs)



# Model 4 =====================
mm<-4
modelName[mm]<-"HM_Cap"
# Specify model fit to evaluate 
BMmodel<-"SR_HierRicker_SurvCap_LowAggPrior"
RunName <- paste("Bin.HierRickerSurvCap_",p*100, sep="")
TMB_Inputs<-TMB_Inputs_HM_priorCap
outputDir <- paste(cohoDir,"DataOut/AnnualRetrospective", RunName, sep="/")

# Extract logistic regression data from .rda file
load(paste(outputDir,"/logisticFitDiagStats_",year,".rda",sep=""))
PearChiSq[mm]<-round(logisticDiagStats$p.PearChiSq,2)
DevChiSq[mm]<-round(logisticDiagStats$p.DevChiSq,2)
quasiR2[mm]<-round(logisticDiagStats$quasiR2,2)
Wald[mm]<-round(logisticDiagStats$p.Wald,2)

# Run leave-one-out logistic model diagnostic:
hitRatio[mm]<-LOO_LRdiagnostics_cohoModel(year=year, p = p, useBern_Logistic=useBern_Logistic,
                                          RunName=RunName, outputDir = outputDir, TMB_Inputs=TMB_Inputs)




# Model 5 =====================
mm<-5
modelName[mm]<-"SubPop"
# Specify model fit to evaluate 
#BMmodel<-"SR_HierRicker_SurvCap_LowAggPrior"
RunName <- paste("Bin.SPopAbundThreshST_",p*100, sep="")
TMB_Inputs<-TMB_Inputs_Subpop
outputDir <- paste(cohoDir,"DataOut/AnnualRetrospective", RunName, sep="/")

# Extract logistic regression data from .rda file
load(paste(outputDir,"/logisticFitDiagStats_",year,".rda",sep=""))
PearChiSq[mm]<-round(logisticDiagStats$p.PearChiSq,2)
DevChiSq[mm]<-round(logisticDiagStats$p.DevChiSq,2)
quasiR2[mm]<-round(logisticDiagStats$quasiR2,2)
Wald[mm]<-round(logisticDiagStats$p.Wald,2)

# Run leave-one-out logistic model diagnostic:
hitRatio[mm]<-LOO_LRdiagnostics_cohoModel(year=year, p = p, useBern_Logistic=useBern_Logistic,
                                          RunName=RunName, outputDir = outputDir, TMB_Inputs=TMB_Inputs)


# Summarize model fit diagnostics

diagStats_byModel<-data.frame(model = modelName, PearChiSq=PearChiSq, DevChiSq=DevChiSq,
                              quasiR2=quasiR2, Wald=Wald, hitRatio=hitRatio)