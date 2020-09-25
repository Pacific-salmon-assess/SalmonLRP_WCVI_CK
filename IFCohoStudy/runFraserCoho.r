# Code by K.Holt & B.Davis
# Started on: February 7, 2020

# This file contains the code required to explore Interior Fraser Coho escapement datasets and run retrospective analyses of 
#   data-based LRPs that use logistic regressions
library(MASS) # dose.p function to get SE around P95
library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)

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
}
sourceAll()

# Load TMB models

compile("TMB_Files/ThreshAbund_Subpop1000.cpp")
dyn.load(dynlib("TMB_Files/ThreshAbund_Subpop1000"))

compile("TMB_Files/SR_HierRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv"))

compile("TMB_Files/SR_IndivRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv"))

compile("TMB_Files/SR_HierRicker_SurvCap.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_SurvCap"))

compile("TMB_Files/SR_IndivRicker_SurvCap.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_SurvCap"))


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
plot_CU_DataObs_Over_Time(CoEscpDat, cohoDir, plotName="Fr_Co_DataByCU")
plot_Num_CUs_Over_Time(CoEscpDat, cohoDir, plotName="Fr_Co_N_CUs")

# Note: these next 2 two escpt plots need to have formatting fixed
plot_CU_Escp_Over_Time(CoEscpDat, cohoDir, plotName="IFC Esc", samePlot = T)
plot_CU_Escp_Over_Time(CoEscpDat, cohoDir, plotName="IFC Esc Separate", samePlot = F)
plot_Subpop_Escp_Over_Time(CoEscpDat_bySubpop, cohoDir, plotName="IFC Esc Separate - by Subpop", samePlot = F)
plot_Subpop_Escp_Over_Time(CoEscpDat_bySubpop, cohoDir, plotName="IFC Esc - by Subpop", samePlot = T)


# ==================================================================================================================
# Run retrospective analyses:
# =====================================================================================================================

# TMB input parameters:
TMB_Inputs_HM <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)


TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)


# Prior means come from running "compareRickerModelTypes.r"
cap_priorMean<-c(10.957092, 5.565526, 11.467815, 21.104274, 14.803877)

TMB_Inputs_HM_priorCap <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
                   logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
                   gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                   cap_mean=cap_priorMean, cap_sig=sqrt(2))

# Prior means come from running "compareRickerModelTypes.r"
cap_priorMean<-c(11.153583,  5.714955, 11.535779, 21.379558, 14.889006)

TMB_Inputs_IM_priorCap <- list(Scale = 1000, logA_Start = 1, Tau_dist = 0.1, 
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                               cap_mean=cap_priorMean, cap_sig=sqrt(2))


TMB_Inputs_Subpop <- list(Scale = 1000)

# TMB_Inputs_Subpop <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
#                           logMuA_sig = 2, Tau_dist = 0.1, Tau_A_dist = 0.1, 
#                           gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)



# Run annual restrospective analyses using CUs ===========================

# Loop over p values and run annual retrospective analyses for each level of p
  ps <- c(seq(0.6, 0.95,.05), 0.99)
  for(pp in 1:length(ps)){
    # Run with Binomial LRP model with hierarchical Ricker 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_HierRicker_Surv", LRPmodel="BinLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bin.HierRickerSurv_",ps[pp]*100, sep=""),
                  bootstrapMode = F, plotLRP=T)
  
    # Run with Bernoulli LRP model with hierarchical Ricker 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_HierRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bern.HierRickerSurv_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T)
    
    # Run with Binomial LRP model with individual model Ricker 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_IndivRicker_Surv", LRPmodel="BinLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bin.IndivRickerSurv_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T)
    
    # Run with Bernoulli LRP model with individual model Ricker 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_IndivRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurv_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T)
    print(ps[pp])
    print("Binomial Hier")
    
    # Run with Binomial LRP model with hierarchical Ricker, with prior on capacity 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_HierRicker_SurvCap", LRPmodel="BinLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_HM_priorCap, outDir=cohoDir, RunName = paste("Bin.HierRickerSurvCap_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T)
    
    print("Bernoulli Hier")
  
    
    # Run with Bernoulli LRP model with hierarchical Ricker, with prior on capacity 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_HierRicker_SurvCap", LRPmodel="BernLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_HM_priorCap, outDir=cohoDir, RunName = paste("Bern.HierRickerSurvCap_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T)
    
    print("Binomial Indiv")
    
    # Run with Binomial LRP model with individual model Ricker, with prior on capacity 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_IndivRicker_SurvCap", LRPmodel="BinLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_IM_priorCap, outDir=cohoDir, RunName = paste("Bin.IndivRickerSurvCap_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T)
    
    print("Bernoulli Indiv")
    
    # Run with Bernoulli LRP model with individual model Ricker, with prior on capacity 
    runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
                   BMmodel = "SR_IndivRicker_SurvCap", LRPmodel="BernLogistic", integratedModel=T,
                   useGenMean=F, TMB_Inputs=TMB_Inputs_IM_priorCap, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurvCap_",ps[pp]*100, sep=""),
                   bootstrapMode = F, plotLRP=T)
    
    
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
BMmodelList<-c("SR_HierRicker_Surv","SR_HierRicker_Surv")
LRPmodelList<-c("BinLogistic", "BernLogistic")

nCUs<-c(5,3,4)

runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, 
             genYrs=3, p = 0.9,
             BMmodelList = BMmodelList, LRPmodelList=LRPmodelList, integratedModel=T,
             useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir)

runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, 
             genYrs=3, p = 0.95,
             BMmodelList = BMmodelList, LRPmodelList=LRPmodelList, integratedModel=T,
             useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir)

runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, 
             genYrs=3, p = 0.8,
             BMmodelList = BMmodelList, LRPmodelList=LRPmodelList, integratedModel=T,
             useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir)

# For sub-population approach ============================================
BMmodelList<-c("ThreshAbund_Subpop1000_ST")
LRPmodelList<-c("BernLogistic", "BinLogistic")

runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2015, endYr=2018, BroodYrLag=2, 
             genYrs=3, p = 0.95,
             BMmodelList = BMmodelList, LRPmodelList=LRPmodelList, integratedModel=F,
             useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir)

runNCUsRetro(nCUList=nCUs, EscpDat=CoEscpDat_bySubpop, SRDat=NULL, startYr=2015, endYr=2018, BroodYrLag=2, 
             genYrs=3, p = 0.80,
             BMmodelList = BMmodelList, LRPmodelList=LRPmodelList, integratedModel=F,
             useGenMean=F, TMB_Inputs=TMB_Inputs_Subpop, outDir=cohoDir)



# =======================================================================
# Make Status Plots for Retrospective Analysis: 
# =======================================================================


# Plot aggregate status by method and number of CUs

yearList<-2015:2018

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_HierRicker_Surv",p=0.90, Dir=cohoDir,
                     inputPrefix = "retroStatus_CU.obj",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_HierRicker_Surv",p=0.95, Dir=cohoDir,
                     inputPrefix = "retroStatus_CU.obj",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "SR_HierRicker_Surv",p=0.90, Dir=cohoDir,
                     inputPrefix = "retroStatus_CU.obj",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_HierRicker_Surv",p=0.80, Dir=cohoDir,
                     inputPrefix = "retroStatus_CU.obj",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "SR_HierRicker_Surv",p=0.80, Dir=cohoDir,
                     inputPrefix = "retroStatus_CU.obj",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "ThreshAbund_Subpop1000_ST", p=0.95, Dir=cohoDir,
                     inputPrefix = "retroStatus_SPop.obj",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "ThreshAbund_Subpop1000_ST", p=0.80, Dir=cohoDir,
                     inputPrefix = "retroStatus_SPop.obj",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "ThreshAbund_Subpop1000_ST", p=0.95, Dir=cohoDir,
                     inputPrefix = "retroStatus_SPop.obj",plotAveLine=TRUE)

plotAggStatus_byNCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "ThreshAbund_Subpop1000_ST", p=0.80, Dir=cohoDir,
                     inputPrefix = "retroStatus_SPop.obj",plotAveLine=TRUE)

# Plot CVs on LRP estimates by method and number of CUs
yearList<-2015:2018

plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_HierRicker_Surv", p=0.80, Dir=cohoDir,
                   inputPrefix = "retroStatus_CU.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))

plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_HierRicker_Surv", p=0.90, Dir=cohoDir,
                   inputPrefix = "retroStatus_CU.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))

plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "SR_HierRicker_Surv", p=0.95, Dir=cohoDir,
                   inputPrefix = "retroStatus_CU.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))

plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "SR_HierRicker_Surv", p=0.80, Dir=cohoDir,
                   inputPrefix = "retroStatus_CU.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))

plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BinLogistic", BMmodel = "SR_HierRicker_Surv", p=0.90, Dir=cohoDir,
                   inputPrefix = "retroStatus_CU.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))

plotLRP.CV_by_nCUs(year=yearList, nCUList=c(5,4,3), LRPmodel="BernLogistic", BMmodel = "ThreshAbund_Subpop1000_ST", p=0.95, Dir=cohoDir,
                   inputPrefix = "retroStatus_SPop.obj", fName = paste("statusPlot_byCUs_",plotYear, sep=""))


# Plot annual retrospective status with bars for various options for LRP metric 

# Make a list of all modelled combinations to be shown in plots
# modelFitList<-c("Bern.HierRickerSurv_75",
#                 "Bern.HierRickerSurv_85",
#                 "Bern.HierRickerSurv_95",
#                 "Bin.HierRickerSurv_60",
#                 "Bin.HierRickerSurv_80",
#                 "Bern.SPopAbundThreshST_75",
#                 "Bern.SPopAbundThreshST_85",
#                 "Bern.SPopAbundThreshST_95",
#                 "Bin.SPopAbundThreshST_60",
#                 "Bin.SPopAbundThreshST_80")

modelFitList<-c("Bin.IndivRickerSurv_60",
                "Bin.IndivRickerSurv_80",
                "Bin.IndivRickerSurv_99",
                "Bin.HierRickerSurv_60",
                "Bin.HierRickerSurv_80",
                "Bin.HierRickerSurv_99",
                "Bin.IndivRickerSurvCap_60",
                "Bin.IndivRickerSurvCap_80",
                "Bin.IndivRickerSurvCap_99",
                "Bin.HierRickerSurvCap_60",
                "Bin.HierRickerSurvCap_80",
                "Bin.HierRickerSurvCap_99",
                "Bern.SPopAbundThreshST_75",
                "Bern.SPopAbundThreshST_85",
                "Bern.SPopAbundThreshST_95")

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



