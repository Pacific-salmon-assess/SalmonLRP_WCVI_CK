# Code by L. Warkentin & K. Holt
# Started on: October 16, 2020

# This file contains the code required to explore South Coast Chum escapement datasets and run retrospective analyses of 
#   data-based LRPs that use logistic regressions

library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(dplyr)
library(tidyr)
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


# ===========================================================================================
# Read-in Chum Data & Run Data Prep Functions:  
# ===========================================================================================
setwd(chumDir)

# Note: The following infilling & brood table reconstruction procedures are taken from: 
    # Holt, C.A., Davis, B., Dobson, D., Godbout, L., Luedke, W., Tadey, J., and Van Will, P. 2018.
      # Evaluating Benchmarks of Biological Status for Data-limited Conservation Units of Pacific
      # Salmon, Focusing on Chum Salmon in Southern BC. DFO Can. Sci. Advis. Sec. Res. Doc.
      # 2018/011. ix + 77 p. Available at: https://waves-vagues.dfo-mpo.gc.ca/Library/40759386.pdf

# All code in this section (including code in chumDataFunctions.r) was written by B. Davis (DFO) for the  
# above paper, and provided to us by Carrie Holt in October 2020 as part of "Retrospective Analysis BD" folder.

source("chumDataFunctions.r")

# Create look-up table for CU names
RawDat <- read.csv("DataIn/Chum Escapement Data With Areas_2013.csv", check.names=F)
CU_short <- c("SCS", "NEVI", "UK", "LB", "BI", "GS", "HSBI")
CU_raw <- unique(RawDat$CU_Name)
CU_names<-c("Southern Coastal Streams", "North East Vancouver Island", "Upper Knight", 
            "Loughborough", "Bute Inlet", "Georgia Strait", "Howe Sound to Burrard Inlet" )
CUdf <- data.frame(CU_short, "CU_raw"=CU_raw[1:7], CU_names)


## Step 1: Infill Escapement Data ===============================================================

# Infill missing values using only wild fish 
#  --> Code taken from "Infilling.R" file written by B. Davis for Holt et al. 2018 Chum CSAS paper
#  --> Infilling.r file was from "Retrospective Analysis BD" folder)

NoQPDatBySite <- RawDat[which(RawDat$Source=="Wild" & RawDat$NME %notin% c("QUALICUM RIVER", "LITTLE QUALICUM RIVER", "PUNTLEDGE RIVER")),]

#need to change to character to preserve names
NoQPDatBySite$NME <- as.character(NoQPDatBySite$NME)
NoQPDatBySite$CU_Name <- as.character(NoQPDatBySite$CU_Name)
NoQPDatBySite$GroupName <- as.character(NoQPDatBySite$GroupName)
NoQPDatBySite$GU_Name <- as.character(NoQPDatBySite$GU_Name)
NoQPDatBySite$Rabcode <- as.character(NoQPDatBySite$Rabcode)
#SummerRun not registering True false arghhhhh 
NoQPDatBySite$SummerRun[which(is.na(NoQPDatBySite$SummerRun))] <- FALSE
# Remove summer run fish
NoQPDatBySite2 <- NoQPDatBySite[which(!NoQPDatBySite$SummerRun), ]

# Change to long form from wide form
LongDatNoQP <- NoQPDatBySite2 %>% gather( "Year", "Escape", 10:70)

#Now Infill
NoQPDat <- Infill(data = LongDatNoQP, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME")

# Remove Fraser, make data frame to feed into infill again
NoQPDatSumm <- as.data.frame(NoQPDat[[2]][which(NoQPDat[[2]]$CU_Name %in% CUdf$CU_raw),])

# Now infill missing years for entire CU
#Infill missing values
NoQPByCU <- Infill(data=NoQPDatSumm, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc")

# Write by CU output
write.csv(NoQPByCU[[1]], "DataOut/WildEscape_w.Infill_ByCU.csv", row.names = F)


## Step 2: Construct Spawner Recruit Brood Tables =============================================

# Prep infilled escapement data
WildEsc <- NoQPByCU[[1]]
  WildEsc$CUinfill <- ifelse(is.na(WildEsc$Escape), TRUE, FALSE) # Flag 'Escape = NA' sites
  
# Read in Return data from P.V.W received 3/7/2016
WildRetWide <- read.csv("DataIn/WildReturnsPVW_2013.csv", check.names=F)
  # Need to get PVW return data into long form:
    # First need to "collapse" areas to CU level
    #WildByCU <- WildRetWide  %>% group_by(CU_Name)   %>%  summarise_each(funs(sum(., na.rm=T)), 3:63)
    # KH changed above line to:
    WildByCU <- WildRetWide  %>% group_by(CU_Name)   %>%  summarise_at(vars("2013": "1953"),sum,na.rm = T)
    WildByCU$CU_Name <- as.character(WildByCU$CU_Name)
    # Now into long form
    WildRetLong <- WildByCU %>% gather("Year", "Return", 2:62)

# Read in age comp data
ACdat <- read.csv("DataIn/AgeComp_2013.csv")

# -> Code taken from "SCChumWildData.R" file written by B. Davis for Holt et al. 2018 Chum CSAS paper:

# Construct Brood table

#Will need to merge EscDat and ECdat to get catch info, make three Brood tables
# Year needs to be integer to join, stupid work around, couldn't figure out other way
WildRetLong$Year <- as.character(WildRetLong$Year)
WildRetLong$Year <- as.integer(WildRetLong$Year)
# Btable1 <- left_join(data.frame(Year=WildEsc$Year, CU=WildEsc$CU_Name, Escape=WildEsc$SiteEsc, CUinfill=WildEsc$CUinfill), 
#                      data.frame(Year=WildRetLong$Year, CU=WildRetLong$CU, Return=WildRetLong$Return), by=c("Year", "CU"))
# KH - changed above line to:
Btable1 <- left_join(data.frame(Year=as.integer(WildEsc$Year), CU=WildEsc$CU_Name, Escape=WildEsc$SiteEsc, CUinfill=WildEsc$CUinfill), 
                     data.frame(Year=WildRetLong$Year, CU=WildRetLong$CU_Name, Return=WildRetLong$Return), by=c("Year", "CU"))

# KH - added following line to change year to factor:
ACdat$Year<-as.integer(ACdat$Year)

Btable <- left_join(Btable1, ACdat, by=c("Year"))

years <- sort(unique(Btable$Year), decreasing=F)
nyears <- length(years)
sites <- unique(Btable$CU)
nsites <- length(sites)

#need to enter age comp ages
ages<-c(3,4,5,6)
nages<-length(ages)

# go through years and calculate recruits by brood year
# Same age comp used for each CU
# will not get recruit estimate for first two years due to missing age comp data (NAs)
# will not get recruit estiamate for last 
Btable$Recruit <- rep(NA, dim(Btable)[1])
for( i in 1:nsites){
  Sdat<- Btable[which(Btable$CU==sites[i]),]
  #cannot estimate returns for last 5 years
  for( j in 3:nyears ){
    # can only calculate up to nyears-6 -- only have age comps up to 2012
    if(j<=nyears-6){
      Rsum <- 0
      for( k in 1:nages ){
        #add up recruits from brood years
        Rsum <- Rsum + Sdat$Return[which(Sdat$Year==(years[j]+ages[k]))] * Sdat[which(Sdat$Year==years[j]+ages[k]), paste("Age", ages[k], sep="")]
      }
    } else {
      Rsum<-NA
    }
    Btable$Recruit[which(Btable$Year==years[j] & Btable$CU==sites[i])] <- Rsum
  }
}


write.csv(Btable, "DataOut/SRdatWild.csv", row.names = F)

# ====================================================================
# Start of L. Warkentin Code 
# ====================================================================
# Read in data and format for using in retrospective analysis

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
# no sub-population analysis for Chum
#plot_Subpop_Escp_Over_Time(CoEscpDat_bySubpop, cohoDir, plotName="IFC Esc Separate - by Subpop", samePlot = F)
#plot_Subpop_Escp_Over_Time(CoEscpDat_bySubpop, cohoDir, plotName="IFC Esc - by Subpop", samePlot = T)


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


TMB_Inputs_Subpop <- list(Scale = 1000)

# TMB_Inputs_Subpop <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
#                           logMuA_sig = 2, Tau_dist = 0.1, Tau_A_dist = 0.1, 
#                           gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)



# Run annual restrospective analyses using CUs ===========================

# Loop over p values and run annual retrospective analyses for each level of p
ps <- 0.95 # for now just use one p value
# ps <- c(seq(0.6, 0.95,.05), 0.99) 

# LW Notes
# Choosing values for runAnnualRetro function:
# - genYrs: average age of return.  runFraserCoho.r uses 3. I think I should use 4 for chum (but check)
# - BroodYrLag: How many years will it take between a fish returning and when you can construct the full recruitment.
#     The BroodYrLag parameter used when prepping for the LRP calculations represents the number of years 
#     it will take to reconstruct recruitment from a single brood year  (i.e., how many years of recruitment 
#     do you need to observe to estimate recruitment from a brood year).  Another way of defining it is the 
#     number of age classes that recruitment (maturation) occurs over. 
#     For the coho example, fish mostly return at age 3 or 4, so BroodYrLag is 2.
#     For the chum example, it looks like fish return at ages 3,4,5,6, so BroodYrLag is 4.
#     This parameter gets used in the runAnnualRetro() function to remove the first few years of data for 
#     which recruitment is NA because the BY has not yet been fully observed.  So, basically it's automating
#     the step that you just mentioned you were doing when prepping the data for model fitting.
#
# The Chum data doesn't have smolt to adult survival (STAS) data, like the coho data does.
#     This means that a new TMB model file needs to be written (e.g., SR_IndivRicker.cpp) that doesn't have the 
#     STAS_Age_<return age> input from <species>SRDat, muLSurv variable (in LRPFunctions.r), and muLSurv parameter 
#     in SR_IndivRicker_Surv.cpp and other TMB models.


for(pp in 1:length(ps)){
  # # Run with Binomial LRP model with hierarchical Ricker 
  # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
  #                BMmodel = "SR_HierRicker_Surv", LRPmodel="BinLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bin.HierRickerSurv_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
  # 
  # # Run with Bernoulli LRP model with hierarchical Ricker 
  # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
  #                BMmodel = "SR_HierRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_HM, outDir=cohoDir, RunName = paste("Bern.HierRickerSurv_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
  # 
  # Run with Binomial LRP model with individual model Ricker
  runAnnualRetro(EscpDat=ChumEscpDat, SRDat=ChumSRDat, startYr=1960, endYr=2017, BroodYrLag=4, genYrs=4, p = ps[pp],
                 BMmodel = "SR_IndivRicker_Surv", LRPmodel="BinLogistic", integratedModel=T,
                 useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRickerSurv_",ps[pp]*100, sep=""),
                 bootstrapMode = F, plotLRP=T)
  # 
  # # Run with Bernoulli LRP model with individual model Ricker 
  # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurv_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
  # 
  # # Run with Binomial LRP model with hierarchical Ricker, with prior on capacity 
  # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
  #                BMmodel = "SR_HierRicker_SurvCap", LRPmodel="BinLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_HM_priorCap, outDir=cohoDir, RunName = paste("Bin.HierRickerSurvCap_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
  # 
  # 
  # # Run with Bernoulli LRP model with hierarchical Ricker, with prior on capacity 
  # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
  #                BMmodel = "SR_HierRicker_SurvCap", LRPmodel="BernLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_HM_priorCap, outDir=cohoDir, RunName = paste("Bern.HierRickerSurvCap_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
  # 
  # # Run with Binomial LRP model with individual model Ricker, with prior on capacity 
  # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_SurvCap", LRPmodel="BinLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM_priorCap, outDir=cohoDir, RunName = paste("Bin.IndivRickerSurvCap_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
  # 
  # 
  # # Run with Bernoulli LRP model with individual model Ricker, with prior on capacity 
  # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_SurvCap", LRPmodel="BernLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM_priorCap, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurvCap_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
  # 
  # 
}
