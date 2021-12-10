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
library(tmbstan)
library(here)
library(zoo)

setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
cohoDir<-paste(rootDir,"/IFCohoStudy",sep="")

setwd(codeDir)

sourceAll <- function(){
  source("plotFunctions.r")
  source("helperFunctions.r")
}

sourceAll()

# ======================================================================
# Read-in Coho data:  
# =====================================================================
setwd(cohoDir)

# CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
# # Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
# colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
# colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
# colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
# colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"

CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")

CoEscpDat_bySubpop<-read.csv("DataIn/IFCoho_escpBySubpop.csv")

# Change column names to yr, CU_Name, Escp, Subpop_Name, CU_ID
CoEscpDat_bySubpop<-CoEscpDat_bySubpop %>% select(MU_Name=MU_Name, yr=Return.Year, CU_Name=Conservation.Unit, Escp=Natural.Returns, Subpop_Name=Sub.Population)
tmp.df<-data.frame(CU_Name=unique(CoEscpDat_bySubpop$CU_Name), CU_ID=seq(1,length(unique(CoEscpDat_bySubpop$CU_Name)),by=1))
CoEscpDat_bySubpop <- left_join(CoEscpDat_bySubpop,tmp.df)


# Summary of CU-level escapements based on Natural.Spawners  
CoEscpDat <- as_tibble(CoEscpDat_bySubpop) %>% group_by(MU_Name, CU_Name, CU_ID, yr) %>% select(MU_Name, CU_Name, CU_ID, yr, Escp) %>%summarize(Escp=sum(Escp))
CoEscpDat <- as.data.frame(CoEscpDat)
# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"




# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)
CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)
CoEscpDat_bySubpop <- CoEscpDat_bySubpop %>% filter(yr >= 1998)

# Roll up escpaments, and get Gen Mean of htat
AggEscp <- CoEscpDat %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Escp, 3, gm_mean, fill = NA, align="right"))

# ===================================================================================
# Plot annual status with bars to show years in which LRP was breached 
# ===================================================================================

# Make a list of all modeled combinations to be shown in plots 
## (note: the output folders/files for each of options must already exist)

modelFitList<-c("Bern.IndivRickerSurv_50",
                "Bern.IndivRickerSurvCap_50",
                "Bern.SPopAbundThreshST_50")

projLRPList<-c("Ricker_50", "Ricker_priorCap_50")

multiDimList<-c("Ricker_50", "Ricker_priorCap_50")



# Specify which thresholds should be tested in the data-based proportional thresholds
ps_Prop<-c(1.0)

# Provide WSP assessment results
WSP_estYr<-2014
WSP_AboveLRP<-TRUE

LRP_estYr<-2020
retroYears<-2000:LRP_estYr

plotStatusBarsCoho_byYear(LRP_estYr, retroYears, Dir=cohoDir, genYrs=3, AggEscp, EscpDat=CoEscpDat, 
                  modelFitList=modelFitList, projLRPList=projLRPList, multiDimList=multiDimList, ps_Prop=ps_Prop,
                  WSP_estYr=WSP_estYr, WSP_AboveLRP=WSP_AboveLRP,
                  outDir = cohoDir, fName = paste("coho-statusPlot_withBars",LRP_estYr,sep=""))


