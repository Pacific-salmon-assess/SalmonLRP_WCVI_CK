
library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
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
  source("getMultiDimStatus.r")
}

sourceAll()

# ======================================================================
# (1) Read-in and format IF Coho data:  
# =====================================================================
setwd(cohoDir)

CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")

# LW: Add wild escapement from State of the Salmon Rapid Assessment repo
# wild <- read.csv("https://raw.githubusercontent.com/BronwynMacDonald/Scanner-Data-Processing/main/DATA_OUT/MERGED_FLAT_FILE_BY_CU.csv")
# wild$CU_Name <- sub(" ", "_", wild$CU_Name) # replace space with underscore for merge
# CoEscpDat$Species <- "Coho" # add species column for merge
# # Merge total escapement data with wild data
# CoEscpDat <- merge(CoEscpDat, wild[,names(wild) %in% c("Species", "Year", "SpnForAbd_Wild", "CU_Name")], by.x=c("Species", "CU_Name", "ReturnYear"), by.y=c("Species", "CU_Name", "Year"))
# # Rename to wild escapement
# names(CoEscpDat)[grep("SpnForAbd_Wild", names(CoEscpDat))] <- "Wild_Escapement"
# CoEscpDat <- CoEscpDat[ , -grep("^Escapement$", names(CoEscpDat))] # remove total escapement column
# write.csv(CoEscpDat, "DataIn/IFCoho_escpByCU_SOTSwild.csv")

# Read in data with wild escapement
CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU_SOTSwild.csv")

# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"
colnames(CoEscpDat)[colnames(CoEscpDat)=="Wild_Escapement"] <- "Escp"

# Run annual restrospective analyses over various levels of p, restrict escapement dataset to 1998+ ============================================================
# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)


# ==================================================================================
# (2) Get multi-dimensional status
# ====================================================================================



#BMsource<-"Ricker"
BMsource<-"Ricker_priorCap"
BMsource<-"WSP2015"


statusDat<- CoEscpDat %>% group_by(CU_Name) %>%
  mutate(Gen_Mean = rollapply(Escp, 3, gm_mean, fill = NA, align="right"))

if (BMsource == "Ricker" | BMsource == "WSP2015") {
  SR_pars<-read.csv(paste(cohoDir,"DataOut/ModelFits/AllEsts_Indiv_Ricker_Surv.csv",sep="/"))
}
  
if (BMsource == "Ricker_priorCap") {
  SR_pars<-read.csv(paste(cohoDir,"DataOut/ModelFits/AllEsts_Indiv_Ricker_Surv_priorCap.csv",sep="/"))
}
  
# Get Sgen estimates for IM model
SR_pars<-as_tibble(SR_pars)
Sgen <- SR_pars %>% filter(Param == "Sgen") %>% select(Estimate, CU_Name)
Sgen <- rename(Sgen, Sgen = Estimate)
Smsy <- SR_pars %>% filter(Param == "Smsy") %>% select(Estimate, CU_Name)
Smsy<- rename(Smsy, Smsy = Estimate)
Smsy<- Smsy %>% add_column(Smsy0.8 = Smsy$Smsy*0.8) %>% select(CU_Name, Smsy0.8)


if (BMsource == "WSP2015") {
# # Test: Sub in Sgen estimates from WSP 2015
 Sgen$Sgen<-c(1585, 741, 1405, 2546, 2337)
 Smsy$Smsy0.8<-c(2785, 1562, 3052, 5285, 4462)
}
 
statusDat <- left_join(statusDat, Sgen, by="CU_Name")
statusDat <- left_join(statusDat, Smsy, by="CU_Name")
statusDat<-na.omit(statusDat)


statusEst<-getMultiDimStatus(statusDat = statusDat, outDir=cohoDir, filename=paste("multiDimStatusEsts_",BMsource,".csv", sep=""))


plot_CU_Escp_withMultiStatus(statusEst, outDir=cohoDir, plotName=paste("coho-CU-EscpSeries-wMultiStatus",BMsource, sep="_"))
