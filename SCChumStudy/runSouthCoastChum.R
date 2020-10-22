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
CUdf$CU_ID <- substr(CUdf$CU_raw, 1,1) # LW added: add a column with just CU ID, match style of IFCoho case study 


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

##################################
### Start of L. Warkentin Code ###
##################################
# Read in chum wild escapement data (infilled) by CU
ChumEscDat <- read.csv("DataOut/WildEscape_w.Infill_ByCU.csv")
# Add MU column - FLAG this may not be necessary - can't find MU in any of the function scripts. Check with Kendra
ChumEscDat$MU_Name <- "SC Chum"
head(ChumEscDat)


# merge escapement data with CU name/id lookup table
ChumEscDat <- merge(ChumEscDat, CUdf, by.x=(""))

# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
colnames(ChumEscDat)[] <- c("CU", "yr", "Ecsp")

  colnames(ChumEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
  colnames(ChumEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
  colnames(ChumEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
  colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"
















