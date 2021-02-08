############################################
# Inner South Coast Chum benchmark analysis
# March 7, 2016
# Estimating recruitment of wild fish only
# Using new estimated returns from P.V.W 
#############################################

# Set working directory
setwd("C:/Users/DavisBr/Documents/Salmon Benchmarks/SCChum")

# load packages
library(dplyr)
library(tidyr)

# Read in wild only escapement
WildEsc <- read.csv("DataOut/Wild_ByCU.csv")

# Only need "SiteEsc" column, but will flag Escape=NA sites
WildEsc$CUinfill <- ifelse(is.na(WildEsc$Escape), TRUE, FALSE)

# Read in Return data from P.V.W received 3/7/2016
WildRetWide <- read.csv("DataIn/WildReturnsPVW.csv", check.names=F)

# Read in age comp data
ACdat <- read.csv("DataIn/AgeComp.csv")

##############################################################

# Need to get PVW return data into long form
# First need to "collapse" areas to CU level
WildByCU <- WildRetWide  %>% group_by(CU_Name)   %>%  summarise_each(funs(sum(., na.rm=T)), 3:63)
WildByCU$CU_Name <- as.character(WildByCU$CU_Name)

# Now into long form
WildRetLong <- WildByCU %>% gather("Year", "Return", 2:62)

# Now construct Brood table

#Will need to merge EscDat and ECdat to get catch info, make three Brood tables
# Year needs to be integer to join, stupid work around, couldn't figure out other way
WildRetLong$Year <- as.character(WildRetLong$Year)
WildRetLong$Year <- as.integer(WildRetLong$Year)
Btable1 <- left_join(data.frame(Year=WildEsc$Year, CU=WildEsc$CU_Name, Escape=WildEsc$SiteEsc, CUinfill=WildEsc$CUinfill), 
                    data.frame(Year=WildRetLong$Year, CU=WildRetLong$CU, Return=WildRetLong$Return), by=c("Year", "CU"))
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


write.csv(Btable, "DataOut/SRdatWild.csv")

