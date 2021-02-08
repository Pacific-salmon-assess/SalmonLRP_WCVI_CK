#################################################
# South Coast Chum infilling 
# Brooke Davis 
# started 9/9/2015
#################################################

#Run this first
setwd("C:/github/RetrospectiveAnalysis_BD/SCChum")

RawDat <- read.csv("DataIn/Chum Escapement Data With Areas.csv", check.names=F)

library(dplyr)
library(tidyr)
FuncPath <- "C:/github/RetrospectiveAnalysis_BD/Functions/Functions.R"
source(FuncPath)

SRdat <- read.csv("DataOut/SRdata.csv")
#need to match CU names
CU_short <- unique(SRdat$CU)
CU_raw <- unique(RawDat$CU_Name)
CU_names<-c("Southern Coastal Streams", "North East Vancouver Island", "Upper Knight", 
            "Loughborough", "Bute Inlet", "Georgia Strait", "Howe Sound to Burrard Inlet" )
CUdf <- data.frame(CU_short, "CU_raw"=CU_raw[1:7], CU_names)

##########################################################################################

## All Data -- wild/enhanced/brood/rack
 
#need to group by site name to collapse any brood/rack/enhanced/rack etc categories
#only want to collapse year columns
AllDatBySite <- RawDat  %>% group_by(CU_Name, GroupName, GU_Name, NME, SummerRun, Rabcode, Area)   %>%  
      summarise_each(funs(sum(., na.rm=T)), 3:63)

#need to turn all 0 Back to NA
AllDatBySite[AllDatBySite==0] <- NA
#need to change to character to preserve names
AllDatBySite$NME <- as.character(AllDatBySite$NME)
AllDatBySite$CU_Name <- as.character(AllDatBySite$CU_Name)
AllDatBySite$GroupName <- as.character(AllDatBySite$GroupName)
AllDatBySite$GU_Name <- as.character(AllDatBySite$GU_Name)
AllDatBySite$Rabcode <- as.character(AllDatBySite$Rabcode)
#SummerRun not registering True false arghhhhh
AllDatBySite$SummerRun[which(is.na(AllDatBySite$SummerRun))] <- FALSE
# Remove summer run fish
AllDatBySite2 <- AllDatBySite[which(!AllDatBySite$SummerRun), ]

# Change to long form from wide form
LongDat <- AllDatBySite2 %>% gather( "Year", "Escape", 8:68)

# 3/29/2016 figure out what proportion of sites are actually observed for each CU/year
ObsDF <- data.frame(CU=character(), Year=numeric(), Nobs=numeric(), Ninfilled=numeric())
CUs <- unique(LongDat$CU_Name)
for(i in 1:length(CUs)){
  CUdat <- LongDat[which(LongDat$CU_Name==CUs[i]),]
  years <- unique(CUdat$Year)
  for(j in 1:length(years)){
    Ydat <- CUdat[which(CUdat$Year==years[j]),]
    Nobs <- sum(is.na(Ydat$Escape)==F)
    Ninfilled <- sum(is.na(Ydat$Escape))
    ObsDF <- rbind(ObsDF, data.frame(CU=CUs[i], Year=years[j], Nobs, Ninfilled))
  }
}
write.csv(ObsDF, "DataOut/InfillObs.csv")

#Now Infill
AllDat <- Infill(data = LongDat, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME")

# Remove Fraser, make data frame to feed into infill again
AllDatSumm <- as.data.frame(AllDat[[2]][which(AllDat[[2]]$CU_Name %in% CUdf$CU_raw),])
AllDatAll <- AllDat[[1]][which(AllDat[[1]]$CU_Name %in% CUdf$CU_raw),c("CU_Name", "NME", "Year", "Props", "SiteEsc", "Area", "Rabcode")]

# Now infill missing years for entire CU

#Infill missing values for all wild/hatchery by CU

AllByCU <- Infill(data=AllDatSumm, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc")

## 9/25/15 Also want all infilled by site
xx <- left_join(AllDatAll, data.frame(CU_Name=AllByCU[[1]]$CU_Name, Year=AllByCU[[1]]$Year, CUEsc = AllByCU[[1]]$SiteEsc))
xx$Escape <- ifelse(is.nan(xx$SiteEsc), xx$CUEsc*xx$Props, xx$SiteEsc)
xx$CU <- CUdf$CU_short[match(xx$CU_Name, CUdf$CU_raw)]
write.csv(xx, "DataOut/AllInfilledBySite.csv", row.names=F)

#Also want sumamrized by CU
write.csv(AllByCU[[1]], "DataOut/AllInfilledByCU.csv")

# And sumamrized by CU x area
AllBySite <- xx
AllByArea <- AllBySite %>% group_by(CU_Name, Area, Year) %>% summarize( AreaTot = sum(Escape))
write.csv(AllByArea, "DataOut/AllInfilledByArea.csv")

#Join with report data
AllByCUSumm <- AllByCU[[1]]
AllByCUSumm$CU <- CUdf$CU_short[match(AllByCUSumm$CU_Name, CUdf$CU_raw)]
#Now export to use in benchmark work
DFout <- data.frame(Year=AllByCUSumm$Year, CU=AllByCUSumm$CU, Escape=AllByCUSumm$SiteEsc)
write.csv(DFout, "DataOut/InfilledAllDat2_15.csv")

# dumb work around, but need years to not be factors to link data sets
AllByCUSumm$Year <- as.character(AllByCUSumm$Year)
AllByCUSumm$Year <- as.numeric(AllByCUSumm$Year)

#now merge two together, need to remove "Escape" Column from output
Infilled <- left_join(AllByCUSumm[,-which(names(AllByCUSumm)=="Escape")],SRdat[,1:4], by=c("CU", "Year"))

pdf("Figures/InFilled_2_26_All.pdf")

par(mfrow=c(2,1), mar=c(3,3,2,2), oma=c(2,2,2,2), mgp=c(2,0.8,0))
loc=c("topright", "topleft")
lind <- c(1,2,1,1,1,2,2)
for(i in 1:length(unique(CU_short))){
  Sdat <- Infilled[which(Infilled$CU==CU_short[i]),]
  ylim=c(min(Sdat$Escape, Sdat$CuEsc, na.rm=T),max(Sdat$Escape, Sdat$SiteEsc, na.rm=T))
  
  plot(x=Sdat$Year, y=Sdat$Escape, col="grey", type="o", pch=19, ylim=ylim, xlab="Year", ylab="Escapement", lwd=1.4, cex=0.7)
  lines(x=Sdat$Year, y=Sdat$SiteEsc, type="o", pch=19, lwd=1.4, cex=0.7)
  legend( loc[lind[i]], col=c("grey", "black"), lty=rep(1,2), pch=rep(19,2), legend=c("Report Data", "Infilled Data"), 
          cex=0.9, bty="n")
  mtext(side=3, text=as.character(CUdf$CU_names[match(CU_short[i], CUdf$CU_short)]))
}

dev.off()

#now my data mostly matches PVW's report, except NEVI 1981, and missing values

########################################################################
# What if use arithmatic mean instead of geometric??

AmeanAll <- Infill(data = LongDat, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME", avg="Amean")

# Remove Fraser, make data frame to feed into infill again
AmeanBySiteAll <- as.data.frame(AmeanAll[[2]][which(AmeanAll[[2]]$CU_Name %in% CUdf$CU_raw),])

# Now infill missing years for entire CU
AmeanAllByCU <- Infill(data=AmeanBySiteAll, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc", avg="Amean")

#Join with report data
AmeanSumm <- AmeanAllByCU[[1]]
AmeanSumm$CU <- CUdf$CU_short[match(AmeanSumm$CU_Name, CUdf$CU_raw)]

# dumb work around, but need years to not be factors to link data sets
AmeanSumm$Year <- as.character(AmeanSumm$Year)
AmeanSumm$Year <- as.numeric(AmeanSumm$Year)

#now merge two together, need to remove "Escape" Column from output
AmeanInfilled <- left_join(AmeanSumm[,-which(names(AmeanSumm)=="Escape")],SRdat[,1:4], by=c("CU", "Year"))


#Now plot report, GeoMean, Amean

pdf("Figures/InFilled9_17_CompMethod.pdf")

par(mfrow=c(2,1), mar=c(3,3,2,2), oma=c(2,2,2,2), mgp=c(2,0.8,0))
loc=c("topright", "topleft")
lind <- c(1,2,1,1,1,2,2)
for(i in 1:length(unique(CU_short))){
  Sdat <- Infilled[which(Infilled$CU==CU_short[i]),]
  AmeanSdat <- AmeanInfilled[which(AmeanInfilled$CU==CU_short[i]),]
  ylim=c( min(Sdat$Escape, Sdat$CuEsc, na.rm=T), max(Sdat$Escape, Sdat$SiteEsc, AmeanSdat$SiteEsc, na.rm=T))
  
  plot(x=Sdat$Year, y=Sdat$Escape, col="grey", type="o", pch=19, ylim=ylim, xlab="Year", ylab="Escapement", lwd=1.4, cex=0.7)
  lines(x=Sdat$Year, y=Sdat$SiteEsc, type="o", pch=19, lwd=1.4, cex=0.7, col="blue")
  lines(x=AmeanSdat$Year, y=AmeanSdat$SiteEsc, type="o", pch=19, lwd=1.4, cex=0.7, col="red")
  legend( loc[lind[i]], col=c("grey", "blue", "red"), lty=rep(1,3), pch=rep(19,3), 
          legend=c("Report Data", "Infilled Geo. Mean", "Infilled Arith. Mean"), 
          cex=0.9, bty="n")
  mtext(side=3, text=as.character(CUdf$CU_names[match(CU_short[i], CUdf$CU_short)]))
}

dev.off()


################################################################################################################

# Now start from beginning and take out all non-wild, re-do infilling

# Start from RawDat, remove non-wild

WildDatBySite <- RawDat[which(RawDat$Source=="Wild"),] # 2/22/16 what about brood/rack? should also be considered enhanced?-- YES

#need to change to character to preserve names
WildDatBySite$NME <- as.character(WildDatBySite$NME)
WildDatBySite$CU_Name <- as.character(WildDatBySite$CU_Name)
WildDatBySite$GroupName <- as.character(WildDatBySite$GroupName)
WildDatBySite$GU_Name <- as.character(WildDatBySite$GU_Name)
WildDatBySite$Rabcode <- as.character(WildDatBySite$Rabcode)
#SummerRun not registering True false arghhhhh 
WildDatBySite$SummerRun[which(is.na(WildDatBySite$SummerRun))] <- FALSE
# Remove summer run fish
WildDatBySite2 <- WildDatBySite[which(!WildDatBySite$SummerRun), ]

# Change to long form from wide form
LongDatWild <- WildDatBySite2 %>% gather( "Year", "Escape", 10:70)

#Now Infill
WildDat <- Infill(data = LongDatWild, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME")
WildDatAmean <- Infill(data = LongDatWild, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME", avg="Amean")
#AllData
#write.csv(WildDat[[1]], "DataOut/InfilledWild9_17_all.csv")
#Sumamarised data
#write.csv(WildDat[[2]], "DataOut/InfilledWild9_17.csv")

# Remove Fraser, make data frame to feed into infill again
WildDatSumm <- as.data.frame(WildDat[[2]][which(WildDat[[2]]$CU_Name %in% CUdf$CU_raw),])
WildDatSummAmean <- as.data.frame(WildDatAmean[[2]][which(WildDatAmean[[2]]$CU_Name %in% CUdf$CU_raw),])
WildDatAll <- WildDat[[1]][which(WildDat[[1]]$CU_Name %in% CUdf$CU_raw),c("CU_Name", "NME", "Year", "Props", "SiteEsc", "Area", "Rabcode")]


# Now infill missing years for entire CU

#Infill missing values
WildByCU <- Infill(data=WildDatSumm, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc")
# also compare arithmatic mean
WildByCUAmean <-  Infill(data=WildDatSummAmean, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc", avg="Amean")

## 9/25/15 Also want all infilled by site
yy <- left_join(WildDatAll, data.frame(CU_Name=WildByCU[[1]]$CU_Name, Year=WildByCU[[1]]$Year, CUEsc = WildByCU[[1]]$SiteEsc))
yy$Escape <- ifelse(is.nan(yy$SiteEsc), yy$CUEsc*yy$Props, yy$SiteEsc)
yy$CU <- CUdf$CU_short[match(yy$CU_Name, CUdf$CU_raw)]
write.csv(yy, "DataOut/WildInfilledBySite.csv")


#Add CU column
WildSumm <- WildByCU[[1]]
WildSumm$CU <- CUdf$CU_short[match(WildSumm$CU_Name, CUdf$CU_raw)]
WildSumm$Year <- as.character(WildSumm$Year)
WildSumm$Year <- as.numeric(WildSumm$Year)
WildSummAmean <- WildByCUAmean[[1]]
WildSummAmean$CU <- CUdf$CU_short[match(WildSummAmean$CU_Name, CUdf$CU_raw)]
WildSummAmean$Year <- as.character(WildSummAmean$Year)
WildSummAmean$Year <- as.numeric(WildSummAmean$Year)

#make simple data frame to export with geo. mean infilling
WildDFout <- data.frame(Year=WildSumm$Year, CU=WildSumm$CU, Escape=WildSumm$SiteEsc)
write.csv(WildDFout, "DataOut/InfilledWild2_15.csv")


#Join All Infilled
# Start with "Infilled"
InfillComp <- data.frame("Year"=Infilled$Year, "CU"=Infilled$CU, "ReportEsc"=Infilled$Escape, "AllGmean" =Infilled$SiteEsc)
# Now add Amean data
InfillComp <- left_join(InfillComp, data.frame("Year"=AmeanSumm$Year, "CU"=AmeanSumm$CU, "AllAmean"=AmeanSumm$SiteEsc), by=c("Year", "CU"))
#Now add Wild Gmean
InfillComp <- left_join(InfillComp, data.frame("Year"=WildSumm$Year, "CU"=WildSumm$CU, "WildGmean"=WildSumm$SiteEsc), by=c("Year", "CU"))
# Wild Amean
InfillComp <- left_join(InfillComp, data.frame("Year"=WildSummAmean$Year, "CU"=WildSummAmean$CU, "WildAmean"=WildSummAmean$SiteEsc), by=c("Year", "CU"))

Cols <- c("#17AB69", "#2120B7",  "#FF314B", "#E8B338")
Trans <- 80
ColsTrans <- paste(Cols, Trans, sep="")

pdf("Figures/InFilled_Comp_All.pdf")

par(mfrow=c(2,1), mar=c(3,3,2,2), oma=c(2,2,2,2), mgp=c(2,0.8,0))
loc=c("topright", "topleft")
lind <- c(1,2,1,1,1,2,2)
for(i in 1:length(unique(CU_short))){
  Sdat <- InfillComp[which(InfillComp$CU==CU_short[i]),]
  ylim=c(min(Sdat$WildGmean, Sdat$WildAmean, na.rm=T), max(Sdat$ReportEsc, Sdat$AllGmean, Sdat$AllAmean,  na.rm=T))
  plot(x=Sdat$Year, y=Sdat$ReportEsc, col="grey", type="o", pch=19, ylim=ylim, xlab="Year", ylab="Escapement", lwd=1.4, cex=0.7)
    lines(x=Sdat$Year, y=Sdat$AllGmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[1])
    lines(x=Sdat$Year, y=Sdat$AllAmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[2])
    lines(x=Sdat$Year, y=Sdat$WildGmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[3])
    lines(x=Sdat$Year, y=Sdat$WildAmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[4])
  #legend is big, only put on one plot per page
  if(i %in% c(1,4,5,7)){
   legend( loc[lind[i]], col=c("grey", ColsTrans), lty=rep(1,5), pch=rep(19,5), legend=c("Report Data", "All Geo. Mean", 
          "All Arith. Mean", "Wild Geo. Mean", "Wild Arith. Mean"), 
          cex=0.8, bty="n")
  }
  mtext(side=3, text=as.character(CUdf$CU_names[match(CU_short[i], CUdf$CU_short)]))
}

dev.off()


################################################################################
#Same as above but also remove Little Qualicum, Qualicum, and Puntledge Rivers entirely

# Start from RawDat, remove non-wild

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
#AllData
#write.csv(NoQPDat[[1]], "DataOut/InfilledNoQP9_17_all.csv")
#Sumamarised data
#write.csv(NoQPDat[[2]], "DataOut/InfilledNoQP9_17.csv")

# Remove Fraser, make data frame to feed into infill again
NoQPDatSumm <- as.data.frame(NoQPDat[[2]][which(NoQPDat[[2]]$CU_Name %in% CUdf$CU_raw),])
NoQPDatAll <- NoQPDat[[1]][which(NoQPDat[[1]]$CU_Name %in% CUdf$CU_raw),c("CU_Name", "NME", "Year", "Props", "SiteEsc", "Area", "Rabcode")]

# Now infill missing years for entire CU

#Infill missing values
NoQPByCU <- Infill(data=NoQPDatSumm, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc")
NoQPByCUAmean <-  Infill(data=NoQPDatSummAmean, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc", avg="Amean")

## 9/25/15 Also want all infilled by site
zz <- left_join(NoQPDatAll, data.frame(CU_Name=NoQPByCU[[1]]$CU_Name, Year=NoQPByCU[[1]]$Year, CUEsc = NoQPByCU[[1]]$SiteEsc))
zz$Escape <- ifelse(is.nan(zz$SiteEsc), zz$CUEsc*zz$Props, zz$SiteEsc)
zz$CU <- CUdf$CU_short[match(zz$CU_Name, CUdf$CU_raw)]
write.csv(zz, "DataOut/NoQPInfilledBySite.csv")

# Also by CU
write.csv(NoQPByCU[[1]], "DataOut/Wild_ByCU.csv")

#Now by CU, Area
AllBySiteWild <- zz
AllByAreaWild <- AllBySiteWild %>% group_by(CU_Name, Area, Year) %>% summarize( AreaTot = sum(Escape))
write.csv(AllByAreaWild, "DataOut/WildBYarea.csv")


#Add CU column
NoQPSumm <- NoQPByCU[[1]]
NoQPSumm$CU <- CUdf$CU_short[match(NoQPSumm$CU_Name, CUdf$CU_raw)]
NoQPSumm$Year <- as.character(NoQPSumm$Year)
NoQPSumm$Year <- as.numeric(NoQPSumm$Year)
NoQPSummAmean <- NoQPByCUAmean[[1]]
NoQPSummAmean$CU <- CUdf$CU_short[match(NoQPSummAmean$CU_Name, CUdf$CU_raw)]
NoQPSummAmean$Year <- as.character(NoQPSummAmean$Year)
NoQPSummAmean$Year <- as.numeric(NoQPSummAmean$Year)

#simple data frame for export
NoQPDFout <- data.frame(Year=NoQPSumm$Year, CU=NoQPSumm$CU, Escape=NoQPSumm$SiteEsc)
write.csv(NoQPDFout, "DataOut/InfilledWildNoQP.csv")

#make plot showing all, wild, this
#Join data together
# Start with "Infilled"
InfillWild <- data.frame("Year"=Infilled$Year, "CU"=Infilled$CU, "ReportEsc"=Infilled$Escape, "AllGmean" =Infilled$SiteEsc)
#Now add Wild Gmean
InfillWild <- left_join(InfillWild, data.frame("Year"=WildSumm$Year, "CU"=WildSumm$CU, "WildGmean"=WildSumm$SiteEsc), by=c("Year", "CU"))
# Wild No QP
InfillWild <- left_join(InfillWild, data.frame("Year"=NoQPSumm$Year, "CU"=NoQPSumm$CU, "NoQPGmean"=NoQPSumm$SiteEsc), by=c("Year", "CU"))


Cols <- c("#17AB69", "#2120B7",  "#FF314B", "#E8B338")
Trans <- 80
ColsTrans <- paste(Cols, Trans, sep="")

pdf("Figures/InFilled_CompWild.pdf")

par(mfrow=c(2,1), mar=c(3,3,2,2), oma=c(2,2,2,2), mgp=c(2,0.8,0))
loc=c("topright", "topleft")
lind <- c(1,2,1,1,1,2,2)
for(i in 1:length(unique(CU_short))){
  Sdat <- InfillWild[which(InfillWild$CU==CU_short[i]),]
  ylim=c(min(Sdat$NoQPGmean, na.rm=T), max( Sdat$AllGmean,  na.rm=T))
  plot(x=Sdat$Year, y=Sdat$AllGmean, col=ColsTrans[1], type="o", pch=19, ylim=ylim, xlab="Year", ylab="Escapement", lwd=1.4, cex=0.7)
  lines(x=Sdat$Year, y=Sdat$WildGmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[2])
  lines(x=Sdat$Year, y=Sdat$NoQPGmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[3])
  #legend is big, only put on one plot per page
  if(i %in% c(1,4,5,7)){
    legend( loc[lind[i]], col=ColsTrans[1:3], lty=rep(1,3), pch=rep(19,3), legend=c("All Data", 
           "Wild Data", "Wild No Q,LQ,P"), 
            cex=0.8, bty="n")
  }
  mtext(side=3, text=as.character(CUdf$CU_names[match(CU_short[i], CUdf$CU_short)]))
}

dev.off()



#############################################################################

# Same but by GU's ( Not for retrospective analysis, at P.V.W.'s request)

#need to group by site name to collapse any brood/rack/enhanced/rack etc categories
#only want to collapse year columns
AllDatBySite <- RawDat  %>% group_by(CU_Name, GroupName, GU_Name, NME, SummerRun, Rabcode)   %>%  summarise_each(funs(sum(., na.rm=T)), 9:69)
#need to turn all 0 Back to NA
AllDatBySite[AllDatBySite==0] <- NA
#need to change to character to preserve names
AllDatBySite$NME <- as.character(AllDatBySite$NME)
AllDatBySite$CU_Name <- as.character(AllDatBySite$CU_Name)
AllDatBySite$GroupName <- as.character(AllDatBySite$GroupName)
AllDatBySite$GU_Name <- as.character(AllDatBySite$GU_Name)
#SummerRun not registering True false arghhhhh 
AllDatBySite$SummerRun[which(is.na(AllDatBySite$SummerRun))] <- FALSE
# Remove summer run fish
AllDatBySite2 <- AllDatBySite[which(!AllDatBySite$SummerRun), ]

# Change to long form from wide form
LongDat <- AllDatBySite2 %>% gather( "Year", "Escape", 7:67)

# Infill 
GUEsc <- Infill(data=LongDat,  groupby=c("GU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME" )
# write to csv's
write.csv(GUEsc[[1]], "DataOut/InfilledByGUAll9_16.csv")
write.csv(GUEsc[[2]], "DataOut/InfilledEscByGU_9_16.csv")

# haven't done yet
pdf("Figures/InFilledDataByGU_9_16.pdf")

par(mfrow=c(2,1), mar=c(3,3,2,2), oma=c(2,2,2,2), mgp=c(2,0.8,0))
loc=c("bottomright", "bottomleft")
lind <- c(1,2,1,1,1,2,2)
for(i in 1:length(unique())){
  Sdat <- Infilled[which(Infilled$CU==CU_short[i]),]
  ylim=c(max(Sdat$Escape, Sdat$CuEsc, na.rm=T), min(Sdat$Escape, Sdat$CuEsc, na.rm=T))
  
  plot(x=Sdat$Year, y=Sdat$Escape, col="grey", type="o", pch=19, ylim=ylim, xlab="Year", ylab="Escapement", lwd=1.4, cex=0.7)
  lines(x=Sdat$Year, y=Sdat$CuEsc, type="o", pch=19, lwd=1.4, cex=0.7)
  legend( loc[lind[i]], col=c("grey", "black"), lty=rep(1,2), pch=rep(19,2), legend=c("Report Data", "Infilled Data"), 
          cex=0.9, bty="n")
  mtext(side=3, text=as.character(CUdf$CU_names[match(CU_short[i], CUdf$CU_short)]))
}

dev.off()







