# Infill Inside South Coast Chum escapement data, to supply Pieter Van Will
# for him to do run reconstruction separately (using catch data and proportion wild)

# This code infills at multiple levels 
# Total Spawners - by CU, area and stream
# Wild Spawners - by CU, area and stream

# Note: The following infilling & brood table reconstruction procedures are taken from: 
# Holt, C.A., Davis, B., Dobson, D., Godbout, L., Luedke, W., Tadey, J., and Van Will, P. 2018.
# Evaluating Benchmarks of Biological Status for Data-limited Conservation Units of Pacific
# Salmon, Focusing on Chum Salmon in Southern BC. DFO Can. Sci. Advis. Sec. Res. Doc.
# 2018/011. ix + 77 p. Available at: https://waves-vagues.dfo-mpo.gc.ca/Library/40759386.pdf

# All code in this section (including code in chumDataFunctions.r) was written by B. Davis (DFO) for the  
# above paper, and provided to us by Carrie Holt in October 2020 as part of "Retrospective Analysis BD" folder.

# Some modifications by Luke Warkentin (DFO) in 2020-2021

setwd("C:/github/SalmonLRP_RetroEval/SCChumStudy")

library(dplyr)
library(tidyr)
library(ggplot2)

# source functions for infilling, etc.
source("chumDataFunctions.r")

# read in raw escapement data, with years up to 2018, from Pieter Van Will
rawdat <- read.csv("DataIn/Chum Escapement Data With Areas_2018.csv", stringsAsFactors = FALSE, check.names=FALSE, strip.white = TRUE) # strip.white for leading and trailing white spaces in Source and SummerRun columns

# ----------------------------------------------------#
# Notes on data variables (From Pieter Van Will): 
# ----------------------------------------------------#
#
#   Source:
#     -Wild are wild spawners
#     -Rack is the # of fish that are harvested at a facility and do not contribute to the system 
#           (part of the total return but not contributing as spawners or enhanced)
#     -Enhanced was an assignment we gave given the typically high magnitude of the enhancement, if it was 
#           a large facility on the system with a large chum production we tended to call the system as 
#           a whole enhanced (like Puntledge or the systems with large spawning channels like Qualicum)
#     -Brood are the fish from the system that are taken to support the supplementation of that stock 
#           (do contribute to the subsequent returns).
#
#   SummerRun: We removed the summer run fish as all the data in regards the reconstruction work is 
#       associated with Fall timed stocks.  
#
#   NME: Individual streams. Note for Qualicum River, Little Qualicum River, Puntledge River - Historically we assume 
#       these three stocks are 100% enhanced at least since enhancement began at those locations.  
#       We have little data in the enhanced contribution found in the returns but for the purposes of 
#       pulling out wild we make the assumption they were 100% enhanced and not included. 



# ----------------------------------------------------#
# Infill total spawners (wild, enhanced, brood, rack)
# ----------------------------------------------------#
#need to match CU names
CU_raw <- unique(rawdat$CU_Name)
CU_short <- c("SCS", "NEVI", "UK", "LB", "BI", "GS", "HSBI")
CU_names<-c("Southern Coastal Streams", "North East Vancouver Island", "Upper Knight",
            "Loughborough", "Bute Inlet", "Georgia Strait", "Howe Sound to Burrard Inlet" )
CUdf <- data.frame(CU_short, "CU_raw"=CU_raw[1:7], CU_names)


rawdat_f <- rawdat[rawdat$SummerRun==FALSE, ] # Remove summer run fish

## All Data -- wild/enhanced/brood/rack
# wide to long format. maintain NA values (uncounted streams)
ldat <- rawdat_f %>% pivot_longer(cols=grep("[[:digit:]]{4}", names(rawdat_f)), names_to="Year", values_to="Escape")

# Note that NAs are unobserved creeks.
# add column that is TRUE if source is wild and escapement is actually a 0 value (as opposed to NA)
ldat$true_0_wild <- ifelse(ldat$Source=="Wild" & ldat$Escape==0, TRUE, FALSE)
ldat$true_0_wild[is.na(ldat$true_0_wild)] <- FALSE # make NA values FALSE

# summarise by stream, to collapse any brood/rack/enhanced/rack etc categories
ldat_s <- ldat %>% group_by(CU_Name, GroupName, GU_Name, NME, SummerRun, Rabcode, Area, Year, true_0_wild) %>%
  summarise(Escape=sum(Escape, na.rm=TRUE))
# for sum escapement values that are 0 and not true 0 wild counts, make back into NA
ldat_s$Escape[ ldat_s$Escape == 0 & ldat_s$true_0_wild == FALSE ] <- NA

# Figure out what proportion of sites are actually observed for each CU/year
ObsDF <- data.frame(CU=character(), Year=numeric(), Nobs=numeric(), Ninfilled=numeric())
CUs <- unique(ldat_s$CU_Name)
for(i in 1:length(CUs)){
  CUdat <-ldat_s[which(ldat_s$CU_Name==CUs[i]),]
  years <- unique(CUdat$Year)
  for(j in 1:length(years)){
    Ydat <- CUdat[which(CUdat$Year==years[j]),]
    Nobs <- sum(is.na(Ydat$Escape)==F)
    Ninfilled <- sum(is.na(Ydat$Escape))
    ObsDF <- rbind(ObsDF, data.frame(CU=CUs[i], Year=years[j], Nobs, Ninfilled))
  }
}
write.csv(ObsDF, "DataOut/summary_n_infilled_by_year.csv")

#Now Infill for total escapement, no summer run
AllDat <- Infill(data = ldat_s, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME")
# The output is a list of two elements. 
# Element 1: infilled total escapement, by stream
# Element 2: results of Element 1, grouped by CU and year and summed (raw escapement and infilled escapement)
#     For example, for CUs that were not counted in a year, GroupEsc is NA and SumRawEsc = 0.
AllDat_el1 <- AllDat[[1]]
AllDat_el2 <- AllDat[[2]]

# # check infilling
# temp <- AllDat[[2]]
# ggplot(temp, aes(y=GroupEsc, x=Year, group=CU_Name)) +
#   geom_point( colour="dodgerblue") +
#   geom_path( colour="dodgerblue") +
#   geom_point(aes(y=SumRawEsc, x=Year)) +
#   facet_wrap(~CU_Name, scales="free_y") +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle=90, vjust=0.5))

# Remove Fraser, make data frame to feed into infill again
# This is the by-stream infilling, not summarized, no Fraser CUs
AllDatAll <- AllDat[[1]][which(AllDat[[1]]$CU_Name %in% CUdf$CU_raw),c("CU_Name", "NME", "Year", "Props", "SiteEsc", "Area", "Rabcode")]

# this is the by-stream infilling, summarized by CU, no Fraser CUs
AllDatSumm <- as.data.frame(AllDat[[2]][which(AllDat[[2]]$CU_Name %in% CUdf$CU_raw),])

# Now infill missing years for entire CUs for years when there were no visits to streams in that CU

# CU-level infilling: Infill missing values for total escapement by CU (for CUs that did not get counts in certain years)
AllByCU <- Infill(data=AllDatSumm, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc")

AllByCU_el1 <- AllByCU[[1]]
AllByCU_el2 <- AllByCU[[2]]


## Also want all infilled by site
# join non-fraser by-stream infilled data with by-CU infilled totals
xx <- left_join(AllDatAll, data.frame(CU_Name=AllByCU[[1]]$CU_Name, Year=AllByCU[[1]]$Year, CUEsc = AllByCU[[1]]$SiteEsc))
# For stream x CU combinations without any observations in a given year, estimate 
#     stream escapement from CU-level infilled escapement using stream proportion from
#     observed years and by-CU infilling escapement
xx$Escape <- ifelse(is.nan(xx$SiteEsc), xx$CUEsc*xx$Props, xx$SiteEsc)
#xx$CU <- CUdf$CU_short[match(xx$CU_Name, CUdf$CU_raw)] # don't need CU abbreviation
# Write total escapement by stream to csv, infilled by stream and CU
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
#write.csv(DFout, "DataOut/InfilledAllDat2_15.csv")

# dumb work around, but need years to not be factors to link data sets
AllByCUSumm$Year <- as.character(AllByCUSumm$Year)
AllByCUSumm$Year <- as.numeric(AllByCUSumm$Year)

#now merge two together, need to remove "Escape" Column from output
#Infilled <- left_join(AllByCUSumm[,-which(names(AllByCUSumm)=="Escape")],SRdat[,1:4], by=c("CU", "Year"))
# 
# pdf("Figures/InFilled_2_26_All.pdf")
# 
# par(mfrow=c(2,1), mar=c(3,3,2,2), oma=c(2,2,2,2), mgp=c(2,0.8,0))
# loc=c("topright", "topleft")
# lind <- c(1,2,1,1,1,2,2)
# for(i in 1:length(unique(CU_short))){
#   Sdat <- Infilled[which(Infilled$CU==CU_short[i]),]
#   ylim=c(min(Sdat$Escape, Sdat$CuEsc, na.rm=T),max(Sdat$Escape, Sdat$SiteEsc, na.rm=T))
#   
#   plot(x=Sdat$Year, y=Sdat$Escape, col="grey", type="o", pch=19, ylim=ylim, xlab="Year", ylab="Escapement", lwd=1.4, cex=0.7)
#   lines(x=Sdat$Year, y=Sdat$SiteEsc, type="o", pch=19, lwd=1.4, cex=0.7)
#   legend( loc[lind[i]], col=c("grey", "black"), lty=rep(1,2), pch=rep(19,2), legend=c("Report Data", "Infilled Data"), 
#           cex=0.9, bty="n")
#   mtext(side=3, text=as.character(CUdf$CU_names[match(CU_short[i], CUdf$CU_short)]))
# }
# 
# dev.off()

#now my data mostly matches PVW's report, except NEVI 1981, and missing values

# ########################################################################
# # What if use arithmatic mean instead of geometric??
# 
# AmeanAll <- Infill(data = LongDat, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME", avg="Amean")
# 
# # Remove Fraser, make data frame to feed into infill again
# AmeanBySiteAll <- as.data.frame(AmeanAll[[2]][which(AmeanAll[[2]]$CU_Name %in% CUdf$CU_raw),])
# 
# # Now infill missing years for entire CU
# AmeanAllByCU <- Infill(data=AmeanBySiteAll, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc", avg="Amean")
# 
# #Join with report data
# AmeanSumm <- AmeanAllByCU[[1]]
# AmeanSumm$CU <- CUdf$CU_short[match(AmeanSumm$CU_Name, CUdf$CU_raw)]
# 
# # dumb work around, but need years to not be factors to link data sets
# AmeanSumm$Year <- as.character(AmeanSumm$Year)
# AmeanSumm$Year <- as.numeric(AmeanSumm$Year)
# 
# #now merge two together, need to remove "Escape" Column from output
# AmeanInfilled <- left_join(AmeanSumm[,-which(names(AmeanSumm)=="Escape")],SRdat[,1:4], by=c("CU", "Year"))
# 
# 
# #Now plot report, GeoMean, Amean
# 
# pdf("Figures/InFilled9_17_CompMethod.pdf")
# 
# par(mfrow=c(2,1), mar=c(3,3,2,2), oma=c(2,2,2,2), mgp=c(2,0.8,0))
# loc=c("topright", "topleft")
# lind <- c(1,2,1,1,1,2,2)
# for(i in 1:length(unique(CU_short))){
#   Sdat <- Infilled[which(Infilled$CU==CU_short[i]),]
#   AmeanSdat <- AmeanInfilled[which(AmeanInfilled$CU==CU_short[i]),]
#   ylim=c( min(Sdat$Escape, Sdat$CuEsc, na.rm=T), max(Sdat$Escape, Sdat$SiteEsc, AmeanSdat$SiteEsc, na.rm=T))
#   
#   plot(x=Sdat$Year, y=Sdat$Escape, col="grey", type="o", pch=19, ylim=ylim, xlab="Year", ylab="Escapement", lwd=1.4, cex=0.7)
#   lines(x=Sdat$Year, y=Sdat$SiteEsc, type="o", pch=19, lwd=1.4, cex=0.7, col="blue")
#   lines(x=AmeanSdat$Year, y=AmeanSdat$SiteEsc, type="o", pch=19, lwd=1.4, cex=0.7, col="red")
#   legend( loc[lind[i]], col=c("grey", "blue", "red"), lty=rep(1,3), pch=rep(19,3), 
#           legend=c("Report Data", "Infilled Geo. Mean", "Infilled Arith. Mean"), 
#           cex=0.9, bty="n")
#   mtext(side=3, text=as.character(CUdf$CU_names[match(CU_short[i], CUdf$CU_short)]))
# }
# 
# dev.off()


################################################################################################################

# Now start from beginning and take out all non-wild, re-do infilling

# Remove non-wild

# Change to long form from wide form
LongDatWild <- ldat[ldat$Source=="Wild",]

#Now Infill
WildDat <- Infill(data = LongDatWild, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME")
#WildDatAmean <- Infill(data = LongDatWild, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME", avg="Amean")
#AllData
#write.csv(WildDat[[1]], "DataOut/InfilledWild_all.csv")
#Sumamarised data
# write.csv(WildDat[[2]], "DataOut/InfilledWild.csv")

# Remove Fraser, make data frame to feed into infill again
WildDatSumm <- as.data.frame(WildDat[[2]][which(WildDat[[2]]$CU_Name %in% CUdf$CU_raw),])
#WildDatSummAmean <- as.data.frame(WildDatAmean[[2]][which(WildDatAmean[[2]]$CU_Name %in% CUdf$CU_raw),])
WildDatAll <- WildDat[[1]][which(WildDat[[1]]$CU_Name %in% CUdf$CU_raw),c("CU_Name", "NME", "Year", "Props", "SiteEsc", "Area", "Rabcode")]


# Now infill missing years for entire CU

#Infill missing values
WildByCU <- Infill(data=WildDatSumm, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc")
# also compare arithmatic mean
#WildByCUAmean <-  Infill(data=WildDatSummAmean, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc", avg="Amean")

## Also want all infilled by site
yy <- left_join(WildDatAll, data.frame(CU_Name=WildByCU[[1]]$CU_Name, Year=WildByCU[[1]]$Year, CUEsc = WildByCU[[1]]$SiteEsc))
yy$Escape <- ifelse(is.nan(yy$SiteEsc), yy$CUEsc*yy$Props, yy$SiteEsc)
#yy$CU <- CUdf$CU_short[match(yy$CU_Name, CUdf$CU_raw)]
write.csv(yy, "DataOut/WildInfilledBySite.csv")


#Add CU column
WildSumm <- WildByCU[[1]]
WildSumm$CU <- CUdf$CU_short[match(WildSumm$CU_Name, CUdf$CU_raw)]
WildSumm$Year <- as.character(WildSumm$Year)
WildSumm$Year <- as.numeric(WildSumm$Year)
# WildSummAmean <- WildByCUAmean[[1]]
# WildSummAmean$CU <- CUdf$CU_short[match(WildSummAmean$CU_Name, CUdf$CU_raw)]
# WildSummAmean$Year <- as.character(WildSummAmean$Year)
# WildSummAmean$Year <- as.numeric(WildSummAmean$Year)

#make simple data frame to export with geo. mean infilling
# WildDFout <- data.frame(Year=WildSumm$Year, CU=WildSumm$CU, Escape=WildSumm$SiteEsc)
# write.csv(WildDFout, "DataOut/InfilledWild.csv")


# #Join All Infilled
# # Start with "Infilled"
# InfillComp <- data.frame("Year"=Infilled$Year, "CU"=Infilled$CU, "ReportEsc"=Infilled$Escape, "AllGmean" =Infilled$SiteEsc)
# # Now add Amean data
# InfillComp <- left_join(InfillComp, data.frame("Year"=AmeanSumm$Year, "CU"=AmeanSumm$CU, "AllAmean"=AmeanSumm$SiteEsc), by=c("Year", "CU"))
# #Now add Wild Gmean
# InfillComp <- left_join(InfillComp, data.frame("Year"=WildSumm$Year, "CU"=WildSumm$CU, "WildGmean"=WildSumm$SiteEsc), by=c("Year", "CU"))
# # Wild Amean
# InfillComp <- left_join(InfillComp, data.frame("Year"=WildSummAmean$Year, "CU"=WildSummAmean$CU, "WildAmean"=WildSummAmean$SiteEsc), by=c("Year", "CU"))
# 
# Cols <- c("#17AB69", "#2120B7",  "#FF314B", "#E8B338")
# Trans <- 80
# ColsTrans <- paste(Cols, Trans, sep="")
# 
# pdf("Figures/InFilled_Comp_All.pdf")
# 
# par(mfrow=c(2,1), mar=c(3,3,2,2), oma=c(2,2,2,2), mgp=c(2,0.8,0))
# loc=c("topright", "topleft")
# lind <- c(1,2,1,1,1,2,2)
# for(i in 1:length(unique(CU_short))){
#   Sdat <- InfillComp[which(InfillComp$CU==CU_short[i]),]
#   ylim=c(min(Sdat$WildGmean, Sdat$WildAmean, na.rm=T), max(Sdat$ReportEsc, Sdat$AllGmean, Sdat$AllAmean,  na.rm=T))
#   plot(x=Sdat$Year, y=Sdat$ReportEsc, col="grey", type="o", pch=19, ylim=ylim, xlab="Year", ylab="Escapement", lwd=1.4, cex=0.7)
#   lines(x=Sdat$Year, y=Sdat$AllGmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[1])
#   lines(x=Sdat$Year, y=Sdat$AllAmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[2])
#   lines(x=Sdat$Year, y=Sdat$WildGmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[3])
#   lines(x=Sdat$Year, y=Sdat$WildAmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[4])
#   #legend is big, only put on one plot per page
#   if(i %in% c(1,4,5,7)){
#     legend( loc[lind[i]], col=c("grey", ColsTrans), lty=rep(1,5), pch=rep(19,5), legend=c("Report Data", "All Geo. Mean", 
#                                                                                           "All Arith. Mean", "Wild Geo. Mean", "Wild Arith. Mean"), 
#             cex=0.8, bty="n")
#   }
#   mtext(side=3, text=as.character(CUdf$CU_names[match(CU_short[i], CUdf$CU_short)]))
# }
# 
# dev.off()
# 

################################################################################
#Same as above but also remove Little Qualicum, Qualicum, and Puntledge Rivers entirely

# Remove non-wild and Little Qualicum, Qualicum, and Puntledge Rivers entirely
LongDatNoQP <- ldat[which(ldat$Source =="Wild" & 
                      ldat$NME %notin% c("QUALICUM RIVER", "LITTLE QUALICUM RIVER", "PUNTLEDGE RIVER")),]

#Now Infill
NoQPDat <- Infill(data = LongDatNoQP, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME")
#AllData
write.csv(NoQPDat[[1]], "DataOut/InfilledNoQP_all.csv")
#Sumamarised data
write.csv(NoQPDat[[2]], "DataOut/InfilledNoQP.csv")

# Remove Fraser, make data frame to feed into infill again
NoQPDatSumm <- as.data.frame(NoQPDat[[2]][which(NoQPDat[[2]]$CU_Name %in% CUdf$CU_raw),])
NoQPDatAll <- NoQPDat[[1]][which(NoQPDat[[1]]$CU_Name %in% CUdf$CU_raw),c("CU_Name", "NME", "Year", "Props", "SiteEsc", "Area", "Rabcode")]

# Now infill missing years for entire CU

#Infill missing values
NoQPByCU <- Infill(data=NoQPDatSumm, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc")
#NoQPByCUAmean <-  Infill(data=NoQPDatSummAmean, groupby=NULL, Uid = NULL , unit="CU_Name", EscCol="GroupEsc", avg="Amean")

## Also want all infilled by site
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
# NoQPSummAmean <- NoQPByCUAmean[[1]]
# NoQPSummAmean$CU <- CUdf$CU_short[match(NoQPSummAmean$CU_Name, CUdf$CU_raw)]
# NoQPSummAmean$Year <- as.character(NoQPSummAmean$Year)
# NoQPSummAmean$Year <- as.numeric(NoQPSummAmean$Year) 

#simple data frame for export
# NoQPDFout <- data.frame(Year=NoQPSumm$Year, CU=NoQPSumm$CU, Escape=NoQPSumm$SiteEsc)
# write.csv(NoQPDFout, "DataOut/InfilledWildNoQP.csv")

#make plot showing all, wild, this
#Join data together
# Start with "Infilled"
InfillWild <- data.frame("Year"=Infilled$Year, "CU"=Infilled$CU, "ReportEsc"=Infilled$Escape, "AllGmean" =Infilled$SiteEsc)
#Now add Wild Gmean
InfillWild <- left_join(InfillWild, data.frame("Year"=WildSumm$Year, "CU"=WildSumm$CU, "WildGmean"=WildSumm$SiteEsc), by=c("Year", "CU"))
# Wild No QP
InfillWild <- left_join(InfillWild, data.frame("Year"=NoQPSumm$Year, "CU"=NoQPSumm$CU, "NoQPGmean"=NoQPSumm$SiteEsc), by=c("Year", "CU"))

# 
# Cols <- c("#17AB69", "#2120B7",  "#FF314B", "#E8B338")
# Trans <- 80
# ColsTrans <- paste(Cols, Trans, sep="")
# 
# pdf("Figures/InFilled_CompWild.pdf")
# 
# par(mfrow=c(2,1), mar=c(3,3,2,2), oma=c(2,2,2,2), mgp=c(2,0.8,0))
# loc=c("topright", "topleft")
# lind <- c(1,2,1,1,1,2,2)
# for(i in 1:length(unique(CU_short))){
#   Sdat <- InfillWild[which(InfillWild$CU==CU_short[i]),]
#   ylim=c(min(Sdat$NoQPGmean, na.rm=T), max( Sdat$AllGmean,  na.rm=T))
#   plot(x=Sdat$Year, y=Sdat$AllGmean, col=ColsTrans[1], type="o", pch=19, ylim=ylim, xlab="Year", ylab="Escapement", lwd=1.4, cex=0.7)
#   lines(x=Sdat$Year, y=Sdat$WildGmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[2])
#   lines(x=Sdat$Year, y=Sdat$NoQPGmean, type="o", pch=19, lwd=1.4, cex=0.7, col=ColsTrans[3])
#   #legend is big, only put on one plot per page
#   if(i %in% c(1,4,5,7)){
#     legend( loc[lind[i]], col=ColsTrans[1:3], lty=rep(1,3), pch=rep(19,3), legend=c("All Data", 
#                                                                                     "Wild Data", "Wild No Q,LQ,P"), 
#             cex=0.8, bty="n")
#   }
#   mtext(side=3, text=as.character(CUdf$CU_names[match(CU_short[i], CUdf$CU_short)]))
# }
# 
# dev.off()

# 
# 
# #############################################################################
# 
# # Same but by GU's ( Not for retrospective analysis, at P.V.W.'s request)
# 
# #need to group by site name to collapse any brood/rack/enhanced/rack etc categories
# #only want to collapse year columns
# AllDatBySite <- RawDat  %>% group_by(CU_Name, GroupName, GU_Name, NME, SummerRun, Rabcode)   %>%  summarise_each(funs(sum(., na.rm=T)), 9:69)
# #need to turn all 0 Back to NA
# AllDatBySite[AllDatBySite==0] <- NA
# #need to change to character to preserve names
# AllDatBySite$NME <- as.character(AllDatBySite$NME)
# AllDatBySite$CU_Name <- as.character(AllDatBySite$CU_Name)
# AllDatBySite$GroupName <- as.character(AllDatBySite$GroupName)
# AllDatBySite$GU_Name <- as.character(AllDatBySite$GU_Name)
# #SummerRun not registering True false arghhhhh 
# AllDatBySite$SummerRun[which(is.na(AllDatBySite$SummerRun))] <- FALSE
# # Remove summer run fish
# AllDatBySite2 <- AllDatBySite[which(!AllDatBySite$SummerRun), ]
# 
# # Change to long form from wide form
# LongDat <- AllDatBySite2 %>% gather( "Year", "Escape", 7:67)
# 
# # Infill 
# GUEsc <- Infill(data=LongDat,  groupby=c("GU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME" )
# # write to csv's
# write.csv(GUEsc[[1]], "DataOut/InfilledByGUAll9_16.csv")
# write.csv(GUEsc[[2]], "DataOut/InfilledEscByGU_9_16.csv")
# 
# # haven't done yet
# pdf("Figures/InFilledDataByGU_9_16.pdf")
# 
# par(mfrow=c(2,1), mar=c(3,3,2,2), oma=c(2,2,2,2), mgp=c(2,0.8,0))
# loc=c("bottomright", "bottomleft")
# lind <- c(1,2,1,1,1,2,2)
# for(i in 1:length(unique())){
#   Sdat <- Infilled[which(Infilled$CU==CU_short[i]),]
#   ylim=c(max(Sdat$Escape, Sdat$CuEsc, na.rm=T), min(Sdat$Escape, Sdat$CuEsc, na.rm=T))
#   
#   plot(x=Sdat$Year, y=Sdat$Escape, col="grey", type="o", pch=19, ylim=ylim, xlab="Year", ylab="Escapement", lwd=1.4, cex=0.7)
#   lines(x=Sdat$Year, y=Sdat$CuEsc, type="o", pch=19, lwd=1.4, cex=0.7)
#   legend( loc[lind[i]], col=c("grey", "black"), lty=rep(1,2), pch=rep(19,2), legend=c("Report Data", "Infilled Data"), 
#           cex=0.9, bty="n")
#   mtext(side=3, text=as.character(CUdf$CU_names[match(CU_short[i], CUdf$CU_short)]))
# }
# 
# dev.off()
# 
# 
