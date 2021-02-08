# Infill Inside South Coast Chum escapement data, to supply Pieter Van Will
# for him to do run reconstruction separately (using catch data and proportion wild), 
# to supply the input returns data for retrospective analysis

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
library(readxl)
library(ggplot2)

# source functions for infilling, etc.
source("R/chumDataFunctions.r")

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
#
# ------------------#
# Further data notes - from Pieter Van Will
# ------------------#
#   Note that NAs are unobserved creeks.
#
#  	Regarding blank escapement count values vs. 0 values, specifically where wild was NA but brood was 0:
#     Pieter Van Will: I would consider any of those as Blanks and would require infill 
#     (perhaps if the Escapement is 0, convert that to a blank and I doubt there any real 0 
#     escapement numbers, even if there are some fish or a 0 in brood stock.)
#
#   Little Qualicum, Qualicum, and Puntledge Rivers: assume these are 100% enhanced, and 
#     remove entirely for infilling for wild escapement

# ----------------------------------------------------#
# Infill total spawners (wild, enhanced, brood, rack)
# ----------------------------------------------------#
#need to match CU names
CU_raw <- unique(rawdat$CU_Name)
CU_short <- c("SCS", "NEVI", "UK", "LB", "BI", "GS", "HSBI")
CU_names<-c("Southern Coastal Streams", "North East Vancouver Island", "Upper Knight",
            "Loughborough", "Bute Inlet", "Georgia Strait", "Howe Sound to Burrard Inlet" )
CUdf <- data.frame(CU_short, "CU_raw"=CU_raw[1:7], CU_names)


rawdat_f <- rawdat[rawdat$SummerRun==FALSE, ] # Remove summer run fish - earlier run timing means they are not intercepted in same fisheries, can't do run reconstruction with them

## Infill total spawners (wild/enhanced/brood/rack)
# wide to long format. maintain NA values (uncounted streams)
ldat <- rawdat_f %>% pivot_longer(cols=grep("[[:digit:]]{4}", names(rawdat_f)), names_to="Year", values_to="Escape")

# summarise by stream, to collapse any brood/rack/enhanced/rack etc categories
ldat_s <- ldat %>% group_by(CU_Name, GroupName, GU_Name, NME, SummerRun, Rabcode, Area, Year) %>%
  summarise(Escape=sum(Escape, na.rm=TRUE))
# # for sum escapement values that are 0, make back into NA (unobserved)
# since na.rm was TRUE, these became 0 values; summing NA values with na.rm=TRUE returns 0
ldat_s$Escape[ ldat_s$Escape == 0 ] <- NA


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

## Also want all infilled by site
# join non-fraser by-stream infilled data with by-CU infilled totals
xx <- left_join(AllDatAll, data.frame(CU_Name=AllByCU[[1]]$CU_Name, Year=AllByCU[[1]]$Year, CUEsc = AllByCU[[1]]$SiteEsc))
# For stream x CU combinations without any observations in a given year, estimate 
#     stream escapement from CU-level infilled escapement using stream proportion from
#     observed years and by-CU infilling escapement
xx$Escape <- ifelse(is.nan(xx$SiteEsc), xx$CUEsc*xx$Props, xx$SiteEsc)
#xx$CU <- CUdf$CU_short[match(xx$CU_Name, CUdf$CU_raw)] # don't need CU abbreviation
# Write total escapement by stream to csv, infilled by stream and CU
write.csv(xx, "DataOut/total_spawners_infilled_by_site.csv", row.names=F)

#Also want sumamrized by CU
write.csv(AllByCU[[1]], "DataOut/total_spawners_infilled_by_CU.csv")

# And sumamrized by CU x area
AllBySite <- xx
AllByArea <- AllBySite %>% group_by(CU_Name, Area, Year) %>% summarize( AreaTot = sum(Escape))
write.csv(AllByArea, "DataOut/total_spawners_infilled_by_area.csv")

# ----------------------------------------------------#
# Infill wild spawners only, and remove Little Qualicum, 
#               Qualicum, and Puntledge Rivers entirely
# ----------------------------------------------------#


# Remove non-wild and Little Qualicum, Qualicum, and Puntledge Rivers entirely
LongDatNoQP <- ldat[which(ldat$Source =="Wild" & 
                      ldat$NME %notin% c("QUALICUM RIVER", "LITTLE QUALICUM RIVER", "PUNTLEDGE RIVER")),]

#Now Infill
NoQPDat <- Infill(data = LongDatNoQP, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME")
#AllData
#write.csv(NoQPDat[[1]], "DataOut/InfilledNoQP_all.csv")
#Sumamarised data
#write.csv(NoQPDat[[2]], "DataOut/InfilledNoQP.csv")

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
write.csv(zz, "DataOut/wild_spawners_infilled_by_site.csv")

# Also by CU
write.csv(NoQPByCU[[1]], "DataOut/wild_spawners_infilled_by_CU.csv")

#Now by Area
AllBySiteWild <- zz
AllByAreaWild <- AllBySiteWild %>% group_by(CU_Name, Area, Year) %>% summarize( AreaTot = sum(Escape))
write.csv(AllByAreaWild, "DataOut/wild_spawners_infilled_by_area.csv")

#===============================================#
#   Check against infilling of 2013 data =======
#===============================================#

# 2018 infilled data by site
total_spawners_infilled_by_site_2018 <- AllBySite
wild_spawners_infilled_by_site_2018 <- AllBySiteWild

# 2013 data infilled by site
total_spawners_infilled_by_site_2013 <- read_excel("DataIn/2013_infilled_data_PieterVanWill_for_checking.xlsx", sheet=1)
wild_spawners_infilled_by_site_2013 <- read_excel("DataIn/2013_infilled_data_PieterVanWill_for_checking.xlsx", sheet=2)

# merge totals 2013/2018 data 
total_infill <- merge(total_spawners_infilled_by_site_2013, total_spawners_infilled_by_site_2018, 
                      by=c("CU_Name", "NME", "Year", "Area", "Rabcode"), suffixes=c("_2013d", "_2018d"), all=TRUE)
# merge wild 2013/2018 data
wild_infill <- merge(wild_spawners_infilled_by_site_2013, wild_spawners_infilled_by_site_2018, 
                      by=c("CU_Name", "NME", "Year", "Area", "Rabcode"), suffixes=c("_2013d", "_2018d"), all=TRUE)
# compare for totals
ggplot(total_infill, aes(x=Escape_2013d, y=Escape_2018d, colour=Year)) +
  geom_point() +
  facet_wrap(~CU_Name, scales="free") + 
  geom_abline(aes(intercept=0, slope=1))
# compare for wild
ggplot(wild_infill, aes(x=Escape_2013d, y=Escape_2018d, colour=Year)) +
  geom_point() +
  facet_wrap(~CU_Name, scales="free") + 
  geom_abline(aes(intercept=0, slope=1))
# Something going weird for southern coastal streams

# Note that everything seems to match up okay except Southern Coastal Stream. 
# It appears that Viner Sound Creek is classified as fall run in the 2013 data and summer run in 2018.
# Since this is a larger run, it changes the proportions of the infilling, leading to pretty different
# infilling values between the 2013 data and the updated 2018 data.
# From Pieter Van Will: 
# Viner should be assigned as Summer.  It is a later timing than most of the summer runs but would not
# really be impacted by the typical fall chum targeted fisheries.  It should not be included in this analysis.

#Compare infilling by stream, 2013 and 2018 data
CUs <- unique(total_infill$CU_Name)
compare_stream_infilling_fig <- function(CU, dat) {
  d <- dat
  d <- d[d$CU_Name==CU, ]
  ggplot(d, aes(y=Escape_2013d, x=Year)) +
    geom_point() +
    geom_point(aes(y=Escape_2018d, x=Year), colour="dodgerblue", shape=1, stroke=2) +
    facet_wrap(~NME, scales="free_y") +
    scale_x_discrete(breaks=seq(1960,2020,10)) +
    ggtitle(CU) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))
}
# for total spawners
fig_list <- as.list(CUs)
fig_list <- purrr::map(fig_list, compare_stream_infilling_fig, dat=total_infill)
fig_list[[1]] # check
pdf("Figures/fig_check_infill_by_stream_total_spawners.pdf", width=18, height=12, pointsize=8)
for (i in 1:length(fig_list)){
  print(fig_list[[i]])
}
dev.off()
# For wild spawners
fig_list_wild <- as.list(CUs)
fig_list_wild <- purrr::map(fig_list_wild, compare_stream_infilling_fig, dat=wild_infill)
fig_list_wild[[1]] # check
pdf("Figures/fig_check_infill_by_stream_wild_spawners.pdf", width=18, height=12, pointsize=8)
for (i in 1:length(fig_list_wild)){
  print(fig_list_wild[[i]])
}
dev.off()

# Compare infilling by CU, 2013 and 2018 data
total_spawners_infilled_by_CU_2013
total_spawners_infilled_by_CU_2018 


