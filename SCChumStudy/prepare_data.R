# Read in raw Chum data, infill and prepare for LRP Analysis, save as csv files
# Plot data 

library(ggplot2)


# Only need if running independently from runSouthCoastChum.R
setwd("C:/github/SalmonLRP_RetroEval/SCChumStudy")

# ===========================================================================================
# Read-in Chum Data & Run Data Prep Functions:  
# ===========================================================================================

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

#Now Infill. Gives a list of 2 data frames. First is infilling by stream, second is a summary that is summarized by CU. 
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
# will not get recruit estimate for last 
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




# =====================================#
# Explore data with figures - LW
# =====================================#

# format data
yrcols <- grep("[[:digit:]]{4}", names(RawDat)) # get position of yr columns
rdl <- RawDat %>% pivot_longer( cols= yrcols, names_to="year", values_to="escapement") # wide to long format for plotting
#check_names <- rdl %>% group_by(CU_Name, GroupName, GU_Name) %>% summarise(n=n()) # check correspondence of CU and GroupName, GU_Name

rdl <- rdl[rdl$Source=="Wild", ] # remove non-wild fish
rdl <- rdl[!(rdl$NME %in% c("QUALICUM RIVER", "LITTLE QUALICUM RIVER", "PUNTLEDGE RIVER")), ] # FLAG - why? remove three rivers removed from analysis below
rdl <- rdl[ rdl$SummerRun==FALSE, ] # remove summer run fish

#-----------#
# Look at raw escapement data - LW
#-----------#
# Appears that GroupName and GU_Name are not nested within CUs, but different categories that overlap irregularly

# Check "QUALICUM RIVER", "LITTLE QUALICUM RIVER", "PUNTLEDGE RIVER"
rdl[rdl$NME %in% c("QUALICUM RIVER", "LITTLE QUALICUM RIVER", "PUNTLEDGE RIVER"),] %>%
  ggplot( aes(y=escapement, x=year, colour=Source, colour=NME, group=Source)) +
  geom_point() +
  geom_path() +
  scale_y_log10() +
  facet_wrap(~NME) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))

# Plot raw escapement, normal vs. summer run, 1 page per CU
options(scipen = 100000)
CUs <- unique(rdl$CU_Name)
make_CU_figs <- function(CU) {
  rdl[rdl$CU_Name==CU,] %>% ggplot( aes(y=escapement, x=year, colour=SummerRun, group=NME)) +
    geom_point() +
    geom_path() +
    scale_y_log10() +
    #facet_grid(NME~.) +
    scale_colour_manual(values=c("black", "red")) +
    ggtitle(CU) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))
}
fig_list <- as.list(CUs)
fig_list <- lapply(fig_list, make_CU_figs)
fig_list[[1]] # check
pdf("Figures/fig_raw_escapement_by_CU.pdf", width=7, height=4, pointsize=8)
for (i in 1:length(fig_list)){
  print(fig_list[[i]])
}
dev.off()


#------------#
# Look at by-stream infilled data - LW
#------------#

# Plot actual vs. infilled data by stream, one page per CU
infill_by_stream <- NoQPDat[[1]] # get by-steam infilled data

CUs <- unique(infill_by_stream$CU_Name)
make_CU_figs <- function(CU) {
  infill_by_stream[infill_by_stream$CU_Name==CU,] %>% ggplot( aes(y=SiteEsc, x=Year, group=NME)) +
    geom_point(colour="dodgerblue", shape=1) +
    #geom_path() +
    #geom_point(aes(y=ContrEsc, x=Year), colour="dodgerblue") +
    geom_point(aes(y=Escape, x=Year), colour="black", shape=1) +
    ggtitle(CU) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))
}
fig_list <- as.list(CUs)
fig_list <- lapply(fig_list, make_CU_figs)
fig_list[[1]] # check
pdf("Figures/fig_by-CU_infill_escapement.pdf", width=7, height=4, pointsize=8)
for (i in 1:length(fig_list)){
  print(fig_list[[i]])
}
dev.off()

# Plot observed and infilled escapement for each CU, one panel per stream
CUs <- unique(infill_by_stream$CU_Name)
make_stream_figs <- function(CU) {
  infill_by_stream[infill_by_stream$CU_Name==CU,] %>% 
    ggplot( aes(y=SiteEsc, x=Year, group=NME)) +
    geom_point(colour="dodgerblue") +
    geom_path() +
    geom_point(aes(y=Escape, x=Year), colour="black") +
    #scale_y_log10() +
    facet_wrap(~NME, scales="free_y") +
    scale_x_discrete(breaks=seq(1960,2010,10)) +
    ggtitle(CU) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))
}
fig_list <- as.list(CUs)
fig_list <- lapply(fig_list, make_stream_figs)
fig_list[[1]] # check
pdf("Figures/fig_by-stream_infill_escapement_streams.pdf", width=18, height=12, pointsize=8)
for (i in 1:length(fig_list)){
  print(fig_list[[i]])
}
dev.off()

# Plot correlations between actual stream escapements by CU. 
# Infilling by stream assumes escapements in streams within a CU have correlation. 
# change campbell river duplicate
rdl$NME[which(rdl$NME=="CAMPBELL RIVER" & rdl$GU_Name=="5 - Fraser")] <- "LITTLE CAMPBELL RIVER"
corCU <- function(CU) {
  df1 <- rdl %>% filter(CU_Name==CU) %>% select(NME, year, escapement) %>%
    pivot_wider(names_from=NME, values_from=escapement)
  cormat <- cor(df1[,-1], use="pairwise.complete.obs")
  cormat1 <- cormat[upper.tri(cormat, diag=FALSE)]
  hist(cormat1, main=CU, xlim=c(-1,1), xlab="Correlations between observed stream escapements")
  #PerformanceAnalytics::chart.Correlation(cor(df1[,-1], use="pairwise.complete.obs") )
}
png("Figures/fig_observed_spawners_by_stream_correlations.png", width=10, height=6, units="in", res=300)
layout(matrix(seq_along(CUs), ncol=3, byrow = TRUE))
for ( i in seq_along(CUs)) {
  corCU(CUs[i])
}
dev.off()

# ----------#
# Plot infilling by CU, to compare actual, infilled by stream, and infilling by CU - LW
# ----------#
infill_by_stream_no_fraser <- NoQPDat[[2]] %>% filter(!(CU_Name %in% c("8 - Lower Fraser", "9 - Fraser Canyon"))) 
png("Figures/fig_compare_actual_infill_by_stream_and_CU.png", height=5, width=10, res=300, units="in")
ggplot(data=NoQPByCU[[1]], aes(y= SiteEsc, x=Year )) + # infilled by CU
  geom_point( colour="red") + # infill by CU is in red
  geom_point(data = infill_by_stream_no_fraser, aes(y=GroupEsc, x=Year), colour="dodgerblue") + # infill by stream is blue
  geom_point(data = infill_by_stream_no_fraser, aes(y=SumRawEsc, x=Year)) + # raw escapement data
  geom_path(data = infill_by_stream_no_fraser, aes(y=SumRawEsc, x=Year, group=CU_Name)) +
  ylab("Escapement") +
  facet_wrap(~CU_Name, scales="free_y") +
  #scale_y_log10() +
  theme_classic() +
  scale_x_discrete(breaks=seq(1960,2010,10)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5))
dev.off()

# Plot escapement and R/S time series on same x axis for each CU --------
# Merge escapement and recruitment data by CU and year
cdat <- merge(NoQPByCU[[1]], Btable, by.x=c("CU_Name", "Year"), by.y=c("CU", "Year"), all=TRUE)
str(cdat)
cdat$Year <- as.numeric(cdat$Year)
cdat$RS <- cdat$Recruit/cdat$SiteEsc

# cdat$RS_pos <- ifelse(cdat$RS >1, 1, 2) # binary colours
# cols <- c("dodgerblue", "firebrick")
ncols <- 20
# get colours
colors <- paletteer::paletteer_c( palette = "scico::roma", n = ncols)
#colors <- paletteer::paletteer_c( palette = "pals::warmcool", n = ncols)

# Transform the numeric variable in bins
rank <- as.factor( as.numeric( cut(log(cdat$RS), ncols)))


options(scipen = 100000)

png("Figures/fig_escapement_RS.png", height=10, width=18, res=200, pointsize=20, units="in")
layout(mat=matrix(1:16, nrow=4, byrow=FALSE))
  yrs <- range(cdat$Year)
  CUs <- unique(cdat$CU_Name)
for(i in 1:14) {
  dat <- cdat[cdat$CU_Name==CUs[i], ] 
  par(mar=c(1,4,1,0)+0.3, bty="l", las=1)
  plot(y=dat$SiteEsc/1000, x=dat$Year, type="b", xlab="year", pch=16, ylab="Escapement (thousands)", 
       main=CUs[i], xlim =c(yrs[1], yrs[2]))
     grid(ny=0)
  par(mar=c(2,4,1,0)+0.3, bty="l", las=1)
  plot(y=dat$RS, x=dat$Year, log="y", type="l", ylab="Recruits/Spawner", xlim =c(yrs[1], yrs[2]), 
      xlab="Year", ylim=c(min(cdat$RS, na.rm=TRUE), max(cdat$RS, na.rm=TRUE)), yaxt="n")
  axis(side=2, at=c(0.01, 0.1, 1, 10, 100), labels=c(0.01, 0.1, 1, 10, 100))
  points(y=dat$RS, x=dat$Year, pch=16, col= colors[ rank[cdat$CU_Name==CUs[i]] ])
    grid(ny=0)
    abline(h=1, lty=3)
    abline(h=median(dat$RS, na.rm=TRUE),lty=2, col="orange")
    text(x=2013, y= median(dat$RS, na.rm=TRUE), col="orange", label=paste0("Median\nR/S=", round(median(dat$RS, na.rm=TRUE), 1)), adj=c(1,-0.5))
}
dev.off()


# -----------#
# Compare sum of max escapement summed across streams to infilling data - LW 
# -----------#

head(Btable)

head(rdl)
CUs <- unique(rdl$CU_Name)

sum_escp <- rdl %>% group_by(NME, CU_Name) %>% summarise(max_escape = max(escapement, na.rm=TRUE), min_escape = min(escapement, na.rm=TRUE)) %>% 
  group_by(CU_Name) %>% summarise(sum_max_escape= sum(max_escape, na.rm=TRUE), sum_min_escape=sum(min_escape, na.rm=TRUE))

Btable_max_escp <- merge(Btable, sum_escp, by.x="CU", by.y="CU_Name", all.x=TRUE)
head(Btable_max_escp)

Btable_max_escp$RS <- Btable_max_escp$Recruit/ Btable_max_escp$Escape
Btable_max_escp$logRS <-log(Btable_max_escp$RS)

# Plot recruits vs. spawners using infilled data, add vertical lines at 
# sum(max observed escapement for each stream) and sum(min observed escapement for each stream)
ggplot(Btable_max_escp, aes(x=Escape, y=Recruit,  fill=logRS)) + 
  geom_point(aes(size=logRS), shape=21, alpha=0.5) + 
  scale_fill_viridis_c() +
  geom_vline( aes(xintercept=sum_max_escape), colour="dodgerblue") +
  geom_vline(aes(xintercept=sum_min_escape), colour="firebrick") +
  geom_text(aes(x =sum_min_escape, y=0, label="Min obs spawners", angle=90, hjust=0)) + 
  geom_text(aes(x =sum_max_escape, y=0, label="Max obs spawners", angle=90, hjust=0)) + 
  geom_abline(aes(intercept=0, slope=1), colour="gray", lty=2) +
  facet_wrap(~CU, scales="free")

ggplot(Btable_max_escp, aes(x=Year, y=RS)) + 
#ggplot(Btable_max_escp, aes(x=Year, y=RS, fill=logRS)) + 
  geom_path() +
  geom_point(shape=21) + 
  #geom_point(shape=21, alpha=0.5, aes( size=logRS)) + 
  scale_fill_viridis_c() +  geom_hline(aes(yintercept=1)) +
  coord_cartesian(ylim=c(0,5)) +
  #scale_y_log10() +
  facet_wrap(~CU ) 

med_RS <- Btable_max_escp %>% group_by(CU) %>% summarise(median_RS = median(RS, na.rm=TRUE))
ggplot(Btable_max_escp, aes(x=RS)) +
  geom_density(colour="dodgerblue") + 
  geom_vline(aes(xintercept=1), lty=2) +
  scale_x_log10() +
  geom_vline(data=med_RS, aes(xintercept=median_RS), colour="dodgerblue") +
  facet_wrap(~CU) 


# Look at difference between recruits and escapement (Fisheries harvest)
png("Figures/fig_escapement_returns_harvest.png", height=8, width=12, units="in", res=300)
ggplot(Btable, aes(x=Year, y=Return)) +
  geom_point() +
  geom_path() +
  geom_linerange(aes(x=Year, ymin=Escape, ymax=Return), colour="red", lwd=1.3) +
  geom_ribbon(aes(x=Year, ymin=0, ymax=Escape), fill="dodgerblue", colour="dodgerblue") +
  #geom_text(aes(x=Year, y=Return+5000, label=round((Return-Escape)/Return*100,0)), size=4) +
  facet_wrap(~CU, scales="free_y")
dev.off()

# How would an alternative to Sgen that is spawners required to get to 
# a benchmark like 50% of max observed spawners compare to Sgen?
# incorporates alpha. 
# read in alpha and beta estimates from LRP integrated model
ests <- read.csv("DataIn/ricker_est_to_compare_with_beta.csv")
head(ests)
sum_escp <- sum_escp[1:7, ]

CUs <- unique(ests$CU_name)


curve( ests$A[1] * x * exp( - ests$B[1] * x), xlim=c(0,100000))
abline(h=sum_escp$sum_max_escape[1]*0.2 )
abline(h=sum_escp$sum_max_escape[1]*0.2 )
abline(v= ests$Sgen[1], col="grey")


function()

layout(mat= matrix(1:9, byrow=TRUE))
for (i in seq_along(CUs)) {
  plot(x=0,y=0,xlim=c(0, sum_escp$sum_max_escape[sum_escp$CU_Name==i]*1.1), ylim=c(0,sum_escp$sum_max_escape[sum_escp$CU_Name==i]*6))
  abline(v=sum_escp$sum_max_escape[sum_escp$CU_Name==i], col="dodgerblue")
  abline(v=sum_escp$sum_min_escape[sum_escp$CU_Name==i], col="firebrick")
  abline(v= ests$Sgen[sum_escp$CU_Name==i], col="grey")
  curve( ests$A[sum_escp$CU_Name==i] * x * exp( - ests$B[sum_escp$CU_Name==i] * x))
  main(i) 
}

# ----------#
# Look at age composition - LW
# ----------#
head(ACdat)

ACdat %>% pivot_longer(cols=grep("Age", names(ACdat)), names_to="Age", values_to="proportion") %>%
  ggplot(., aes(y=proportion, x=Year, fill=Age)) +
  geom_col()







