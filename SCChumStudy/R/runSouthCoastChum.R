# Code by Luke Warkentin & Kendra Holt
# This file contains the code required to explore South Coast Chum escapement datasets and run retrospective analyses of 
#   data-based LRPs that use logistic regressions

library(rsample)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(TMB)
library(here)

options(scipen=1000000)

chumDir <- here() # needed for some plot functions
codeDir <- file.path(dirname(here()), "Code") # get code directory in project parent folder
r_to_source <- list.files(codeDir, pattern="\\.R|\\.r") # get r files to source from parent code directory
# function to source all files required
sourceAll <- function(dir, files){
  fpaths <- paste0( dir, "/", files)
  sapply(fpaths, source)
}
sourceAll(dir=codeDir, files= r_to_source) # source all r files from parent directory

# Load TMB models
# make vector of TMB cpp file names without .cpp extension
TMB_for_chum <- c("SR_IndivRicker_NoSurv",
                  "SR_IndivRicker_NoSurv_LowAggPrior",
                  "LRP_Logistic_Only",
                  "LRP_Logistic_Only_LowAggPrior")
# compile and load the TMB .cpp files
compile_load_TMB(dir=codeDir, files=TMB_for_chum)

# ====================================================================#
# Read in data and format for using in retrospective analysis------------
# ====================================================================#
# Run to re-do infilling and brood table
#source("R/make_brood_table.r") 

# Read in chum wild escapement data (infilled) by CU
ChumEscpDat <- read.csv("DataOut/wild_spawners_CU_infilled_by_site_CU.csv")
ChumEscpDat$MU <- "SC Chum" # Add MU column - FLAG this may not be necessary - can't find MU in any of the function scripts. Check with Kendra
#ChumEscpDat$CU_ID <- as.integer(substr(ChumEscpDat$CU_Name, 1,1)) # pull out CU ID from raw CU column name. Need CU_ID column for percentile benchmarks
ChumEscpDat$CU <- as.integer(substr(ChumEscpDat$CU_Name, 1,1))  # need CU column for some functions
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

# Call functions to plot data availability (note that this is for infilled data, plots hide missing observations):
# plot_CU_DataObs_Over_Time(ChumEscpDat, chumDir, plotName="SC_Chum_DataByCU")
# plot_Num_CUs_Over_Time(ChumEscpDat, chumDir, plotName="SC_Chum_N_CUs")
# 
# # Note: these next 2 two escpt plots need to have formatting fixed
# plot_CU_Escp_Over_Time(ChumEscpDat, chumDir, plotName="SCChum Esc", samePlot = T)
# plot_CU_Escp_Over_Time(ChumEscpDat, chumDir, plotName="SCChum Esc Separate", samePlot = F)

# Setup for retrospective analysis
# Data prep

# FLAG: Should probably limit stock-recruit data to year > 1959/1960 to allow for full brood year returns up to age 6. 
# This may be done automatically, see retroFunctions.r line 20 
# This is done automatically using the BroodYrLag variable 
# First make full time series to use with percentile benchmarks (don't need recruit parameters so can use 1953-1957)
ChumEscpDat_full <- ChumEscpDat
ChumSRDat_full <- ChumSRDat

# remove years without full recruitment  - FLAG: still needing to run this, otherwise have NA recruitment going into model optimization
ChumEscpDat <- ChumEscpDat[ChumEscpDat$yr >= 1958 ,] # remove years without full recruitment
ChumSRDat <- ChumSRDat[ChumSRDat$BroodYear >= 1958 ,] # remove years without full recruitment

# make data frames that do not include any observations from CUs that have any amount of CU-level infilling (Upper Knight and Bute Inlet)
CUs_not_use <- unique(ChumEscpDat$CU_Name[which(is.na(ChumEscpDat$Escape))]) # get CUs that have NA values for escapement in some years (means that they were infilled at CU level)
ChumEscpDat_no_CU_infill <- ChumEscpDat[!(ChumEscpDat$CU_Name %in% CUs_not_use), ]
ChumSRDat_no_CU_infill <- ChumSRDat[!(ChumSRDat$CU_Name %in% CUs_not_use), ]

# make data frames that do not include any observations from years that have any amount of CU-level infilling
years_not_use <- unique(ChumEscpDat$yr[which(is.na(ChumEscpDat$Escape))])
# need to adjust this so that it also doesn't have brood year
# We know which years use CU-level infilling for escapement. 
# Now, which years use those CU-level infilled escapement data 
# to reconstruct the brood table? 
# Any years that use these escapement values in the reconstruction of 
# spawners or recruits. 
# Fish return at 3,4,5, or 6 years old. 
# Each brood year x with CU-level infilling cannot be used for spawners or recruits.
# Cannot use years x-(3:6) which are the brood years with recruits returning as 3-6 yr olds in year x
# cannot use years x+(3:6) which are the years with 3-6 year old recruits from brood year x
linked_yrs <- function(no_count_yrs) {
  y <- as.vector(sapply(no_count_yrs, function(x){
      c((x-6):(x-3),x,(x+3):(x+6))
  }))
  unique(y)
}

linked_yrs(2000)

years_not_use2 <- linked_yrs(years_not_use)

ChumEscpDat_no_CU_infill_yrs <- ChumEscpDat[!(ChumEscpDat$yr %in% years_not_use2), ]
ChumSRDat_no_CU_infill_yrs <- ChumSRDat[!(ChumSRDat$BroodYear %in% years_not_use2), ]
length(unique(ChumEscpDat_no_CU_infill_yrs$yr))
# Note that there are only 15 years of data after removing
# years with CU-level infilling. May make it not feasible to
# check effect. 

# Get data frame that only drops the least monitored CU (Upper Knight) and the years 
# (and associated spawner/recruit years as above) with CU-level infilling in Bute Inlet
# removed
# get years that Bute wasn't monitored
# years_not_use_bute <- ChumEscpDat$yr[intersect(which(is.na(ChumEscpDat$Escape)), which(ChumEscpDat$CU_Name=="Bute Inlet"))]
# expand these years to include recruits/spawners that are linked to these years
# years_not_use_bute2 <- linked_yrs(years_not_use_bute)
# remove Upper Knight, and any years linked to uncounted years in Bute Inlet
# ChumEscpDat_no_knight <- ChumEscpDat[ !(ChumEscpDat$CU_Name=="Upper Knight" | ChumEscpDat$yr %in% years_not_use_bute2), ]
# ChumSRDat_no_knight <- ChumSRDat[  !(ChumSRDat$CU_Name=="Upper Knight" | ChumSRDat$BroodYear %in% years_not_use_bute2),  ]

# Specify p value for logistic regression
ps <- c(seq(0.6, 0.95,.05), 0.99) 

# ========================================================================#
# Run annual restrospective analysis using stock recruit parameter based Sgen----------
# ========================================================================#

TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1, 
                      extra_eval_iter=TRUE)

# Loop over p values and run annual retrospective analyses for each level of p

# LW Notes
# Choosing values for runAnnualRetro function:
# - genYrs: average age of return.  runFraserCoho.r uses 3. I use 4 for chum (age 4 returns made up 69% of returns based on 2013 data). Update for 2018 age composition data.
# - BroodYrLag: How many years will it take between a fish returning and when you can construct the full recruitment.
#     The BroodYrLag parameter used when prepping for the LRP calculations represents the number of years 
#     it will take to reconstruct recruitment from a single brood year  (i.e., how many years of recruitment 
#     do you need to observe to estimate recruitment from a brood year).  Another way of defining it is the 
#     number of age classes that recruitment (maturation) occurs over. 
#     For the coho example, fish mostly return at age 3 or 4, so BroodYrLag is 2.
#     For the chum example, it looks like fish return at ages 3,4,5,6, so BroodYrLag is 4.
#     This parameter gets used in the runAnnualRetro() function to remove the first few years of data for 
#     which recruitment is NA because the BY has not yet been fully observed.

# Run retrospective analysis
for(pp in 1:length(ps)){
  # Run with Binomial LRP model with individual model Ricker
  # runAnnualRetro(EscpDat=ChumEscpDat, SRDat=ChumSRDat, startYr=1975, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_NoSurv", LRPmodel="BinLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)

  # Run with CUs with CU-level infilling removed
  runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill, SRDat=ChumSRDat_no_CU_infill, startYr=1975, endYr=2012, BroodYrLag=4, genYrs=4, p = ps[pp],
                BMmodel = "SR_IndivRicker_NoSurv", LRPmodel="BinLogistic", integratedModel=T,
                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_noCUinfill_",ps[pp]*100, sep=""),
                bootstrapMode = F, plotLRP=T, runLogisticDiag = TRUE)
  # for some reason this gives errors if startYr is set to less than 1975. Perhaps there aren't enough stock-recruit 
  # observations for the model to fit? The residuals are large even when there are more points.
  
  # # Run with years with CU-level infilling removed
  # runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill_yrs, SRDat=ChumSRDat_no_CU_infill_yrs, startYr=1967, endYr=1972, BroodYrLag=4, genYrs=4, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_NoSurv", LRPmodel="BinLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_noCUinfill_yrs_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T, runLogisticDiag=TRUE)
  # 
  # # Run with no Upper Knight CU, no CU-level infilled yrs
  # runAnnualRetro(EscpDat=ChumEscpDat_no_knight, SRDat=ChumSRDat_no_knight, startYr=1967, endYr=1998, BroodYrLag=4, genYrs=4, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_NoSurv", LRPmodel="BinLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_no_knight",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T, runLogisticDiag=TRUE)
  # 
  # Run with Bernoulli LRP model with individual model Ricker
  runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill, SRDat=ChumSRDat_no_CU_infill, startYr=1975, endYr=2012, BroodYrLag=4, genYrs=4, p = ps[pp],
                 BMmodel = "SR_IndivRicker_NoSurv", LRPmodel="BernLogistic", integratedModel=T,
                 useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bern.IndivRicker_NoSurv_noCUinfill_",ps[pp]*100, sep=""),
                 bootstrapMode = F, plotLRP=T, runLogisticDiag=TRUE)
  # for some reason this gives errors if startYr is set to less than 1975. Perhaps there aren't enough stock-recruit 
  # observations for the model to fit? The residuals are large even when there are more points.
}

# -----------------------------------------------------#
# Run Sgen LRP with low aggregate penalty 
# -----------------------------------------------------#
# Get values for low aggregate likelihood penalty model 
# Idea is to parameterize (mu and sigma) the aggregate abundance 
# associated with a very low proportion (essentially 0; p=0.01) of
# CUs being above their benchmark. The idea is to have 95% of this 
# estimate between the CU Sgen and the sum of all their Sgen 
# benchmarks or upper estimates of Sgens (e.g., all CUs are just 
# below their lower benchmarks). The mean of this distribution is 
# halfway between these lower (2.5% quantile) and upper (97.5% quantile) values. 
# 
# # get Sgen estimates from running integrated model without penalty
# ests <- read.csv("DataOut/AnnualRetrospective/Bin.IndivRicker_NoSurv_noCUinfill_90/annualRetro_SRparsByCU.csv", stringsAsFactors = FALSE)
# 
# ests1 <- ests[ests$retroYr ==max(ests$retroYr),] # get just last retro year for estimates
# # make lower limit the lowest CU Sgen (this gives essentially the same results as 
# # using the average abundance of the smallest CU)
# low_lim <- min(ests1$est_Sgen) 
# # make upper limit the sum of CU benchmarks
# #hi_lim <- sum(ests1$est_Sgen) # sum Sgen estimates to get upper limit for penalty mu
# # sum upper CI of Sgen estimates to get upper limit for penalty mu
# hi_lim <- sum(na.omit(ests1$up_Sgen)) 
# 
# # make mu of penalty value mean of these two values, divide by scale
# B_penalty_mu <- mean(c(low_lim, hi_lim))
# # get SD that gives 95% density between lower and upper limits (using getSD helper function)
# dum<-optim(par=200, fn = getSD, method="Brent",lower=1, upper=50000000, low_lim=low_lim, hi_lim=hi_lim)
# B_penalty_sigma<-dum$par
# 
# # get SD for prior penalty so that 95% density is between lower and upper limits
# # plot to check
# plot( x = seq( 0, max(ChumSRDat$Spawners), 100), y=dnorm( seq(0, max(ChumSRDat$Spawners), 100), mean=mean(c(low_lim, hi_lim)), sd=B_penalty_sigma), 
#       xlim=c(0, max(ChumSRDat$Spawners)), type="l", ylab="density", xlab="aggregate adundance")
# abline(v=c(low_lim, hi_lim, mean(c(low_lim, hi_lim))), col="dodgerblue", lty=c(2,2,1)) # plot upper and lower values and mean
# 
# # Add penalty values to TMB_Inputs_IM
# TMB_Inputs_IM_LowAggPrior <- list(Scale = 1000, logA_Start = 1,
#                       Tau_dist = 0.1,
#                       gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
#                       B_penalty_mu = B_penalty_mu, B_penalty_sigma = B_penalty_sigma, extra_eval_iter=TRUE)
# 
# # OBSOLETE ->
# #B_penalty_sigmas <- c(B_penalty_sigma, B_penalty_sigma* 1.5, B_penalty_sigma * 2) # test with 1.5 and 2 times SD
# # Effect of increasing SD:
# #   Without CUs with CU infilling:
# #       SD times 1.5 and 2 make very minor differences in LRP (very small increases). No flipping of logistic curve
# #   With using CUs with CU infilling:
# #    Using 1.5 times SD makes logistic curve flip in one year, and 2 times SD makes it flip several years (with p=0.9)
# # <- OBSOLETE
#
# # Run retrospective analysis with likelihood penalty on aggregate abundance at low proportion CUs above benchmark 
# #       (to bring logistic curve intercept down)
# 
# for(pp in 1:length(ps)){
#   #for (i in 1:length(B_penalty_sigmas)){ # for testing sensitivity to B_penalty_sigma
#   # Run with Binomial LRP model with individual model Ricker without CU-level infilling
#   runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill, SRDat=ChumSRDat_no_CU_infill, startYr=1967, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
#                 BMmodel = "SR_IndivRicker_NoSurv_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
#                 useGenMean=F, TMB_Inputs=TMB_Inputs_IM_LowAggPrior, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_LowAggPrior_noCUinfill_", ps[pp]*100, sep=""),
#                 bootstrapMode = F, plotLRP=T , runLogisticDiag=F)
#   # Run without Upper Knight CU
#   # runAnnualRetro(EscpDat=ChumEscpDat_no_knight, SRDat=ChumSRDat_no_knight, startYr=1967, endYr=1998, BroodYrLag=4, genYrs=4, p = ps[pp],
#   #                BMmodel = "SR_IndivRicker_NoSurv_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
#   #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM_LowAggPrior, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_LowAggPrior_no_knight_", ps[pp]*100, sep=""),
#   #                bootstrapMode = F, plotLRP=T , runLogisticDiag=F)
#   #}
# }

# ========================================================================#
# Run aggregate LRP with Percentile benchmarks-----------
# ========================================================================#

# Currently there is a discrepancy between CU_ID variable in percentile/subpop and Sgen 
# based retrospective. Need to sort this out. 

# First, make data frame of CUs and which percentile to use for their lower benchmark,
# based on their productivity (Ricker alpha) and mean exploitation rate (ER), 
# according to Table 6 in Holt et al. 2018
# Here, I do this manually, taking alpha values from data outputs from retro year 2010 and from mean of exploitation rates
#
# Southern Coastal Streams: alpha = 1.3, ER = 0.04, benchmark = NA (further evaluation required) 
#
# Northeast Vancouver Island: alpha = 1.7, ER = 0.05, benchmark = 50%
#
# Upper Knight: alpha = 2.6 (note this is higher than Holt et al. 2018, likely due to recent CU-level infilling escapement;
#       choosing alpha = 2.22 value from Holt et al. 2018 table 2 to be conservative. It doesn't really matter for the logistic 
#       regression because this CU is not used in the logistic regression because of CU-level infilling), ER = 0.07, benchmark = 50%
#
# Loughborough: alpha = 2.4, ER = 0.18, benchmark = 50%
#
# Bute Inlet: alpha = 2.7 (note this is higher than Holt et al. 2018, likely due to recent CU-level infilling escapement;
#       choosing alpha = 2.46 value from Holt et al. 2018 table 2 to be conservative. It doesn't really matter for the logistic 
#       regression because this CU is not used in the logistic regression because of CU-level infilling), ER = 0.22, 
#       benchmark = NA (further evaluation required) 
#
# Georgia Strait: alpha = 3.3, ER = 0.38, benchmark = 25%
#
# Howe Sound-Burrard Inlet: alpha = 2.6, ER = 0.32, benchmark = 25%

# Make data frame with CU names and which percentile to use for 
which_perc_benchmark <- data.frame("CU" = c("Southern Coastal Streams", "North East Vancouver Island", "Upper Knight", 
  "Loughborough", "Bute Inlet", "Georgia Strait", "Howe Sound-Burrard Inlet"),
  "percentile" = c(NA, 0.5, 0.5, 0.5, NA, 0.25, 0.25))
# also make one with 0.75 for the two CUs where percentiles are not recommended, to make figure with all 7 CUs
which_perc_benchmark_hack <- data.frame("CU" = c("Southern Coastal Streams", "North East Vancouver Island", "Upper Knight", 
                                            "Loughborough", "Bute Inlet", "Georgia Strait", "Howe Sound-Burrard Inlet"),
                                   "percentile" = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.25, 0.25))

# Remove any CUs that have NA percentile benchmark (percentile benchmarks not recommended in Holt et al. 2018 Table 6)
ChumEscpDat_perc <- ChumEscpDat_no_CU_infill[ !(ChumEscpDat_no_CU_infill$CU_Name %in% which_perc_benchmark$CU[is.na(which_perc_benchmark$percentile)] ), ]
ChumSRDat_perc <- ChumSRDat_no_CU_infill[ !(ChumSRDat_no_CU_infill$CU_Name %in% which_perc_benchmark$CU[is.na(which_perc_benchmark$percentile)] ), ]

TMB_Inputs_Percentile <- list(Scale = 1000, logA_Start = 1,
                              Tau_dist = 0.1,
                              gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                              B_penalty_mu = NA, B_penalty_sigma = NA, perc_benchmark= which_perc_benchmark)
TMB_Inputs_Percentile_hack <- list(Scale = 1000, logA_Start = 1,
                              Tau_dist = 0.1,
                              gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                              B_penalty_mu = NA, B_penalty_sigma = NA, perc_benchmark= which_perc_benchmark_hack)
# Run retrospective analysis using percentile benchmarks
for(pp in 1:length(ps)){
  # Run with Binomial LRP with CUs with CU-level infilling removed # changed endYr to 2018 as recruits not needed
  # runAnnualRetro(EscpDat=ChumEscpDat_perc, SRDat=ChumSRDat_perc, startYr=1967, endYr=2018, BroodYrLag=4, genYrs=4, p = ps[pp],
  #                BMmodel = "Percentile", LRPmodel="BinLogistic", LRPfile="LRP_Logistic_Only",integratedModel=F,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_Percentile, outDir=chumDir, RunName = paste("Bin.Percentile_noCUinfill_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T,  runLogisticDiag=T)
  # # Run with CU-level infilling (to make plots with percentile benchmarks over time for all 7 CUs)
  # # This is just for making the plot with status for all 7 CUs
  runAnnualRetro(EscpDat=ChumEscpDat, SRDat=ChumSRDat, startYr=1967, endYr=2018, BroodYrLag=4, genYrs=4, p = ps[pp],
                 BMmodel = "Percentile", LRPmodel="BinLogistic", LRPfile="LRP_Logistic_Only",integratedModel=F,
                 useGenMean=F, TMB_Inputs=TMB_Inputs_Percentile_hack, outDir=chumDir, RunName = paste("Bin.Percentile_",ps[pp]*100, sep=""),
                 bootstrapMode = F, plotLRP=T, runLogisticDiag=T)
  # Run without Upper Knight CU
  # runAnnualRetro(EscpDat=ChumEscpDat_no_knight, SRDat=ChumSRDat_no_knight, startYr=1967, endYr=1998, BroodYrLag=4, genYrs=4, p = ps[pp],
  #                BMmodel = "Percentile", LRPmodel="BinLogistic", LRPfile="LRP_Logistic_Only",integratedModel=F,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_Percentile, outDir=chumDir, RunName = paste("Bin.Percentile_no_knight_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
  # 
  # Run with Bernouli regression with CUs with CU-level infilling removed # changed endYr to 2018 as recruits not needed
  # runAnnualRetro(EscpDat=ChumEscpDat_perc, SRDat=ChumSRDat_perc, startYr=1967, endYr=2018, BroodYrLag=4, genYrs=4, p = ps[pp],
  #                BMmodel = "Percentile", LRPmodel="BernLogistic", LRPfile="LRP_Logistic_Only",integratedModel=F,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_Percentile, outDir=chumDir, RunName = paste("Bern.Percentile_noCUinfill_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T, runLogisticDiag=T)
  
}

# ---------------------------------------------#
# Run percentile LRP with low aggregate penalty values
# ---------------------------------------------#
# Get low aggregate penalty values
# pdat <- read.csv("DataOut/AnnualRetrospective/Bin.Percentile_noCUinfill_60/annualRetro_perc_benchmarks.csv", stringsAsFactors = FALSE)
# pdat1 <- pdat[pdat$retro_year ==max(pdat$retro_year),]
# pdat2 <- pdat1 %>% group_by(CU_Name, CU, benchmark_perc_25) %>% summarise(n())
# # make lower limit=0, upper limit sum of 25% benchmarks for retro year 2010
# low_lim_perc <- 0
# hi_lim_perc <- sum(pdat2$benchmark_perc_25)
# 
# B_penalty_perc_mu <- mean(c(low_lim_perc, hi_lim_perc)) 
# dum_perc<-optim(par=200, fn = getSD, method="Brent",lower=1, upper=1000000, low_lim=low_lim_perc, hi_lim=hi_lim_perc)
# B_penalty_perc_sigma<-dum_perc$par
# 
# TMB_Inputs_Percentile <- list(Scale = 1000, logA_Start = 1,
#                               Tau_dist = 0.1,
#                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
#                               B_penalty_mu = B_penalty_perc_mu, B_penalty_sigma = B_penalty_perc_sigma)
# 
# # Run retrospective analysis using percentile benchmarks, with likelihood penalty
# for(pp in 1:length(ps)){
#   # with prior penalty on low aggregate abundance
#   runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill, SRDat=ChumSRDat_no_CU_infill, startYr=1967, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
#                  BMmodel = "Percentile", LRPmodel="BinLogistic", LRPfile="LRP_Logistic_Only_LowAggPrior", integratedModel=F,
#                  useGenMean=F, TMB_Inputs=TMB_Inputs_Percentile, outDir=chumDir, RunName = paste("Bin.Percentile_LowAggPrior_noCUinfill_",ps[pp]*100, sep=""),
#                  bootstrapMode = F, plotLRP=T, runLogisticDiag=T)
# }

# ========================================================================#
# Make figures from outputs ---------- 
# ========================================================================#

# ---------------------------------------------#
# Plot percentile benchmark over time----------
# ---------------------------------------------#
pdat <- read.csv("DataOut/AnnualRetrospective/Bin.Percentile_60/annualRetro_perc_benchmarks.csv", stringsAsFactors = FALSE)
#pdat <- read.csv("DataOut/AnnualRetrospective/Bin.Percentile_noCUinfill_60/annualRetro_perc_benchmarks.csv", stringsAsFactors = FALSE)

#cudf <- data.frame("CU_Name_num" = unique(pdat$CU_Name), "b" = unique(pdat$CU_Name))
#pdat2 <- merge(pdat, cudf, by.x="CU_Name", by.y="b") # fix names so consistent with other plots

# add in actual percentile to use (Southern Coastal streams + Bute are NA, used 0.5 as a filler to run retrospective) )
pdat2 <- merge(pdat, which_perc_benchmark, by.x="CU_Name", by.y="CU", all.x=TRUE)
pdat2$perc_appr <- ifelse(is.na(pdat2$percentile)==FALSE, 1,0) # add variable that is 1 if percentile benchmarks are appropriate, 0 if not
# add column to chum escapement data that is true if CU infilling was done
ChumEscpDat_full$CU_infill <- ifelse(is.na(ChumEscpDat_full$Escape), TRUE, FALSE)
pdat3 <- merge(pdat2, ChumEscpDat_full[,names(ChumEscpDat_full) %in% c("CU_Name","yr", "CU_infill")], by=c("CU_Name","yr"), all.x=TRUE)

# Get data frame that has one row per year, and the corresponding retro year benchmarks. 
# Also has the rows from before retro years (years used up to first retro year)
pbm <- pdat3[union(grep(min(pdat3$retro_year), pdat3$retro_year), which(pdat3$retro_year==pdat3$yr)), ]

png("Figures/fig_perc_benchmarks_annual_retro.png", width=6.5, height=7, res=300, units="in", pointsize=30)
ggplot(pbm, aes(y=benchmark_perc_25, x=retro_year, shape=as.factor(perc_appr))) +
  geom_line(colour="gray", linetype=2) +
  geom_line(aes(y=Escp, x=yr), size=0.2) +
  geom_line(aes(y=benchmark_perc_50, x=retro_year), linetype=2) +
  geom_point(aes(y=Escp, x=yr, fill=as.factor(AboveBenchmark),  alpha=CU_infill)) +
  scale_fill_manual( values=c("red", "darkgreen"), guide=NULL) +
  scale_shape_manual(values=c(1, 21), guide="none") +
  scale_alpha_manual(values=c(1,0.3), guide="none") +
  facet_wrap(~paste0(CU_Name, " (", ifelse(is.na(percentile),"NA)", paste0(use_perc*100,"%)"))), scales="free_y", ncol=2) +
  ylab("Escapement with 25% and 50% benchmarks") + xlab("Year") +
  coord_cartesian(expand=FALSE, clip="off") +
  geom_hline(aes(yintercept=0)) +
  theme_classic() +
  theme(strip.background = element_blank())
dev.off()

# --------------------------------------------------------------#
# Plot alpha, beta, SMSY and Sgen over time---------
# --------------------------------------------------------------#
mdat <- read.csv("DataOut/AnnualRetrospective/Bin.IndivRicker_NoSurv_noCUinfill_60/annualRetro_SRparsByCU.csv", stringsAsFactors = FALSE)
mdat1 <- mdat %>% pivot_longer(cols=est_B:up_Sgen, names_to="param", values_to="est")

# Plot alpha, beta, SMSY and Sgen on one plot
png("Figures/fig_a_b_SMSY_Sgen_retro.png", width=6, height=8, res=300, units="in", pointsize=12)
layout(mat=matrix(1:18, byrow=FALSE, ncol=3)) # 15 panels
CUs <- unique(mdat$CU_Name)
brk_a <- seq(0,max(mdat$est_A)*1.05, 1)
#brk_b <- seq(0,max(mdat$est_B)*1.05, 0.00001)
par(las=1)
for(i in 1:2) {
  d <- mdat[ mdat$CU_Name == CUs[i], ]# get CU data
  par(mar=c(0,5,2,0)+0.2)
  plot(x=d$retroYr, y = d$est_A, main=CUs[i], cex.main=1, yaxs="i", ylab="", type="l", xaxt="n",  bty="n", ylim=c(min(brk_a), max(brk_a)) )
  mtext(text=expression(alpha),side=2, line=4.4 )
  #axis(side=2, at=brk_a, labels=brk_a)
  par(mar=c(1,5,1,0)+0.2)
  plot(x=d$retroYr, y = d$est_B, ylab="", yaxs="i", type="l", xaxt="n",bty="n", ylim=c(0, max(d$est_B)) )
  mtext(text=expression(beta),side=2, line=4.4 )
  #axis(side=2, at=brk_b, labels=brk_b)
  par(mar=c(2,5,0,0)+0.2)
  plot(x=d$retroYr, y = d$est_Smsy, col="dodgerblue", xlab="", ylab=expression(S[gen]*" and SMSY"), type="l", ylim=c(min(c(d$est_Smsy, d$low_Sgen)), max(c(d$est_Smsy, d$up_Sgen))), bty="l" )
  lines(x=d$retroYr, y = d$est_Sgen)
  polygon(x = c(d$retroYr, rev(d$retroYr)), y=c(d$up_Sgen, rev(d$low_Sgen)), col=adjustcolor(col="gray",alpha=0.3), border=NA)
}
for(i in 3:length(CUs)) {
  d <- mdat[ mdat$CU_Name == CUs[i], ]# get CU data
  par(mar=c(0,4,2,0)+0.2)
  plot(x=d$retroYr, y = d$est_A, main=CUs[i],cex.main=1, ylab="", yaxs="i",type="l", xaxt="n", bty="n" , ylim=c(min(brk_a), max(brk_a)))
  par(mar=c(1,4,1,0)+0.2)
  plot(x=d$retroYr, y = d$est_B, ylab="", type="l", yaxs="i",xaxt="n", bty="n" , ylim=c(0, max(d$est_B)) )
  par(mar=c(2,4,0,0)+0.2)
  plot(x=d$retroYr, y = d$est_Smsy, col="dodgerblue", ylab="", type="l", xlab="", ylim=c(min(c(d$est_Smsy, d$low_Sgen)), max(c(d$est_Smsy, d$up_Sgen))), bty="l" )
  lines(x=d$retroYr, y = d$est_Sgen)
  polygon(x = c(d$retroYr, rev(d$retroYr)), y=c(d$up_Sgen, rev(d$low_Sgen)), col=adjustcolor(col="gray",alpha=0.3), border=NA)
}
dev.off()
# 
# # Plot Sgen estimates over time, with SMSY
# png("Figures/fig_Sgen_annual_retro.png", width=12, height=3, res=300, units="in")
# ggplot(mdat, aes(y=est_Sgen, x=retroYr)) +
#   geom_line() +
#   geom_ribbon(aes(ymin=low_Sgen, ymax=up_Sgen, x=retroYr),colour=NA, alpha=0.2) +
#   geom_line(aes(y=est_Smsy, x=retroYr), colour='dodgerblue') +
#   facet_wrap(~CU_Name, scales="free_y", nrow=1) +
#   ylab("Sgen and SMSY") + xlab("Year") +
#   geom_hline(aes(yintercept=0), linetype=2) +
#   theme_bw()
# dev.off()
# 
# # Plot Alpha and beta on same plot
# png("Figures/fig_a_b_annual_retro.png", width=10, height=4, res=300, units="in")
# mdat1 %>% filter(param %in% c("est_A", "est_B")) %>%
#   ggplot(., aes(y=est, x=retroYr)) +
#   geom_line() +
#   geom_hline(aes(yintercept=0)) +
#   facet_wrap(param~CU_Name, scales="free_y", ncol=5)+
#   theme_bw()
# dev.off()

# --------------------------------------------------------------#
# Plot LRP status by year, compare different methods -----------
# --------------------------------------------------------------#

# Roll up escpaments, and get Gen Mean of that
AggEscp <- ChumEscpDat_full %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Escp, 3, gm_mean, fill = NA, align="right"))


# Plot annual status with bars to show years in which LRP was breached

# Read in data from multi-dimensional/ decision tree and simple percentile approach for comparison

md <- read.csv("DataIn/LRP_compare_methods.csv", header=TRUE)

md$AboveLRP <- ifelse(md$lrp_status=="above", TRUE, FALSE)

# Make new names for scenarios
snames <- c("(a) CUs with data in all years & perc. benchmarks",
            "(b) All years with data & perc. benchmarks",
            "(a) CUs with data in all years & perc. benchmarks",
            "(c) CUs with data in all years",
            "(d) All years with data",
            "(b) All years with data & perc. benchmarks")
skey <- data.frame(id = sort(unique(md$scenario)), name=snames) 
md$data_name <- skey$name[match(md$scenario, skey$id)]
md$data_name <- paste0(md$scenario_name, " - ", md$data_name)
plotStatusBarsChum_byYear(Status_DF = md, AggEscp=AggEscp, fName="fig_compare_LRP_methods")

# --------------------------------------------------------------#
# Plot ricker, SMSY, Sgen estimates from integrated model
# --------------------------------------------------------------#
ests <- read.csv("DataOut/AnnualRetrospective/Bin.IndivRicker_NoSurv_noCUinfill_60/annualRetro_SRparsByCU.csv", stringsAsFactors = FALSE)
ests1 <- ests[ests$retroYr ==max(ests$retroYr),] # get just one retro year for estimates
t <- merge(ests1, ChumSRDat[, names(ChumSRDat) %in% c("BroodYear", "Spawners", "Recruits", "CU_Name")], by=c("CU_Name"))

CUs <- unique(t$CU_Name)

# --------------------------------------------------------------#
# Plot recruits~spawner with ricker and Sgen and SMSY
# --------------------------------------------------------------#
png("Figures/fig_ricker_integrated_model.png", width=10, height=6, res=300, units="in")
layout(mat=matrix(1:8, byrow = TRUE, ncol=4))
par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
for(i in 1:length(unique(ests$CU_ID))) {
  dat <- t[t$CU_Name==CUs[i],]
  plot(dat$Recruits ~ dat$Spawners, type="p", main=CUs[i], xlab="Spawners", ylab="Recruits")
  curve( (dat$est_A[1] * x * exp(- dat$est_B[1] * x)), add=TRUE, lty=2)
  abline(v=dat$est_Sgen, lty=2, col="orange")
  abline(v=dat$est_Smsy, lty=2, col="dodgerblue")
}
  plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')
  legend(x=1,y=5, lty=2, col=c("orange", "dodgerblue"), legend=c("Sgen", "SMSY"))
dev.off()

# --------------------------------------------------------------#
# Plot log(recruits/spawner)~spawner with linear ricker and Sgen and SMSY
# --------------------------------------------------------------#
png("Figures/fig_linear_ricker_integrated_model.png", width=10, height=6, res=300, units="in")
layout(mat=matrix(1:8, byrow = TRUE, ncol=4))
par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
for(i in 1:length(unique(ests$CU_ID))) {
  dat <- t[t$CU_Name==CUs[i],]
  plot(log(dat$Recruits/dat$Spawners) ~ dat$Spawners, type="p", main=CUs[i], xlab="Spawners", ylab="log(Recruits/Spawner)")
  curve( log(dat$est_A[1]) - dat$est_B[1] * x, add=TRUE, lty=2)
  abline(v=dat$est_Sgen, lty=2, col="orange")
  abline(v=dat$est_Smsy, lty=2, col="dodgerblue")
}
plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')
legend(x=1,y=5, lty=2, col=c("black", "orange", "dodgerblue"), legend=c("linear Ricker", "Sgen", "SMSY"))
dev.off()

# --------------------------------------------------------------#
# Residual plots
# --------------------------------------------------------------#

# get residuals
t$pred_rec <- t$Spawners * t$est_A * exp(-t$est_B * t$Spawners) 
t$resid_rec <- t$Recruits - t$pred_rec

t$pred_RS <- t$est_A * exp(-t$est_B * t$Spawners) 
t$resid_RS <- t$Recruits/t$Spawners - t$pred_RS

t$pred_logRS <- log(t$est_A) - t$est_B * t$Spawners
t$resid_logRS <- log(t$Recruits/t$Spawners) - t$pred_logRS

# Plot residuals recruits ~ spawners
layout(mat=matrix(1:8, byrow = TRUE, ncol=4))
par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
for(i in 1:length(unique(ests$CU_ID))) {
  dat <- t[t$CU_Name==CUs[i],]
  plot(dat$resid_rec ~ dat$Spawners, type="p", main=CUs[i], xlab="Spawners", ylab="Recruits, obs-preds")
  abline(v=dat$est_Sgen, lty=2, col="orange")
  abline(v=dat$est_Smsy, lty=2, col="dodgerblue")
  abline(h=0, lty=2)
}
plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')
legend(x=1,y=5, lty=2, col=c("orange", "dodgerblue"), legend=c("Sgen", "SMSY"))

# Plot residuals recruits/spawner ~ spawners
layout(mat=matrix(1:8, byrow = TRUE, ncol=4))
par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
for(i in 1:length(unique(ests$CU_ID))) {
  dat <- t[t$CU_Name==CUs[i],]
  plot( (dat$Recruits/dat$Spawners - dat$pred_rec/dat$Spawners) ~ dat$Spawners, type="p", main=CUs[i], xlab="Spawners", ylab="R/S, obs-preds")
  abline(v=dat$est_Sgen, lty=2, col="orange")
  abline(v=dat$est_Smsy, lty=2, col="dodgerblue")
  abline(h=0, lty=2)
}
plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')
legend(x=1,y=5, lty=2, col=c("orange", "dodgerblue"), legend=c("Sgen", "SMSY"))

# # Plot residuals recruits/spawner ~ recruits/spawner
# layout(mat=matrix(1:8, byrow = TRUE, ncol=4))
# par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
# for(i in 1:length(unique(ests$CU_ID))) {
#   dat <- t[t$CU_Name==CUs[i],]
#   plot( y= dat$resid_RS, x= dat$Recruits/dat$Spawners, type="p", main=CUs[i], xlab="observed R/S", ylab="resid R/S")
#   abline(h=1, lty=2)
# }
# plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')

# Plot resid log(R/S) ~ spawners
layout(mat=matrix(1:8, byrow = TRUE, ncol=4))
par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
for(i in 1:length(unique(ests$CU_ID))) {
  dat <- t[t$CU_Name==CUs[i],]
  plot(dat$resid_logRS ~ dat$Spawners, type="p", main=CUs[i], xlab="Spawners", ylab="resid log(R/S)")
  abline(h=0, lty=2)
  abline(v=dat$est_Sgen, lty=2, col="orange")
  abline(v=dat$est_Smsy, lty=2, col="dodgerblue")
}
plot(x=1:10, y=1:10, type='n', ann = FALSE, xaxt = 'n', yaxt = 'n', bty = 'n')
legend(x=1,y=5, lty=2, col=c("orange", "dodgerblue"), legend=c("Sgen", "SMSY"))

# Plot residuals ~ observed recruits
layout(mat=matrix(1:8, byrow = TRUE, ncol=4))
par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
for(i in 1:length(unique(ests$CU_ID))) {
  dat <- t[t$CU_Name==CUs[i],]
  plot(dat$resid_rec ~ dat$Recruits, type="p", main=CUs[i], xlab="Recruits", ylab="Recruits, obs-pred")
}

# Plot residuals  

# Plot log(R/S) residuals ~ observed recruits
layout(mat=matrix(1:8, byrow = TRUE, ncol=4))
par(mfrow=c(2,4), mar=c(4,4,4,0)+0.2)
for(i in 1:length(unique(ests$CU_ID))) {
  dat <- t[t$CU_Name==CUs[i],]
  plot(dat$resid_logRS ~ dat$Recruits, type="p", main=CUs[i])
  abline(a=0, b=1, lty=2, col="orange")
}

