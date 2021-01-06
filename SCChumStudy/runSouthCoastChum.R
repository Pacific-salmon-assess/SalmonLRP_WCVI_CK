# Code by L. Warkentin & K. Holt
# Started on: October 16, 2020

# This file contains the code required to explore South Coast Chum escapement datasets and run retrospective analyses of 
#   data-based LRPs that use logistic regressions

library(rsample)
library(tidyverse)
#library(ggplot2)
#library(tidyr)
#library(dplyr)
options(scipen=1000000)
library(gridExtra)
library(reshape2)
library(TMB)

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

# Load TMB models
compile("TMB_Files/SR_IndivRicker_NoSurv.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_NoSurv"))

compile("TMB_Files/SR_IndivRicker_NoSurv_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_NoSurv_LowAggPrior"))

# Switch to chum directory
setwd(chumDir)

# Run to re-do infilling and brood table
# source("prepare_data.r") 

# ====================================================================
# Read in data and format for using in retrospective analysis
# ====================================================================

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

# ====================================================================================
# Call functions to plot data availability:
# ====================================================================================
plot_CU_DataObs_Over_Time(ChumEscpDat, chumDir, plotName="SC_Chum_DataByCU")
plot_Num_CUs_Over_Time(ChumEscpDat, chumDir, plotName="SC_Chum_N_CUs")

# Note: these next 2 two escpt plots need to have formatting fixed
plot_CU_Escp_Over_Time(ChumEscpDat, chumDir, plotName="SCChum Esc", samePlot = T)
plot_CU_Escp_Over_Time(ChumEscpDat, chumDir, plotName="SCChum Esc Separate", samePlot = F)

# =====================================================================================================================
# Run retrospective analyses:
# =====================================================================================================================

# Get priors for likelihood penalty on aggregate abundance at low probability


TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)

# ========================================================================
# Run annual restrospective analyses using CUs ===========================
# ========================================================================

# Loop over p values and run annual retrospective analyses for each level of p
ps <- 0.90 # for now just use one p value
#ps <- c(seq(0.6, 0.95,.05), 0.99) 

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

# not sure if this step is needed
ChumEscpDat <- ChumEscpDat[ChumEscpDat$yr >= 1958 ,] # remove years without full recruitment
ChumSRDat <- ChumSRDat[ChumSRDat$BroodYear >= 1958 ,] # remove years without full recruitment

source(paste0(codeDir, "/LRPFunctions.r"))


# Run retrospective analysis
for(pp in 1:length(ps)){
  # Run with Binomial LRP model with individual model Ricker
  runAnnualRetro(EscpDat=ChumEscpDat, SRDat=ChumSRDat, startYr=1970, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
                 BMmodel = "SR_IndivRicker_NoSurv", LRPmodel="BinLogistic", integratedModel=T,
                 useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_",ps[pp]*100, sep=""),
                 bootstrapMode = F, plotLRP=T)
  # 
  # # Run with Bernoulli LRP model with individual model Ricker 
  # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurv_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
}


# Get values for low aggregate likelihood penalty model
# Idea is to parameterize (mu and sigma) the aggregate abundnace associated with a very low proportion (essentially 0; p=0.01) of
# CUs being above their benchmark. The idea is to have 95% of this estimate between the lowest CU abundance (e.g., all CUs are gone
# except one) and the sum of all their benchmarks (e.g., all CUs are just below their lower benchmarks). The mean of this
# distribution is halfway between these lower and upper values.

# first get mean size of smallest CU
low_lim <- ChumSRDat %>% group_by(CU_Name) %>% summarise(mean_by_CU = mean(Spawners, na.rm=TRUE)) %>% pull(mean_by_CU) %>% min(., na.rm=TRUE) # FLAG: check that I should use spawners here 
# get sum of CU benchmarks
# get Sgen estimates from running integrated model without penalty
ests <- read.csv("DataOut/AnnualRetrospective/Bin.IndivRicker_NoSurv_90/annualRetro__SRparsByCU.csv", stringsAsFactors = FALSE)
ests1 <- ests[ests$retroYr ==max(ests$retroYr),] # get just last retro year for estimates
hi_lim <- sum(ests1$est_Sgen) # sum Sgen estimates to get upper limit for penalty mu
# make mu of penalty value mean of these two values
B_penalty_mu <- mean(c(low_lim, hi_lim))/TMB_Inputs_IM$Scale
# get sd for prior penalty
sum( dnorm( seq( low_lim, hi_lim, 1), mean=mean(c(low_lim, hi_lim)), sd = 56800)) # FLAG: Should give 95% density.
# plot to check
plot( x = seq( 0, max(ChumSRDat$Spawners), 100), y=dnorm( seq(0, max(ChumSRDat$Spawners), 100), mean=mean(c(low_lim, hi_lim)), sd=56800), xlim=c(0, max(ChumSRDat$Spawners)), type="l")
abline(v=c(low_lim, hi_lim, mean(c(low_lim, hi_lim)))) # plot upper and lower values and mean
B_penalty_sigma <- 56800/TMB_Inputs_IM$Scale # FLAG: This gives weird results, logistic curve flipped 



# Run retrospective analysis with likelihood penalty on aggregate abundance at low proportion CUs above benchmark 
#       (to bring logistic curve intercept down)

for(pp in 1:length(ps)){
  # Run with Binomial LRP model with individual model Ricker
  runAnnualRetro(EscpDat=ChumEscpDat, SRDat=ChumSRDat, startYr=1970, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
                 BMmodel = "SR_IndivRicker_NoSurv_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
                 useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_LowAggPrior",ps[pp]*100, sep=""),
                 bootstrapMode = F, plotLRP=T, B_penalty_mu = B_penalty_mu, B_penalty_sigma = B_penalty_sigma)
  # 
  # # Run with Bernoulli LRP model with individual model Ricker 
  # runAnnualRetro(EscpDat=CoEscpDat, SRDat=CoSRDat, startYr=2015, endYr=2018, BroodYrLag=2, genYrs=3, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_Surv", LRPmodel="BernLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=cohoDir, RunName = paste("Bern.IndivRickerSurv_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
}

# ----------------#
# Examine output  
# ----------------#
# Plot ricker, SMSY, Sgen estimates from integrated model
ests <- read.csv("DataOut/AnnualRetrospective/Bin.IndivRicker_NoSurv_90/annualRetro__SRparsByCU.csv", stringsAsFactors = FALSE)
ests1 <- ests[ests$retroYr ==max(ests$retroYr),] # get just one retro year for estimates
t <- merge(ests1, ChumSRDat[, names(ChumSRDat) %in% c("BroodYear", "Spawners", "Recruits", "CU_Name")], by=c("CU_Name"))

CUs <- unique(t$CU_Name)

# Plot recruits~spawner with ricker and Sgen and SMSY
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

# Plot log(recruits/spawner)~spawner with linear ricker and Sgen and SMSY
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

# get residuals
t$pred_rec <- t$Spawners * t$est_A * exp(-t$est_B * t$Spawners) 
t$resid_rec <- t$Recruits - pred_rec

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
  abline(a=0, b=1, lty=2, col="orange")
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
