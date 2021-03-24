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

# ====================================================================
# Read in data and format for using in retrospective analysis
# ====================================================================
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

# ====================================================================================
# Call functions to plot data availability (note that this is for infilled data):
# ====================================================================================
plot_CU_DataObs_Over_Time(ChumEscpDat, chumDir, plotName="SC_Chum_DataByCU")
plot_Num_CUs_Over_Time(ChumEscpDat, chumDir, plotName="SC_Chum_N_CUs")

# Note: these next 2 two escpt plots need to have formatting fixed
plot_CU_Escp_Over_Time(ChumEscpDat, chumDir, plotName="SCChum Esc", samePlot = T)
plot_CU_Escp_Over_Time(ChumEscpDat, chumDir, plotName="SCChum Esc Separate", samePlot = F)

# =====================================================================================================================
# Setup for retrospective analysis
# =====================================================================================================================

# Data prep

# FLAG: Should probably limit stock-recruit data to year > 1959/1960 to allow for full brood year returns up to age 6. 
# This may be done automatically, see retroFunctions.r line 20 
# This is done automatically using the BroodYrLag variable 
# remove years without full recruitment  - FLAG: still needing to run this, otherwise have NA recrtuiment going into model optimization
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
# Spawners = CU-infilled years + max age of fish
# Recruits = CU-infilled years - max age of fish
# max age of fish is 6 years # FLAG change to 3-6 yrs. 
# make vector of years not to use, all cu-level infilled years and all the years 
# Make vector that includes the 6 years before and after each CU-level infilled year
years_not_use2 <- unique(
                     as.vector(
                        sapply( years_not_use, function(x) {
                                 yf <- x - 6
                                 yl <- x + 6
                                 y <- yf:yl 
                                 y })))
ChumEscpDat_no_CU_infill_yrs <- ChumEscpDat[!(ChumEscpDat$yr %in% years_not_use2), ]
ChumSRDat_no_CU_infill_yrs <- ChumSRDat[!(ChumSRDat$BroodYear %in% years_not_use2), ]
range(ChumEscpDat_no_CU_infill_yrs$yr) 
# Note that there are only 15 years of data after removing
# years with CU-level infilling. May make it not feasible to
# check effect. 

# Specify p value for logistic regression
ps <- c(seq(0.6, 0.95,.05), 0.99) 

# ========================================================================
# Run annual restrospective analysis using stock recruit parameter based Sgen
# ========================================================================

TMB_Inputs_IM <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)

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
  # runAnnualRetro(EscpDat=ChumEscpDat, SRDat=ChumSRDat, startYr=1970, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_NoSurv", LRPmodel="BinLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
  # Run with CUs with CU-level infilling removed
  runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill, SRDat=ChumSRDat_no_CU_infill, startYr=1967, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
                BMmodel = "SR_IndivRicker_NoSurv", LRPmodel="BinLogistic", integratedModel=T,
                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_noCUinfill_",ps[pp]*100, sep=""),
                bootstrapMode = F, plotLRP=T)
  
  # Run with years with CU-level infilling removed
  runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill_yrs, SRDat=ChumSRDat_no_CU_infill_yrs, startYr=1967, endYr=1972, BroodYrLag=4, genYrs=4, p = ps[pp],
                 BMmodel = "SR_IndivRicker_NoSurv", LRPmodel="BinLogistic", integratedModel=T,
                 useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_noCUinfill_yrs_",ps[pp]*100, sep=""),
                 bootstrapMode = F, plotLRP=T)
}

# -----------------------------------------------------#
# Get values for low aggregate likelihood penalty model 
# FLAG: need to redo this for the 5 CUs without CU-level infilling.
# -----------------------------------------------------#
# Idea is to parameterize (mu and sigma) the aggregate abundance 
# associated with a very low proportion (essentially 0; p=0.01) of
# CUs being above their benchmark. The idea is to have 95% of this 
# estimate between the CU Sgen and the sum of all their Sgen 
# benchmarks or upper estimates of Sgens (e.g., all CUs are just 
# below their lower benchmarks). The mean of this distribution is 
# halfway between these lower (2.5% quantile) and upper (97.5% quantile) values. 

# get Sgen estimates from running integrated model without penalty
ests <- read.csv("DataOut/AnnualRetrospective/Bin.IndivRicker_NoSurv_noCUinfill_90/annualRetro_SRparsByCU.csv", stringsAsFactors = FALSE)

ests1 <- ests[ests$retroYr ==max(ests$retroYr),] # get just last retro year for estimates
# make lower limit the lowest CU Sgen (this gives essentially the same results as 
# using the average abundance of the smallest CU)
low_lim <- min(ests1$est_Sgen) 
# make upper limit the sum of CU benchmarks
#hi_lim <- sum(ests1$est_Sgen) # sum Sgen estimates to get upper limit for penalty mu
# sum upper CI of Sgen estimates to get upper limit for penalty mu
hi_lim <- sum(na.omit(ests1$up_Sgen)) 


# make mu of penalty value mean of these two values, divide by scale
B_penalty_mu <- mean(c(low_lim, hi_lim))
# get SD that gives 95% density between lower and upper limits (using getSD helper function)
dum<-optim(par=200, fn = getSD, method="Brent",lower=1, upper=50000000, low_lim=low_lim, hi_lim=hi_lim)
B_penalty_sigma<-dum$par

# get SD for prior penalty so that 95% density is between lower and upper limits
# plot to check
plot( x = seq( 0, max(ChumSRDat$Spawners), 100), y=dnorm( seq(0, max(ChumSRDat$Spawners), 100), mean=mean(c(low_lim, hi_lim)), sd=B_penalty_sigma), 
      xlim=c(0, max(ChumSRDat$Spawners)), type="l", ylab="density", xlab="aggregate adundance")
abline(v=c(low_lim, hi_lim, mean(c(low_lim, hi_lim))), col="dodgerblue", lty=c(2,2,1)) # plot upper and lower values and mean

# Add penalty values to TMB_Inputs_IM
TMB_Inputs_IM_LowAggPrior <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                      B_penalty_mu = B_penalty_mu, B_penalty_sigma = B_penalty_sigma)


# Section below is obsolete
#B_penalty_sigmas <- c(B_penalty_sigma, B_penalty_sigma* 1.5, B_penalty_sigma * 2) # test with 1.5 and 2 times SD
# Effect of increasing SD:
#   Without CUs with CU infilling:
#       SD times 1.5 and 2 make very minor differences in LRP (very small increases). No flipping of logistic curve
#   With using CUs with CU infilling:
#    Using 1.5 times SD makes logistic curve flip in one year, and 2 times SD makes it flip several years (with p=0.9)

# -------------------------------------#
# Run retrospective analysis with likelihood penalty on aggregate abundance at low proportion CUs above benchmark 
#       (to bring logistic curve intercept down)
# -------------------------------------#

for(pp in 1:length(ps)){
  #for (i in 1:length(B_penalty_sigmas)){ # for testing sensitivity to B_penalty_sigma
  # Run with Binomial LRP model with individual model Ricker
  # runAnnualRetro(EscpDat=ChumEscpDat, SRDat=ChumSRDat, startYr=2010, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
  #                BMmodel = "SR_IndivRicker_NoSurv_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_IM, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_LowAggPrior_", ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T, B_penalty_mu = B_penalty_mu, B_penalty_sigma = B_penalty_sigma)
  # Without CU-level infilling
  runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill, SRDat=ChumSRDat_no_CU_infill, startYr=1967, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
                BMmodel = "SR_IndivRicker_NoSurv_LowAggPrior", LRPmodel="BinLogistic", integratedModel=T,
                useGenMean=F, TMB_Inputs=TMB_Inputs_IM_LowAggPrior, outDir=chumDir, RunName = paste("Bin.IndivRicker_NoSurv_LowAggPrior_noCUinfill_", ps[pp]*100, sep=""),
                bootstrapMode = F, plotLRP=T )
  #}
}

# ---------------------------#
# Percentile benchmarks
# ---------------------------#

# Currently there is a discrepancy between CU_ID variable in percentile/subpop and Sgen 
# based retrospective. Need to sort this out. 

TMB_Inputs_Percentile <- list(Scale = 1000, logA_Start = 1,
                              Tau_dist = 0.1,
                              gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                              B_penalty_mu = NA, B_penalty_sigma = NA)

# Run retrospective analysis using percentile benchmarks
for(pp in 1:length(ps)){
  # Run with Binomial LRP with CUs with CU-level infilling removed
  runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill, SRDat=ChumSRDat_no_CU_infill, startYr=1967, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
                 BMmodel = "Percentile", LRPmodel="BinLogistic", LRPfile="LRP_Logistic_Only",integratedModel=F,
                 useGenMean=F, TMB_Inputs=TMB_Inputs_Percentile, outDir=chumDir, RunName = paste("Bin.Percentile_noCUinfill_",ps[pp]*100, sep=""),
                 bootstrapMode = F, plotLRP=T)
  # with prior penalty on low aggregate abundance
  # runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill, SRDat=ChumSRDat_no_CU_infill, startYr=1970, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
  #                BMmodel = "LRP_Logistic_Only_LowAggPrior", LRPmodel="BinLogistic", integratedModel=F,
  #                useGenMean=F, TMB_Inputs=TMB_Inputs_Percentile, outDir=chumDir, RunName = paste("Bin.Percentile_LowAggPrior_noCUinfill_",ps[pp]*100, sep=""),
  #                bootstrapMode = F, plotLRP=T)
}


# Get low aggregate penalty values
pdat <- read.csv("DataOut/AnnualRetrospective/Bin.Percentile_noCUinfill_60/annualRetro_perc_benchmarks.csv", stringsAsFactors = FALSE)
pdat1 <- pdat[pdat$retro_year ==max(pdat$retro_year),]
pdat2 <- pdat1 %>% group_by(CU_Name, CU, benchmark_perc_25) %>% summarise(n())
# make lower limit=0, upper limit sum of 25% benchmarks for retro year 2010
low_lim_perc <- 0
hi_lim_perc <- sum(pdat2$benchmark_perc_25)

B_penalty_perc_mu <- mean(c(low_lim_perc, hi_lim_perc)) 
dum_perc<-optim(par=200, fn = getSD, method="Brent",lower=1, upper=1000000, low_lim=low_lim_perc, hi_lim=hi_lim_perc)
B_penalty_perc_sigma<-dum_perc$par

TMB_Inputs_Percentile <- list(Scale = 1000, logA_Start = 1,
                              Tau_dist = 0.1,
                              gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1,
                              B_penalty_mu = B_penalty_perc_mu, B_penalty_sigma = B_penalty_perc_sigma)

# Run retrospective analysis using percentile benchmarks, with likelihood penalty
for(pp in 1:length(ps)){
  # with prior penalty on low aggregate abundance
  runAnnualRetro(EscpDat=ChumEscpDat_no_CU_infill, SRDat=ChumSRDat_no_CU_infill, startYr=1967, endYr=2010, BroodYrLag=4, genYrs=4, p = ps[pp],
                 BMmodel = "Percentile", LRPmodel="BinLogistic", LRPfile="LRP_Logistic_Only_LowAggPrior", integratedModel=F,
                 useGenMean=F, TMB_Inputs=TMB_Inputs_Percentile, outDir=chumDir, RunName = paste("Bin.Percentile_LowAggPrior_noCUinfill_",ps[pp]*100, sep=""),
                 bootstrapMode = F, plotLRP=T)
}


# ----------------#
# Examine output  
# ----------------#

# Plot 25% benchmark over time
pdat <- read.csv("DataOut/AnnualRetrospective/Bin.Percentile_noCUinfill_60/annualRetro_perc_benchmarks.csv", stringsAsFactors = FALSE)
png("Figures/fig_perc_benchmarks_annual_retro.png", width=8, height=4, res=300, units="in")
ggplot(pdat, aes(y=benchmark_perc_25, x=retro_year)) +
  geom_line(colour="dodgerblue", lwd=1.2) +
  geom_point(aes(y=Escp, x=yr), shape=1) +
  geom_line(aes(y=Escp, x=yr)) +
  facet_wrap(~CU_Name, scales="free_y") +
  ylab("25% benchmark") + xlab("Retro year") +
  geom_hline(aes(yintercept=0), linetype=2) +
  theme_bw()
dev.off()

# ---------------------------------------------#
# Look at stock recruit parameters over time
# ---------------------------------------------#
# Plot Sgen over time
mdat <- read.csv("DataOut/AnnualRetrospective/Bin.IndivRicker_NoSurv_noCUinfill_95/annualRetro_SRparsByCU.csv", stringsAsFactors = FALSE)
mdat1 <- mdat %>% pivot_longer(cols=est_B:up_Sgen, names_to="param", values_to="est")

# Plot Sgen estimates over time
png("Figures/fig_Sgen_annual_retro.png", width=12, height=3, res=300, units="in")
ggplot(mdat, aes(y=est_Sgen, x=retroYr)) +
  geom_line() +
  geom_ribbon(aes(ymin=low_Sgen, ymax=up_Sgen, x=retroYr),colour=NA, alpha=0.2) +
  facet_wrap(~CU_Name, scales="free_y", nrow=1) +
  ylab("Sgen") + xlab("Year") +
  geom_hline(aes(yintercept=0), linetype=2) +
  theme_bw()
dev.off()

# Plot Alpha and beta on same plot
png("Figures/fig_a_b_annual_retro.png", width=10, height=4, res=300, units="in")
mdat1 %>% filter(param %in% c("est_A", "est_B")) %>%
  ggplot(., aes(y=est, x=retroYr)) +
  geom_line() +
  facet_grid( param~CU_Name, scales="free_y")+
  theme_bw()
dev.off()

# Plot ricker, SMSY, Sgen estimates from integrated model
ests <- read.csv("DataOut/AnnualRetrospective/Bin.IndivRicker_NoSurv_noCUinfill_95/annualRetro_SRparsByCU.csv", stringsAsFactors = FALSE)
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
