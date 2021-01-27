year <- 2018
p<-0.99
useGenMean <- FALSE
genYrs<-3
BroodYrLag<-2
useBern_Logistic <- FALSE


library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(tmbstan)
library(here)

setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
cohoDir<-paste(rootDir,"/IFCohoStudy",sep="")

setwd(codeDir)

sourceAll <- function(){
  source("benchmarkFunctions.r")
  source("LRPFunctions.r")
  source("plotFunctions.r")
  source("retroFunctions.r")
  source("helperFunctions.r")
}
sourceAll()


compile("TMB_Files/SR_IndivRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv"))

compile("TMB_Files/SR_IndivRicker_Surv_LowAggPrior.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_LowAggPrior"))


# ======================================================================
# Read-in Coho data:  
# =====================================================================
setwd(cohoDir)

CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"

CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")

CoEscpDat_bySubpop<-read.csv("DataIn/IFCoho_escpBySubpop.csv")

# Change column names to yr, CU_Name, Escp, Subpop_Name, CU_ID
CoEscpDat_bySubpop<-CoEscpDat_bySubpop %>% select(yr=Return.Year, CU_Name=Conservation.Unit, Escp=Natural.Returns, Subpop_Name=Sub.Population)
tmp.df<-data.frame(CU_Name=unique(CoEscpDat_bySubpop$CU_Name), CU_ID=seq(1,length(unique(CoEscpDat_bySubpop$CU_Name)),by=1))
CoEscpDat_bySubpop <- left_join(CoEscpDat_bySubpop,tmp.df)

setwd(codeDir)


# Get data into required format  ============================================================
# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)
CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)
CoEscpDat_bySubpop <- CoEscpDat_bySubpop %>% filter(yr >= 1998)

# Roll up escpaments, and get Gen Mean of htat
AggEscp <- CoEscpDat %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Escp, 3, gm_mean, fill = NA, align="right"))





# Specify model options ===============================================
LRPmodel <- "BinLogistic"


Dat <- CoSRDat %>%  filter(BroodYear <= (year-BroodYrLag))
EscpDat.yy <- CoEscpDat %>% filter(yr <= year) 
Dat$yr_num <- group_by(Dat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
Dat$CU_ID <- group_by(Dat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
EscDat <- EscpDat.yy %>%  right_join(unique(Dat[,c("CU_ID", "CU_Name")]))



figDir<-paste(cohoDir, "Figures", "testLogisticPen", sep="/")
if (file.exists(figDir) == FALSE){
  dir.create(figDir)
} 


# Model 1: No penalty
BMmodel <- "SR_IndivRicker_Surv"


TMB_Inputs <- list(Scale = 1000, logA_Start = 1,
                      Tau_dist = 0.1,
                      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)

# Call LRP model
LRP_Mod <- Run_Ricker_LRP(SRDat = Dat, EscDat = EscDat, BMmodel = BMmodel, Bern_Logistic = useBern_Logistic, 
                          useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs)

plotLogistic(Data = LRP_Mod$Logistic_Data, Preds = LRP_Mod$Preds, 
             LRP = LRP_Mod$LRP, useGenMean = useGenMean,
             plotName = paste("LogisticMod", year, sep ="_"), outDir = figDir,
             p = p, useBern_Logistic = useBern_Logistic)



# Model 2: With Low Agg Abundance Prior

setwd(cohoDir)
# get Sgen estimates from running integrated model without penalty
ests <- read.csv("DataOut/AnnualRetrospective/Bin.IndivRickerSurv_99/annualRetro__SRparsByCU.csv", stringsAsFactors = FALSE)
#ests <- read.csv(paste(cohoDir,"/DataOut/AnnualRetrospective/Bin.IndivRickerSurv_99/annualRetro__SRparsByCU.csv",sep=""), stringsAsFactors = FALSE)
ests1 <- ests[ests$retroYr ==max(ests$retroYr),] # get just last retro year for estimates
# make lower limit the lowest CU Sgen (this gives essentially the same results as using the average abundance of the smallest CU)
low_lim <- min(na.omit(ests1$est_Sgen)) 
# make upper limit the sum of CU benchmarks
hi_lim <- sum(na.omit(ests1$est_Sgen)) # sum Sgen estimates to get upper limit for penalty mu

# make mu of penalty value mean of these two values, divide by scale
B_penalty_mu <- mean(c(low_lim, hi_lim))

# Function to get SD for prior penalty so that 95% density is between lower and upper limits
getSD<-function(par, low_lim,hi_lim) {
  den<-sum( dnorm( seq( low_lim, hi_lim, 1), mean=mean(c(low_lim, hi_lim)), sd = par))
  return(abs(den - 0.95)) 
}

dum<-optim(par=200, fn = getSD, method="Brent",lower=1, upper=5000, low_lim=low_lim, hi_lim=hi_lim)
B_penalty_sigma<-dum$par


# Plot distribution of penalty
#hist(rnorm(1000,B_penalty_mu, B_penalty_sigma))

# plot to check
plot( x = seq( 0, max(CoSRDat$Spawners), 100), y=dnorm( seq(0, max(CoSRDat$Spawners), 100), mean=mean(c(low_lim, hi_lim)), sd=B_penalty_sigma), 
      xlim=c(0, max(CoSRDat$Spawners)), type="l", ylab="density", xlab="aggregate adundance")
abline(v=c(low_lim, hi_lim, mean(c(low_lim, hi_lim))), col="dodgerblue", lty=c(2,2,1)) # plot upper and lower values and mean

TMB_Inputs_LowAggPrior <- list(Scale = 1000, logA_Start = 1,
                   Tau_dist = 0.1,
                   gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1, 
                   B_penalty_mu=B_penalty_mu/1000, B_penalty_sigma = B_penalty_sigma / 1000)


BMmodel <- "SR_IndivRicker_Surv_LowAggPrior"

LRP_Mod_LowAggPrior <- Run_Ricker_LRP(SRDat = Dat, EscDat = EscDat, BMmodel = BMmodel, Bern_Logistic = useBern_Logistic, 
                          useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs=TMB_Inputs_LowAggPrior)

plotLogistic(Data = LRP_Mod_LowAggPrior$Logistic_Data, Preds = LRP_Mod_LowAggPrior$Preds, 
             LRP = LRP_Mod_LowAggPrior$LRP, useGenMean = useGenMean,
             plotName = paste("LogisticMod_LowAggPrior", year, sep ="_"), outDir = figDir,
             p = p, useBern_Logistic = useBern_Logistic)



# Model 3: Low hi_Lim

B_penalty_sigma<-B_penalty_mu*5

# Plot distribution of penalty
#hist(rnorm(1000,B_penalty_mu, B_penalty_sigma))

# plot to check
plot( x = seq( 0, max(CoSRDat$Spawners), 100), y=dnorm( seq(0, max(CoSRDat$Spawners), 100), mean=mean(c(low_lim, hi_lim)), sd=B_penalty_sigma), 
      xlim=c(0, max(CoSRDat$Spawners)), type="l", ylab="density", xlab="aggregate adundance")
abline(v=c(low_lim, hi_lim, mean(c(low_lim, hi_lim))), col="dodgerblue", lty=c(2,2,1)) # plot upper and lower values and mean

TMB_Inputs_LowAggPrior <- list(Scale = 1000, logA_Start = 1,
                               Tau_dist = 0.1,
                               gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1, 
                               B_penalty_mu=B_penalty_mu/1000, B_penalty_sigma = B_penalty_sigma / 1000)


BMmodel <- "SR_IndivRicker_Surv_LowAggPrior"

LRP_Mod_LowAggPrior_hiSig <- Run_Ricker_LRP(SRDat = Dat, EscDat = EscDat, BMmodel = BMmodel, Bern_Logistic = useBern_Logistic, 
                                      useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs=TMB_Inputs_LowAggPrior)

plotLogistic(Data = LRP_Mod_LowAggPrior_hiSig$Logistic_Data, Preds = LRP_Mod_LowAggPrior_hiSig$Preds, 
             LRP = LRP_Mod_LowAggPrior_hiSig$LRP, useGenMean = useGenMean,
             plotName = paste("LogisticMod_LowAggPrior_hiSig", year, sep ="_"), outDir = figDir,
             p = p, useBern_Logistic = useBern_Logistic)

# Create data frame for making LRP comparison plots





modelNames<-c("No Pen","Pen SD = 95%CI", "Pen CV = 500%")

for (i in 1:length(modelNames)) {

  if (i == 1) Mod<-LRP_Mod
  if (i == 2) Mod<-LRP_Mod_LowAggPrior
  if (i == 3) Mod<-LRP_Mod_LowAggPrior_hiSig
  
  df<-Mod$LRP
  df$Name<-modelNames[i]
  
  if (i == 1) {
    LRP.df<-df }
  else {
    LRP.df<-rbind(LRP.df,df)
  }
  
}

LRP.df$Name<-factor(LRP.df$Name, levels = c("No Pen","Pen SD = 95%CI", "Pen CV = 500%"))

g<-ggplot(data=LRP.df, mapping=aes(x=Name, y=fit)) + geom_point() +
  geom_errorbar(mapping=aes(x=Name, ymax=upr, ymin=lwr), width=0) +
  xlab("") + ylab("Aggregate LRP") +
  theme_classic()




