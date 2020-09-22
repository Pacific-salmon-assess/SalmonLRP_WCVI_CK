# Make simplified version of TMB run code for demonstration
library(TMB)
library(dplyr)
library(zoo)

sourceAll <- function(){
  source("benchmarkFunctions.r")
  source("LRPFunctions.r")
  source("plotFunctions.r")
  source("retroFunctions.r")
  source("helperFunctions.r")
}
sourceAll()

codeDir<-getwd()
setwd('..')
tmpDir<-getwd()
cohoDir<-paste(tmpDir,"/IFCohoStudy",sep="")
setwd(cohoDir)
CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")
setwd(codeDir)

# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"

CoSRDat$yr_num <- group_indices(CoSRDat, BroodYear) - 1
CoSRDat$CU_ID <- group_indices(CoSRDat, CU_ID) - 1
SRDat <- CoSRDat
EscDat <- CoEscpDat %>%  right_join(unique(SRDat[,c("CU_ID", "CU_Name")]))

# Base LRP model ===============================================================

# only need to compile if changed model
# dyn.unload(dynlib("TMB_Files/Aggregate_LRPs"))
# compile("TMB_Files/Aggregate_LRPs.cpp")
dyn.load(dynlib("TMB_Files/Aggregate_LRPs"))

# want to scale down all data so close to 0
Scale <- 1000
Bern_Logistic <- TRUE
genYrs <- 3

data <- list()
data$S <- SRDat$Spawners/Scale 
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Bern_Logistic <- as.numeric(Bern_Logistic)
data$Sgen_sig <- 0.1

# also give model year for which will fit logistic model
# only give S_Compare values when all obs for all CUs
# which years do we have obs for all CUs?
Num_CUs_Over_Time <- EscDat %>%  filter(is.na(Escp) == F) %>% group_by(yr) %>% 
                     summarise(n=length(unique((CU_Name))))
Mod_Yrs <- Num_CUs_Over_Time$yr[Num_CUs_Over_Time$n == max(Num_CUs_Over_Time$n)]
Logistic_Dat <- EscDat %>% filter(yr %in% Mod_Yrs)
# need year as index
Logistic_Dat$yr_num <- group_indices(Logistic_Dat, yr) - 1
# use trigger for whether or not use generational means outside
# of TMB to keep TMB model as simple as possible
Agg_Abund <- Logistic_Dat %>% group_by(yr) %>% summarise(Agg_Esc = sum(Escp)) %>%
             mutate(Gen_Mean = rollapply(Agg_Esc, genYrs, gm_mean, fill = NA, align="right"))
# add data for logistic model to data list
data$LM_S <- Logistic_Dat$Escp / Scale
data$LM_Agg_Abund <- Agg_Abund$Agg_Esc / Scale
data$LM_yr <- Logistic_Dat$yr_num
data$LM_stk <- Logistic_Dat$CU_ID

# range of agg abund to predict from
data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
# set probability/proportion value used to calc LRP
data$p <- 0.95

param <- list()
param$logA <- rep(1, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
param$B_0 <- 2
param$B_1 <- 0.1

# Phase 1 estimate SR params
# use map to "phase" parameter estimates -- will not fit Sgen, B0, B1 in first round
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
obj <- MakeADFun(data, param, DLL="Aggregate_LRPs", silent=TRUE, map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj))

# pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
# use SMSYs to set upper bounds on Sgen, leave upper bounds for other params at Inf
upper = c(rep(Inf, 3*N_Stocks), log(SMSYs), Inf, Inf )
# use 1/3 of SMSY as starting Sgen value
pl$logSgen <- log(0.3*SMSYs)

#Phase 2 get Sgen, SMSY etc.
map2 <- list(B_0=factor(NA), B_1=factor(NA))
obj <- MakeADFun(data, pl, DLL="Aggregate_LRPs", silent=TRUE, map=map2)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper )
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
summary(sdreport(obj))

# Phase 3 fit logistic model
# Hold other estimates constant
obj <- MakeADFun(data, pl2, DLL="Aggregate_LRPs", silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
          upper = upper)
  
# Create Table of outputs
All_Ests_Basic <- data.frame(summary(sdreport(obj)))
All_Ests_Basic$Param <- row.names(All_Ests_Basic)

# Basic Hier Model ==========================================================================
# only need to compile if changed model
# dyn.unload(dynlib("TMB_Files/Aggregate_LRPs_Hier"))
# compile("TMB_Files/Aggregate_LRPs_Hier.cpp")

dyn.load(dynlib("TMB_Files/Aggregate_LRPs_Hier"))

TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 2.5, 
                   logMuA_sig = 2, Tau_dist = 0.1, Tau_A_dist = 0.1, 
                   gamma_mean = 0, gamma_sig = 10)

# Add prior values to data inputs
data$logMuA_mean <- TMB_Inputs$logMuA_mean 
data$logMuA_sig <- TMB_Inputs$logMuA_sig
data$Tau_dist <- TMB_Inputs$Tau_dist
data$Tau_A_dist <- TMB_Inputs$Tau_A_dist 

data$Sgen_sig <- 1

# add starting values for parameters
param$logMuA <- 2.5
param$logSigmaA <- 1


# Phase 1 estimate SR params
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
obj <- MakeADFun(data, param, DLL="Aggregate_LRPs_Hier", silent=TRUE, random = "logA", map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj))

# pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
upper = c(rep(Inf, 3*N_Stocks), log(SMSYs), Inf, Inf )
pl$logSgen <- log(0.3*SMSYs)

#Phase 2 get Sgen, SMSY etc.
map2 <- list(B_0=factor(NA), B_1=factor(NA))
obj <- MakeADFun(data, pl, DLL="Aggregate_LRPs_Hier", silent=TRUE, random = "logA", map=map2)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper )
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
summary(sdreport(obj))

# Phase 3 fit logistic model
# Hold other estimates constant
map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
# logMuA and logSigmaA disappear after mask, so need to save here
HyperParams <- data.frame(summary(sdreport(obj))) %>% 
  mutate(Param = row.names(summary(sdreport(obj)))) %>% 
  filter(Param %in% c("logMuA", "logSigmaA"))
obj <- MakeADFun(data, pl2, DLL="Aggregate_LRPs_Hier", silent=TRUE, map=map3)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper)

# Create Table of outputs
All_Ests_Hier <- data.frame(summary(sdreport(obj)))
All_Ests_Hier$Param <- row.names(All_Ests_Hier)
All_Ests_Hier <- full_join(All_Ests_Hier, HyperParams)


# Hierarchical RPA model ============================================================

# only need to compile if changed model
# dyn.unload(dynlib("TMB_Files/Aggregate_LRPs_Hier_Surv"))
 compile("TMB_Files/Aggregate_LRPs_Hier_Surv.cpp")

dyn.load(dynlib("TMB_Files/Aggregate_LRPs_Hier_Surv"))

# add proportions of age 3 recruits, survival to data inputs
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
data$logSurv_3 <- log(SRDat$STAS_Age_3)
data$logSurv_4 <- log(SRDat$STAS_Age_4)

# add prior on gamma
data$gamma_mean <- TMB_Inputs$gamma_mean 
data$gamma_sig <- TMB_Inputs$gamma_sig

# Add mean survival to data inputs
muSurv <- SRDat %>% group_by(CU_ID) %>% 
  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
data$muLSurv <- log(muSurv$muSurv)

# add starting values for parameters
param$gamma <- 0


# Phase 1 estimate SR params
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
obj <- MakeADFun(data, param, DLL="Aggregate_LRPs_Hier_Surv", silent=TRUE, random = "logA", map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj))

# pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
upper = c(rep(Inf, 3*N_Stocks), log(SMSYs), Inf, Inf , Inf, Inf)
pl$logSgen <- log(0.3*SMSYs)

#Phase 2 get Sgen, SMSY etc.
map2 <- list(B_0=factor(NA), B_1=factor(NA))
obj <- MakeADFun(data, pl, DLL="Aggregate_LRPs_Hier_Surv", silent=TRUE, random = "logA", map=map2)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper )
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
summary(sdreport(obj))

# Phase 3 fit logistic model
# Hold other estimates constant
map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
# logMuA and logSigmaA disappear after mask, so need to save here
HyperParams <- data.frame(summary(sdreport(obj))) %>% 
  mutate(Param = row.names(summary(sdreport(obj)))) %>% 
  filter(Param %in% c("logMuA", "logSigmaA"))
obj <- MakeADFun(data, pl2, DLL="Aggregate_LRPs_Hier_Surv", silent=TRUE, map=map3)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
          upper = upper)

# Create Table of outputs
All_Ests_Hier_Surv <- data.frame(summary(sdreport(obj)))
All_Ests_Hier_Surv$Param <- row.names(All_Ests_Hier_Surv)
All_Ests_Hier_Surv <- full_join(All_Ests_Hier_Surv, HyperParams)

# Hier Surv Depensation ============================================================

# dyn.unload(dynlib("TMB_Files/Aggregate_LRPs_Hier_Surv_Dep"))
# compile("TMB_Files/Aggregate_LRPs_Hier_Surv_Dep.cpp")

dyn.load(dynlib("TMB_Files/Aggregate_LRPs_Hier_Surv_Dep"))

# add S_dep to data
data$S_dep <- 1000/Scale;

# Phase 1 estimate SR params
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
obj <- MakeADFun(data, param, DLL="Aggregate_LRPs_Hier_Surv_Dep", silent=TRUE, random = "logA", map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj))

# pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
upper = c(rep(Inf, 3*N_Stocks), log(SMSYs), Inf, Inf , Inf, Inf)
pl$logSgen <- log(0.3*SMSYs)

#Phase 2 get Sgen, SMSY etc.
map2 <- list(B_0=factor(NA), B_1=factor(NA))
obj <- MakeADFun(data, pl, DLL="Aggregate_LRPs_Hier_Surv_Dep", silent=TRUE, random = "logA", map=map2)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper )
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
summary(sdreport(obj))

# Phase 3 fit logistic model
# Hold other estimates constant
map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
# logMuA and logSigmaA disappear after mask, so need to save here
HyperParams <- data.frame(summary(sdreport(obj))) %>% 
  mutate(Param = row.names(summary(sdreport(obj)))) %>% 
  filter(Param %in% c("logMuA", "logSigmaA"))
obj <- MakeADFun(data, pl2, DLL="Aggregate_LRPs_Hier_Surv_Dep", silent=TRUE, map=map3)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper)

# Create Table of outputs
All_Ests_Hier_Surv_Dep <- data.frame(summary(sdreport(obj)))
All_Ests_Hier_Surv_Dep$Param <- row.names(All_Ests_Hier_Surv_Dep)
All_Ests_Hier_Surv_Dep <- full_join(All_Ests_Hier_Surv_Dep, HyperParams)

# Hier Surv No Priors ============================================================

# This version may not converge -- keeps all alphas same

# only need to compile if changed model
# dyn.unload(dynlib("TMB_Files/Aggregate_LRPs_Hier_Surv_NoPriors"))
# compile("TMB_Files/Aggregate_LRPs_Hier_Surv_NoPriors.cpp")

dyn.load(dynlib("TMB_Files/Aggregate_LRPs_Hier_Surv_NoPriors"))


# Add prior values to data inputs
data$logMuA_mean <- data$logMuA_sig <- data$Tau_dist <- data$Tau_A_dist <- NULL
data$gamma_mean <- data$gamma_sig <- data$S_dep <- NULL

# Phase 1 estimate SR params
map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
obj <- MakeADFun(data, param, DLL="Aggregate_LRPs_Hier_Surv_NoPriors", silent=TRUE, random = "logA", map=map)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
summary(sdreport(obj))

# pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
upper = c(rep(Inf, 3*N_Stocks), log(SMSYs), Inf, Inf , Inf, Inf)
pl$logSgen <- log(0.3*SMSYs)

#Phase 2 get Sgen, SMSY etc.
map2 <- list(B_0=factor(NA), B_1=factor(NA))
obj <- MakeADFun(data, pl, DLL="Aggregate_LRPs_Hier_Surv_NoPriors", silent=TRUE, random = "logA", map=map2)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper )
pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
summary(sdreport(obj))

# Phase 3 fit logistic model
# Hold other estimates constant
map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
# logMuA and logSigmaA disappear after mask, so need to save here
HyperParams <- data.frame(summary(sdreport(obj))) %>% 
  mutate(Param = row.names(summary(sdreport(obj)))) %>% 
  filter(Param %in% c("logMuA", "logSigmaA"))
obj <- MakeADFun(data, pl2, DLL="Aggregate_LRPs_Hier_Surv_NoPriors", silent=TRUE, map=map3)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
              upper = upper)

# Create Table of outputs
All_Ests_Hier_Surv_NP <- data.frame(summary(sdreport(obj)))
All_Ests_Hier_Surv_NP$Param <- row.names(All_Ests_Hier_Surv_NP)
All_Ests_Hier_Surv_NP <- full_join(All_Ests_Hier_Surv_NP, HyperParams)
