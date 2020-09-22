# Code used to fit intergrated LRP models using MCMC in stan

# For stantmb example, see: https://github.com/kaskr/tmbstan#examples


nIter<-100000

library(TMB)
library(dplyr)

# read in data
priorDir<-getwd()
setwd('..')
cohoDir<-getwd()
setwd('..')
tmpDir<-getwd()
codeDir<-paste(tmpDir,"/Code",sep="")
setwd(cohoDir)
CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")
setwd(codeDir)

compile("TMB_Files/SR_HierRicker_Surv.cpp")
dyn.load(dynlib("TMB_Files/SR_HierRicker_Surv"))


sourceAll <- function(){
  source("benchmarkFunctions.r")
  source("LRPFunctions.r")
  source("plotFunctions.r")
  source("retroFunctions.r")
  source("helperFunctions.r")
}

sourceAll()



# Change header names to match generic data headers (this will allow generic functions from Functions.r to be used)
colnames(CoEscpDat)[colnames(CoEscpDat)=="CU_ID"] <- "CU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="MU_Name"] <- "MU"
colnames(CoEscpDat)[colnames(CoEscpDat)=="ReturnYear"] <- "yr"
colnames(CoEscpDat)[colnames(CoEscpDat)=="Escapement"] <- "Escp"


# Restrict data set to years 1998+ based on recommendation from Michael Arbeider
CoEscpDat <- CoEscpDat %>% filter(yr >= 1998)
CoSRDat <- CoSRDat %>% filter(BroodYear >= 1998)

# Prep data frame
CoSRDat$yr_num <- CoSRDat$BroodYear - min(CoSRDat$BroodYear)
CoSRDat$CU_ID <- group_indices(CoSRDat, CU_ID) - 1

CoEscpDat$yr_num <- group_indices(CoEscpDat, BroodYear) - 1
CoEscpDat<- CoEscpDat %>% right_join(unique(CoSRDat[,c("CU_ID", "CU_Name")]))


SRDat<-CoSRDat
EscDat <- CoEscpDat



Bern_Logistic <- TRUE
useGenMean <- FALSE
p <- 0.95
genYrs<-3


# *************************************************************************************
# Fit Hier_Ricker_Surv model
# ********************************************************************************

TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1,
                   logMuA_sig = 2, Tau_dist = 0.01, Tau_A_dist = 0.1,
                   gamma_mean = 0, gamma_sig = 100, S_dep = 1000, Sgen_sig = 0.5)

Mod <- "SR_HierRicker_Surv"

# Set-up for call to TMB
Scale <- TMB_Inputs$Scale

data <- list()
data$S <- SRDat$Spawners/Scale
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Bern_Logistic <- as.numeric(Bern_Logistic)

# also give model year for which will fit logistic model
# only give S_Compare values when all obs for all CUs
# give stocks associated with S_Compare
# which years do we have obs for all CUs?
Num_CUs_Over_Time <- EscDat %>%  filter(is.na(Escp) == F) %>% group_by(yr) %>% summarise(n=length(unique((CU_Name))))
Mod_Yrs <- Num_CUs_Over_Time$yr[Num_CUs_Over_Time$n == max(Num_CUs_Over_Time$n)]
Logistic_Dat <- EscDat %>% filter(yr %in% Mod_Yrs)

Agg_Abund <- Logistic_Dat %>% group_by(yr) %>% summarise(Agg_Esc = sum(Escp)) %>%
  mutate(Gen_Mean = rollapply(Agg_Esc, genYrs, gm_mean, fill = NA, align="right"))

# need year as index
Logistic_Dat$yr_num <- group_indices(Logistic_Dat, yr) - 1

if(useGenMean == T){
  GenMean_DF <- EscDat %>% group_by(CU_ID) %>% mutate(Gen_Mean = rollapply(Escp, genYrs, gm_mean, fill = NA, align="right")) %>%
    filter(is.na(Gen_Mean) == F) %>% ungroup() %>% mutate(yr_num = yr - min(yr))
  data$LM_S <- GenMean_DF$Gen_Mean / Scale
  data$LM_yr <- GenMean_DF$yr_num
  data$LM_stk <- GenMean_DF$CU_ID
  Agg_Abund <- Agg_Abund %>% filter(is.na(Gen_Mean) == F)
  data$LM_Agg_Abund <- Agg_Abund$Gen_Mean / Scale
} else {
  data$LM_S <- Logistic_Dat$Escp / Scale
  data$LM_Agg_Abund <- Agg_Abund$Agg_Esc / Scale
  data$LM_yr <- Logistic_Dat$yr_num
  data$LM_stk <- Logistic_Dat$CU_ID
}
# range of agg abund to predict from
data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
# range of spawner abundance to predict recruitment from
data$Pred_Spwn <- rep(seq(0,max(data$S)*1.1,length=100), N_Stocks)
for (i in 1:N_Stocks) {
  if (i == 1) {
    data$stk_predS <- rep(unique(SRDat$CU_ID)[1],100)
  } else {
    data$stk_predS<-c(data$stk_predS,rep(unique(SRDat$CU_ID)[i],100))
  }
}

data$p <- p

# set variance to be used for likelihood for estimating Sgen
data$Sgen_sig <- TMB_Inputs$Sgen_sig

# set-up init params
param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale)
param$B_0 <- 2
param$B_1 <- 0.1
param$logMuA <- TMB_Inputs$logA_Start
param$logSigmaA <- 1
param$gamma <- 0

# specify data
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
data$logSurv_3 <- log(SRDat$STAS_Age_3)
data$logSurv_4 <- log(SRDat$STAS_Age_4)

muSurv <- SRDat %>% group_by(CU_ID) %>%
  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
data$muLSurv <- log(muSurv$muSurv)

data$Tau_dist <- TMB_Inputs$Tau_dist
data$gamma_mean <- TMB_Inputs$gamma_mean
data$gamma_sig <- TMB_Inputs$gamma_sig
data$logMuA_mean <- TMB_Inputs$logMuA_mean
data$logMuA_sig <- TMB_Inputs$logMuA_sig
data$Tau_A_dist <- TMB_Inputs$Tau_A_dist


# Phase 1 estimate SR params ============
  map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)

  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
  pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1


  # -- pull out SMSY values from phase 1
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

  pl$logSgen <- log(0.3*SMSYs)


# Phase 2 get Sgen, SMSY etc. =================
  map2 <- list(B_0=factor(NA), B_1=factor(NA))
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA", map=map2)


   # Create upper bounds vector that is same length and order as start vector that will be given to nlminb
  upper<-unlist(obj$par)
  upper[1:length(upper)]<-Inf
  upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
  upper<-unname(upper)

  lower<-unlist(obj$par)
  lower[1:length(lower)]<--Inf
  lower[names(lower) =="logSgen"] <- log(0.001)
  lower<-unname(lower)


  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
                upper = upper, lower=lower)
  pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2


  # -- pull out SMSY values from phase 2
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

# Phase 3 fit logistic model ======================
# Hold other estimates constant

  map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
  # logMuA and logSigmaA disappear after mask, so need to save here
  HyperParams <- data.frame(summary(sdreport(obj))) %>%
    mutate(Param = row.names(summary(sdreport(obj)))) %>%
    filter(Param %in% c("logMuA", "logSigmaA"))
  obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE, map=map3)

  # Create upper bounds vector that is same length and order as start vector that will be given to nlminb
  upper<-unlist(obj$par)
  upper[1:length(upper)]<-Inf
  upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
  upper<-unname(upper)

  lower<-unlist(obj$par)
  lower[1:length(lower)]<--Inf
  lower[names(lower) =="logSgen"] <- log(0.001)
  lower<-unname(lower)

  opt <- tryCatch(
    {nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
            upper = upper, lower=lower)},
    error=function(cond) {
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      return(NA)
    }
  )


  # -- pull out SMSY values from phase 3
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

  # Phase 4 do it all ....

  pl3<-obj$env$parList(opt$par)

  obj <- MakeADFun(data, pl3, DLL=Mod, silent=TRUE)

  # Create upper bounds vector that is same length and order as start vector that will be given to nlminb
  upper<-unlist(obj$par)
  upper[1:length(upper)]<-Inf
  upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
  upper<-unname(upper)

  lower<-unlist(obj$par)
  lower[1:length(lower)]<--Inf
  lower[names(lower) =="logSgen"] <- log(0.001)
  lower<-unname(lower)


  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
                upper = upper, lower = lower)

  # Fit mcmc with STAN ==================================================================

  library(tmbstan)
  library(rstan)
  
  fitmcmc <- tmbstan(obj, chains=3, iter=nIter, init=opt$par,
                         control = list(adapt_delta = 0.95),upper=upper, lower=lower)


  ## Can get ESS and Rhat from rstan::monitor
  mon <- monitor(fitmcmc)
  max(mon$Rhat)
  min(mon$Tail_ESS)

  write.csv(mon,paste(priorDir,"/mcmcDiag_Hier_Ricker_Surv.csv", sep=""))

  # Marginal posteriors for Sgen
  png(filename = paste(priorDir,"/postHist_Sgens_Hier_Ricker_Surv.png",sep=""))
  par(mfrow=c(3,2))

  # pull out posterior vals
  post<-as.matrix(fitmcmc)

  hist(exp(post[,'logSgen[1]'])*Scale)
  hist(exp(post[,'logSgen[2]'])*Scale)
  hist(exp(post[,'logSgen[3]'])*Scale)
  hist(exp(post[,'logSgen[4]'])*Scale)
  hist(exp(post[,'logSgen[5]'])*Scale)

  dev.off()



  ## Trace plot
  parNames<-names(fitmcmc)[-length(names(fitmcmc))]
  png(filename = paste(priorDir,"/tracePlot_Ricker_Hier_Surv.png",sep=""))
  traceplot(fitmcmc, pars=parNames, inc_warmup=TRUE)
  dev.off()






#   # *************************************************************************************
#   # Fit Hier_Ricker_Surv model - with SR only (i.e., logistic regression removed)
#   # ********************************************************************************
#   
#   TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
#                      logMuA_sig = 2, Tau_dist = 0.01, Tau_A_dist = 0.1, 
#                      gamma_mean = 0, gamma_sig = 100, 
#                      #Sgen_sig = (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/1000*0.02)
#                      Sgen_sig = rep(0.2,5))
#   Mod <- "Aggregate_LRPs_Hier_Surv_onlySR"
#   
#   
#   # Set-up for call to TMB
#   Scale <- TMB_Inputs$Scale
#   
#   data <- list()
#   data$S <- SRDat$Spawners/Scale 
#   data$logR <- log(SRDat$Recruits/Scale)
#   data$stk <- as.numeric(SRDat$CU_ID)
#   N_Stocks <- length(unique(SRDat$CU_Name))
#   data$N_Stks <- N_Stocks
#   data$yr <- SRDat$yr_num
#   data$Sgen_sig <- TMB_Inputs$Sgen_sig # set variance to be used for likelihood for estimating Sgen
#   
#   # set-up init params
#   param <- list()
#   param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
#   param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
#   param$logSigma <- rep(-2, N_Stocks)
#   param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
#   param$logMuA <- TMB_Inputs$logA_Start
#   param$logSigmaA <- 1
#   param$gamma <- 0
#   
#   # specify data
#   data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
#   data$logSurv_3 <- log(SRDat$STAS_Age_3)
#   data$logSurv_4 <- log(SRDat$STAS_Age_4)
#   
#   muSurv <- SRDat %>% group_by(CU_ID) %>% 
#     summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
#   data$muLSurv <- log(muSurv$muSurv)
#   
#   
#   data$Tau_dist <- TMB_Inputs$Tau_dist
#   data$gamma_mean <- TMB_Inputs$gamma_mean 
#   data$gamma_sig <- TMB_Inputs$gamma_sig
#   data$logMuA_mean <- TMB_Inputs$logMuA_mean 
#   data$logMuA_sig <- TMB_Inputs$logMuA_sig
#   data$Tau_A_dist <- TMB_Inputs$Tau_A_dist 
#   
#   
#   # Phase 1 estimate SR params ============
#   map <- list(logSgen=factor(rep(NA, data$N_Stks))) # Fix Sgen
#   obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)
#   
#   opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
#   pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
#   
#   
#   # -- pull out SMSY values
#   All_Ests <- data.frame(summary(sdreport(obj)))
#   All_Ests$Param <- row.names(All_Ests)
#   SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
#   
#   pl$logSgen <- log(0.3*SMSYs)
#   
#   
#   # Phase 2 get Sgen, SMSY etc. =================
#   #map2 <- list(B_0=factor(NA), B_1=factor(NA))
#   obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA")
#   
#   
#   # Create upper bounds vector that is same length and order as start vector that will be given to nlminb
#   upper<-unlist(obj$par)
#   upper[1:length(upper)]<-Inf
#   upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
#   upper<-unname(upper) 
#   
#   lower<-unlist(obj$par)
#   lower[1:length(lower)]<--Inf
#   lower[names(lower) =="logSgen"] <- log(0.001)
#   lower<-unname(lower)
#   
#   opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
#                 upper = upper, lower=lower )
#   
#   pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
#   
#   
#   # Fit mcmc with STAN ==================================================================
#   
#   #fitmcmc_onlySR <- tmbstan(obj, chains=3, iter=nIter, init=unlist(pl2), 
#   #                   control = list(adapt_delta = 0.95),upper=upper, lower=lower)
#   
#   fitmcmc_onlySR <- tmbstan(obj, chains=3, iter=nIter, init=opt$par, 
#                             control = list(adapt_delta = 0.95),upper=upper, lower=lower)
#   
#   ## Can get ESS and Rhat from rstan::monitor
#   mon <- monitor(fitmcmc_onlySR)
#   max(mon$Rhat)
#   min(mon$Tail_ESS)
#   
#   write.csv(mon,paste(priorDir,"/mcmcDiag_Hier_Ricker_Surv_onlySR.csv", sep=""))
#   
#   ## Pairs plot of the fixed effects
#   
#   parNames<-names(fitmcmc_onlySR)[-length(names(fitmcmc_onlySR))]
#   
#   
#   # Marginal posteriors for Sgen
#   png(filename = paste(priorDir,"/postHist_Sgens_Hier_Ricker_Surv_onlySR.png",sep=""))
#   par(mfrow=c(3,2))  
#   
#   # pull out posterior vals
#   post<-as.matrix(fitmcmc_onlySR)
#   
#   hist(exp(post[,'logSgen[1]'])*Scale)
#   hist(exp(post[,'logSgen[2]'])*Scale)
#   hist(exp(post[,'logSgen[3]'])*Scale)
#   hist(exp(post[,'logSgen[4]'])*Scale)
#   hist(exp(post[,'logSgen[5]'])*Scale)
#   
#   dev.off() 
#   
#   
#   
#   ## Trace plot
#   png(filename = paste(priorDir,"/tracePlot_Ricker_Hier_Surv_onlySR.png",sep=""))
#   traceplot(fitmcmc_onlySR, pars=parNames, inc_warmup=TRUE)
#   dev.off()
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
  
  
  # *************************************************************************************
  # Fit Hier_Ricker_Surv_UniformTau model
  # ********************************************************************************

  Mod <- "Aggregate_LRPs_Hier_Surv_UniformTau"

  TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, logMuA_sig = 2,
                             logSigma_lower =log(0.000001), logSigma_upper = log(5),
                             logSigmaA_lower =log(0.000001), logSigmaA_upper = log(5),
                             gamma_mean = 0, gamma_sig = 100, S_dep = 1000, Sgen_sig = 0.5)


  # will define this here, since used so much
  Scale <- TMB_Inputs$Scale

  data <- list()
  data$S <- SRDat$Spawners/Scale
  data$logR <- log(SRDat$Recruits/Scale)
  data$stk <- as.numeric(SRDat$CU_ID)
  N_Stocks <- length(unique(SRDat$CU_Name))
  data$N_Stks <- N_Stocks
  data$yr <- SRDat$yr_num
  data$Bern_Logistic <- as.numeric(Bern_Logistic)



  # Give model year for which will fit logistic model
  # which years do we have obs for all CUs?
  Num_CUs_Over_Time <- EscDat %>%  filter(is.na(Escp) == F) %>% group_by(yr) %>% summarise(n=length(unique((CU_Name))))
  Mod_Yrs <- Num_CUs_Over_Time$yr[Num_CUs_Over_Time$n == max(Num_CUs_Over_Time$n)]
  Logistic_Dat <- EscDat %>% filter(yr %in% Mod_Yrs)

  Agg_Abund <- Logistic_Dat %>% group_by(yr) %>% summarise(Agg_Esc = sum(Escp)) %>%
    mutate(Gen_Mean = rollapply(Agg_Esc, genYrs, gm_mean, fill = NA, align="right"))

  # need year as index
  Logistic_Dat$yr_num <- group_indices(Logistic_Dat, yr) - 1

  if(useGenMean == T){
    GenMean_DF <- EscDat %>% group_by(CU_ID) %>% mutate(Gen_Mean = rollapply(Escp, genYrs, gm_mean, fill = NA, align="right")) %>%
      filter(is.na(Gen_Mean) == F) %>% ungroup() %>% mutate(yr_num = yr - min(yr))

    data$LM_S <- GenMean_DF$Gen_Mean / Scale
    data$LM_yr <- GenMean_DF$yr_num
    data$LM_stk <- GenMean_DF$CU_ID
    Agg_Abund <- Agg_Abund %>% filter(is.na(Gen_Mean) == F)
    data$LM_Agg_Abund <- Agg_Abund$Gen_Mean / Scale
  } else {
    data$LM_S <- Logistic_Dat$Escp / Scale
    data$LM_Agg_Abund <- Agg_Abund$Agg_Esc / Scale
    data$LM_yr <- Logistic_Dat$yr_num
    data$LM_stk <- Logistic_Dat$CU_ID
  }
  # range of agg abund to predict from
  data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
  #range of spawner abundance to predict recruitment from
  data$Pred_Spwn <- rep(seq(0,max(data$S)*1.1,length=100), N_Stocks)
  for (i in 1:N_Stocks) {
    if (i == 1) {
      data$stk_predS <- rep(unique(SRDat$CU_ID)[1],100)
    } else {
      data$stk_predS<-c(data$stk_predS,rep(unique(SRDat$CU_ID)[i],100))
    }
  }
  
  
  data$p <- p

  # set variance to be used for likelihood for estimating Sgen
  data$Sgen_sig <- TMB_Inputs$Sgen_sig

  param <- list()
  param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
  param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
  param$logSigma <- rep(0.5,N_Stocks) #rep(-2, N_Stocks)
  param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale)
  param$B_0 <- 2
  param$B_1 <- 0.1


  data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
  data$logSurv_3 <- log(SRDat$STAS_Age_3)
  data$logSurv_4 <- log(SRDat$STAS_Age_4)
  param$gamma <- 0
  muSurv <- SRDat %>% group_by(CU_ID) %>%
    summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
  data$muLSurv <- log(muSurv$muSurv)
  data$gamma_mean <- TMB_Inputs$gamma_mean
  data$gamma_sig <- TMB_Inputs$gamma_sig
  data$logMuA_mean <- TMB_Inputs$logMuA_mean
  data$logMuA_sig <- TMB_Inputs$logMuA_sig
  param$logMuA <- TMB_Inputs$logA_Start
  param$logSigmaA <- 1

  # Phase 1 estimate SR params
  map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)

  # add upper and lower on logSigma and logSigmaA as a way to create uniform prior
  upper<-unlist(obj$par)
  upper[1:length(upper)]<-Inf
  upper[names(upper) == "logSigma"]<-TMB_Inputs$logSigma_upper
  upper[names(upper) == "logSigmaA"]<-TMB_Inputs$logSigmaA_upper
  upper<-unname(upper)

  lower<-obj$par
  lower[1:length(lower)]<--Inf
  lower[names(lower) == "logSigmaA"]<-TMB_Inputs$logSigmaA_lower
  lower[names(lower) == "logSigma"]<-TMB_Inputs$logSigma_lower
  lower<-unname(lower)

  # now call nlminb with bounds logSigma and logSigmaA
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5), lower=lower, upper=upper)

  pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1


  # pull out SMSY values from phase 1
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

  pl$logSgen <- log(0.3*SMSYs)

  #Phase 2 get Sgen, SMSY etc.
  map2 <- list(B_0=factor(NA), B_1=factor(NA))

  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA", map=map2)

  upper<-unlist(obj$par)
  upper[1:length(upper)]<-Inf
  upper[names(upper) == "logSigma"]<-TMB_Inputs$logSigma_upper
  upper[names(upper) == "logSigmaA"]<-TMB_Inputs$logSigmaA_upper
  upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy
  upper<-unname(upper)

  lower<-unlist(obj$par)
  lower[1:length(lower)]<--Inf
  lower[names(lower) == "logSigma"]<-TMB_Inputs$logSigma_lower
  lower[names(lower) == "logSigmaA"]<-TMB_Inputs$logSigmaA_lower
  lower[names(lower) =="logSgen"] <- log(0.001)
  lower<-unname(lower)


  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
                upper = upper, lower=lower)

  pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2

  # pull out SMSY values from phase 2
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]


  # # Phase 3 fit logistic model
  # # Hold other estimates constant
  # map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
  # # logMuA and logSigmaA disappear after mask, so need to save here
  # HyperParams <- data.frame(summary(sdreport(obj))) %>%
  #   mutate(Param = row.names(summary(sdreport(obj)))) %>%
  #   filter(Param %in% c("logMuA", "logSigmaA"))
  # obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE, map=map3)
  #
  #
  # upper<-unlist(obj$par)
  # upper[1:length(upper)]<-Inf
  # upper[names(upper) == "logSigma"]<-TMB_Inputs$logSigma_upper
  # #upper[names(upper) == "logSigmaA"]<-TMB_Inputs$logSigmaA_upper
  # upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy
  # upper<-unname(upper)
  #
  # lower<-unlist(obj$par)
  # lower[1:length(lower)]<--Inf
  # lower[names(lower) == "logSigma"]<-TMB_Inputs$logSigma_lower
  # #lower[names(lower) == "logSigmaA"]<-TMB_Inputs$logSigmaA_lower
  # lower[names(lower) =="logSgen"] <- log(0.001)
  # lower<-unname(lower)
  #
  # opt <- tryCatch(
  #   {nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
  #           upper = upper, lower=lower)},
  #   error=function(cond) {
  #     # Choose a return value in case of error
  #     return(NA)
  #   },
  #   warning=function(cond) {
  #     return(NA)
  #   }
  # )
  #
  # pl3 <- obj$env$parList(opt$par) # Parameter estimate after phase 3
  #
  # # pull out SMSY values from phase 2
  # All_Ests <- data.frame(summary(sdreport(obj)))
  # All_Ests$Param <- row.names(All_Ests)
  # SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

  # Phase 4 do it all ....

  obj <- MakeADFun(data, pl2, DLL=Mod, random = "logA", silent=TRUE)

  upper<-unlist(obj$par)
  upper[1:length(upper)]<-Inf
  upper[names(upper) == "logSigma"]<-TMB_Inputs$logSigma_upper
  upper[names(upper) == "logSigmaA"]<-TMB_Inputs$logSigmaA_upper
  upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy
  upper<-unname(upper)

  lower<-unlist(obj$par)
  lower[1:length(lower)]<--Inf
  lower[names(lower) == "logSigma"]<-TMB_Inputs$logSigma_lower
  lower[names(lower) == "logSigmaA"]<-TMB_Inputs$logSigmaA_lower
  lower[names(lower) =="logSgen"] <- log(0.001)
  lower<-unname(lower)

  opt <- tryCatch(
    {nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
            upper = upper, lower=lower)},
    error=function(cond) {
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      return(NA)
    }
  )

  library(tmbstan)
  library(rstan)

  fitmcmc_Unif <- tmbstan(obj, chains=3, iter=nIter, init=opt$par,
                     control = list(adapt_delta = 0.95),upper=upper, lower=lower)

    ## Can get ESS and Rhat from rstan::monitor
    mon <- monitor(fitmcmc_Unif)
    max(mon$Rhat)
    min(mon$Tail_ESS)

    write.csv(mon,paste(priorDir,"/mcmcDiag_Hier_Ricker_Surv_Unif.csv", sep=""))

    ## Pairs plot of the fixed effects

    parNames<-names(fitmcmc_Unif)[-length(names(fitmcmc_Unif))]


    # Marginal posteriors for Sgen
    png(filename = paste(priorDir,"/postHist_Sgens_Hier_Ricker_Surv_Unif.png",sep=""))
    par(mfrow=c(3,2))

    # pull out posterior vals
    post<-as.matrix(fitmcmc_Unif)

    hist(exp(post[,'logSgen[1]'])*Scale)
    hist(exp(post[,'logSgen[2]'])*Scale)
    hist(exp(post[,'logSgen[3]'])*Scale)
    hist(exp(post[,'logSgen[4]'])*Scale)
    hist(exp(post[,'logSgen[5]'])*Scale)

    dev.off()



    ## Trace plot
    png(filename = paste(priorDir,"/tracePlot_Ricker_Hier_Surv_Unif.png",sep=""))
    traceplot(fitmcmc_Unif, pars=parNames, inc_warmup=TRUE)
    dev.off()
 