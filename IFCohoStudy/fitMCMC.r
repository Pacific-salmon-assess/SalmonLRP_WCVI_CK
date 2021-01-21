# Code used to fit intergrated LRP models using MCMC in stan

# For stantmb example, see: https://github.com/kaskr/tmbstan#examples


nIter<-10000


library(dplyr)
library(TMB)
library(tmbstan)
library(tidyverse)

setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
cohoDir<-paste(rootDir,"/IFCohoStudy",sep="")

# read in data
setwd(cohoDir)
CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")

setwd(codeDir)
compile("TMB_Files/SR_IndivRicker_Surv_noLRP.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker_Surv_noLRP"))


compile("TMB_Files/SR_IndivRicker.cpp")
dyn.load(dynlib("TMB_Files/SR_IndivRicker"))


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

# Start MCMC fit ==========================================

# Option 1) Fit hierarchical model with survival covariate; no LRP estimation
# TMB_Inputs <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, 
#                       logMuA_sig = sqrt(2), Tau_dist = 0.1, Tau_A_dist = 0.1, 
#                       gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 1)
# 
# Mod <- "SR_HierRicker_Surv_noLRP"

# # Option 2) Fit individual models with survival covariate; no LRP estimation
# TMB_Inputs <- list(Scale = 1000, logA_Start = 1,
#      Tau_dist = 0.1,
#      gamma_mean = 0, gamma_sig = 10, S_dep = 1000, Sgen_sig = 0.4)
# 
# Mod <- "SR_IndivRicker_Surv_noLRP"


# Option 3) Fit individual models without survival covariate; no LRP estimation
TMB_Inputs <- list(Scale = 1000, logA_Start = 1,
     Tau_dist = 0.1, Sgen_sig = 0.95)

Mod <- "SR_IndivRicker"


# Set-up for call to TMB
Scale <- TMB_Inputs$Scale

data <- list()
data$Bayes <- 1
data$S <- SRDat$Spawners/Scale
data$logR <- log(SRDat$Recruits/Scale)
data$stk <- as.numeric(SRDat$CU_ID)
N_Stocks <- length(unique(SRDat$CU_Name))
data$N_Stks <- N_Stocks
data$yr <- SRDat$yr_num
data$Sgen_sig <- TMB_Inputs$Sgen_sig # set variance to be used for likelihood for estimating Sgen

# set-up init params
param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
data$Tau_dist <- TMB_Inputs$Tau_dist


param <- list()
param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
param$logSigma <- rep(-2, N_Stocks)
param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale)

if (Mod == "SR_HierRicker_Surv_noLRP" | Mod == "SR_IndivRicker_Surv_noLRP") {
  param$gamma <- 0
}

if (Mod == "SR_HierRicker_Surv_noLRP") {
  param$logMuA <- TMB_Inputs$logA_Start
  param$logSigmaA <- 1
}

# specify data
data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits

if (Mod == "SR_HierRicker_Surv_noLRP" | Mod == "SR_IndivRicker_Surv_noLRP") {
  data$logSurv_3 <- log(SRDat$STAS_Age_3)
  data$logSurv_4 <- log(SRDat$STAS_Age_4)
}

if (Mod == "SR_HierRicker_Surv_noLRP" | Mod == "SR_IndivRicker_Surv_noLRP") {
muSurv <- SRDat %>% group_by(CU_ID) %>%
  summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
  data$muLSurv <- log(muSurv$muSurv)
}

data$Tau_dist <- TMB_Inputs$Tau_dist

if (Mod == "SR_HierRicker_Surv_noLRP") {
  data$logMuA_mean <- TMB_Inputs$logMuA_mean
  data$logMuA_sig <- TMB_Inputs$logMuA_sig
  data$Tau_A_dist <- TMB_Inputs$Tau_A_dist
}

if (Mod == "SR_HierRicker_Surv_noLRP" | Mod == "SR_IndivRicker_Surv_noLRP") {
  data$gamma_mean <- TMB_Inputs$gamma_mean
  data$gamma_sig <- TMB_Inputs$gamma_sig
}


# range of spawner abundance to predict recruitment from
data$Pred_Spwn <- rep(seq(0,max(data$S)*1.1,length=100), N_Stocks) # vectors of spawner abundance to use for predicted recruits, one vector of length 100 for each stock
for (i in 1:N_Stocks) { # make a vector of same length as Pred_Spwn with CU ids
  if (i == 1) {
    data$stk_predS <- rep(unique(SRDat$CU_ID)[1],100) 
  } else {
    data$stk_predS<-c(data$stk_predS,rep(unique(SRDat$CU_ID)[i],100))
  }
}

# Phase 1 estimate SR params ============
map <- list(logSgen=factor(rep(NA, data$N_Stks))) # Fix Sgen



if (Mod == "SR_HierRicker_Surv_noLRP") {
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)
} else {
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE,  map=map)
}

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1


# -- pull out SMSY values
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

pl$logSgen <- log(0.3*SMSYs)


# Phase 2 get Sgen, SMSY etc. =================
if (Mod == "SR_HierRicker_Surv_noLRP") {
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA")
} else {
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE)
}


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
              upper = upper, lower=lower )

pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2

All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)


# Calculate cor matrix from MPD fits
resids<-as_tibble(data.frame(stock=data$stk,year=data$yr, resid=obj$report()$R_Resid))
cor_matrix <- resids %>% pivot_wider(names_from = stock, names_prefix="stock",values_from = resid)  %>%
  select(-year) %>% cor()
write.table(cor_matrix, paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/", Mod, "_MPDCorrMat.csv", sep=""),row.names=F, col.names=F, sep=",")



# Fit mcmc with STAN ==================================================================

fitmcmc <- tmbstan(obj, chains=3, iter=nIter, init=opt$par,
                          control = list(adapt_delta = 0.99),upper=upper, lower=lower)

## Can get ESS and Rhat from rstan::monitor
mon <- monitor(fitmcmc)
max(mon$Rhat)
min(mon$Tail_ESS)

# pull out posterior vals
post<-as.matrix(fitmcmc)
post<-as_tibble(post)
post<-post %>% add_column(iteration=as.numeric(row.names(post)))

if (Mod == "SR_IndivRicker") {
  post_long_alpha<-post %>% select(starts_with("logA"), iteration) %>% pivot_longer(starts_with("logA"),names_to="stock", values_to="logA")
} else {
  post_long_alpha<-post %>% select(starts_with("logA"), iteration, gamma) %>% pivot_longer(starts_with("logA"),names_to="stock", values_to="logA")
}
post_long_alpha$stock<-rep(1:5,length=nrow(post_long_alpha))

if (Mod == "SR_IndivRicker") {
  post_long_beta<-post %>% select(starts_with("logB"), iteration) %>% pivot_longer(starts_with("logB"),names_to="stock", values_to="logB")
} else {
  post_long_beta<-post %>% select(starts_with("logB"), iteration, gamma) %>% pivot_longer(starts_with("logB"),names_to="stock", values_to="logB")
}
post_long_beta$stock<-rep(1:5,length=nrow(post_long_beta))

if (Mod == "SR_IndivRicker") {  
  post_long_sigma<-post %>% select(starts_with("logSigma") & !starts_with("logSigmaA"), iteration) %>% 
    pivot_longer(starts_with("logSigma"),names_to="stock", values_to="logSigma")
} else {
  post_long_sigma<-post %>% select(starts_with("logSigma") & !starts_with("logSigmaA"), iteration, gamma) %>% 
    pivot_longer(starts_with("logSigma"),names_to="stock", values_to="logSigma")
}
post_long_sigma$stock<-rep(1:5,length=nrow(post_long_sigma))

if (Mod == "SR_IndivRicker") { 
  post_long <- post_long_alpha %>%select(stk=stock, alpha=logA) %>% add_column(beta = exp(post_long_beta$logB)/Scale, sigma=exp(post_long_sigma$logSigma))
} else {
  post_long <- post_long_alpha %>%select(stk=stock, alpha=logA) %>% add_column(beta = exp(post_long_beta$logB)/Scale, sigma=exp(post_long_sigma$logSigma), gamma = post_long_beta$gamma)
}
write.csv(post_long, paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",Mod,"_mcmcPosterior.csv", sep=""),row.names=F)


# ========================================================================================================
# Create plots
# ================================================================================================

CU_list<-unique(SRDat[, "CU_Name"])
CUID_list<-unique(SRDat[, "CU_ID"])
nCUs<-length(CU_list)


plotPostHist<-function(x, post, parName, Scale, applyExp=T, CUNames) {

  margPost<-post[,paste(parName,"[",x,"]", sep="")]
  
  if (applyExp == TRUE) {
    margPost <- exp(margPost)
  }
  if (parName == "logSgen") {
    margPost <- margPost*Scale
    xlab<-"Sgen"
  }
  if (parName == "logA") {
    xlab<-"Alpha(logA)"
  }
  if (parName == "logSigma") {
    margPost <- margPost*Scale
    xlab<-"Sigma"
  }
  if (parName == "logB") {
    margPost <- margPost/Scale
    xlab<-"Beta"
  }
  
  
  margPost<-as.data.frame(margPost)
 
  p<-ggplot(margPost, aes(x=margPost)) + geom_histogram() +
    labs(title=CUNames[x], x=xlab, y = "Count")
 
}

library(gridExtra)

ps<-list()
ps<-lapply(1:nCUs, plotPostHist, post=as.matrix(fitmcmc), parName="logSgen", Scale=1000, applyExp=T, CUNames=CU_list)
pdf(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",Mod,"SgenPost.pdf", sep=""))
do.call(grid.arrange, ps)
dev.off()

ps<-list()
ps<-lapply(1:nCUs, plotPostHist, post=as.matrix(fitmcmc), parName="logA", Scale=1000, applyExp=F, CUNames=CU_list)
pdf(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",Mod,"AlphaPost.pdf", sep=""))
do.call(grid.arrange, ps)
dev.off()

ps<-list()
ps<-lapply(1:nCUs, plotPostHist, post=as.matrix(fitmcmc), parName="logSigma", Scale=1000, applyExp=T, CUNames=CU_list)
pdf(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",Mod,"SigmaPost.pdf", sep=""))
do.call(grid.arrange, ps)
dev.off()

ps<-list()
ps<-lapply(1:nCUs, plotPostHist, post=as.matrix(fitmcmc), parName="logB", Scale=1000, applyExp=T, CUNames=CU_list)
pdf(paste("C:/github/SalmonLRP_RetroEval/IFCohoStudy/DataOut/ModelFits/",Mod,"BetaPost.pdf", sep=""))
do.call(grid.arrange, ps)
dev.off()

# ========================================================================================================
# View mcmc fits
# ================================================================================================



library(shinystan)
launch_shinystan(fitmcmc)

class(fitmcmc)
methods(class="stanfit")
## Pairs plot of the fixed effects
pairs(fitmcmc, pars=names(obj$par))
## Trace plot
traceplot(fitmcmc, pars=names(obj$par), inc_warmup=TRUE)

post<-as.matrix(fitmcmc)

sd0 <- rep(NA, len=nrow(post))

for(i in 1:nrow(post)){
  r <- obj$report(post[i,-ncol(post)])
  sd0[i] <- r$sd0
}
hist(sd0)
