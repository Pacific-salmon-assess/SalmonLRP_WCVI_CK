library(zoo)

# ===== Functions for fitting logistic regression models to data =========================================
#
# Run_Ricker_LRP()    - Run TMB Ricker and LRP estimation, either Hier Ricker or regular (integratedModel = T)
# Run_LRP()           - Run TMB LRP estimation (no Ricker data or estimation) (integratedModel = F)
# Run_LRP_Logistic_inR() - 
# ================================================================================


# Run TMB RIcker and LRP estimation, either Hier Ricker or regular ========================================================
Run_Ricker_LRP <- function(SRDat, EscDat, BMmodel, Bern_Logistic, 
                           useGenMean, genYrs, p,
                           TMB_Inputs) {
  
  # Specify Mod (i.e., TMB model name) specific to each benchmark model type
  Mod <- BMmodel
  
  Scale <- TMB_Inputs$Scale
  
  # Set-up data list for TMB function calls 
  data <- list()
  data$S <- SRDat$Spawners/Scale 
  data$logR <- log(SRDat$Recruits/Scale)
  data$stk <- as.numeric(SRDat$CU_ID)
  N_Stocks <- length(unique(SRDat$CU_Name))
  data$N_Stks <- N_Stocks
  data$yr <- SRDat$yr_num
  data$Bern_Logistic <- as.numeric(Bern_Logistic)
  data$p <- p    # specify p level to use for LRP calculation from logistic regression
  data$Sgen_sig <- TMB_Inputs$Sgen_sig  # set variance to be used for likelihood for estimating Sgen
  
  Num_CUs_Over_Time <- EscDat %>%  filter(is.na(Escp) == FALSE) %>% group_by(yr) %>% summarise(n=length(unique((CU_Name))))
  Mod_Yrs <- Num_CUs_Over_Time$yr[Num_CUs_Over_Time$n == max(Num_CUs_Over_Time$n)]
  Logistic_Dat <- EscDat %>% filter(yr %in% Mod_Yrs)
  
  Agg_Abund <- Logistic_Dat %>% group_by(yr) %>% summarise(Agg_Esc = sum(Escp)) %>%
               mutate(Gen_Mean = rollapply(Agg_Esc, genYrs, gm_mean, fill = NA, align="right"))
  
  # need year as index
  Logistic_Dat$yr_num <- group_by(Logistic_Dat, yr) %>% group_indices() - 1
  
 if(useGenMean == T){
   GenMean_DF <- EscDat %>% group_by(CU_ID) %>% mutate(Gen_Mean = rollapply(Escp, genYrs, gm_mean, fill = NA, align="right")) %>%
                           filter(is.na(Gen_Mean) == FALSE) %>% ungroup() %>% mutate(yr_num = yr - min(yr))
   
   data$LM_S <- GenMean_DF$Gen_Mean / Scale
   data$LM_yr <- GenMean_DF$yr_num
   data$LM_stk <- GenMean_DF$CU_ID
   Agg_Abund <- Agg_Abund %>% filter(is.na(Gen_Mean) == FALSE) 
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
 data$Pred_Spwn <- rep(seq(0,max(data$S)*1.1,length=100), N_Stocks) # vectors of spawner abundance to use for predicted recruits, one vector of length 100 for each stock
 for (i in 1:N_Stocks) { # make a vector of same length as Pred_Spwn with CU ids
   if (i == 1) {
     data$stk_predS <- rep(unique(SRDat$CU_ID)[1],100) 
   } else {
    data$stk_predS<-c(data$stk_predS,rep(unique(SRDat$CU_ID)[i],100))
   }
 }
 
 # Set-up parameter list for TMB function calls
  param <- list()
  param$logA <- rep(TMB_Inputs$logA_Start, N_Stocks)
  param$logSigma <- rep(-2, N_Stocks)
  param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
  param$B_0 <- 2
  param$B_1 <- 0.1
  data$Tau_dist <- TMB_Inputs$Tau_dist
  
  # add data and parameters specific to Ricker model without a survival covariate
  if(Mod %in% c("SR_IndivRicker_NoSurv")) {
        param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
    }
  
  # add data and parameters specific to Ricker model with survival covariate
  if(Mod %in% c("SR_HierRicker_Surv", "SR_IndivRicker_Surv", "SR_HierRicker_SurvCap","SR_IndivRicker_SurvCap")){
    data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
    data$logSurv_3 <- log(SRDat$STAS_Age_3)
    data$logSurv_4 <- log(SRDat$STAS_Age_4)
    muSurv <- SRDat %>% group_by(CU_ID) %>% # get average smolt to adult survival 
       summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
    data$muLSurv <- log(muSurv$muSurv)
    data$gamma_mean <- TMB_Inputs$gamma_mean 
    data$gamma_sig <- TMB_Inputs$gamma_sig
    param$gamma <- 0
    
     # add data and parameters specific to hierarchical models:
     if(Mod %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")){
       data$logMuA_mean <- TMB_Inputs$logMuA_mean 
       data$logMuA_sig <- TMB_Inputs$logMuA_sig
       data$Tau_A_dist <- TMB_Inputs$Tau_A_dist 
       param$logMuA <- TMB_Inputs$logA_Start
       param$logSigmaA <- 1
     }
     # add data and parameters specific to models with carrying capacity prior:
     if(Mod %in% c("SR_HierRicker_SurvCap", "SR_IndivRicker_SurvCap")) {
       param$cap <- TMB_Inputs$cap_mean
       data$cap_mean<-TMB_Inputs$cap_mean
       data$cap_sig<-TMB_Inputs$cap_sig
     }
     # add parameter specific to model without cap on carrying capacity (model parameterized for B):
     if(Mod %in% c("SR_HierRicker_Surv", "SR_IndivRicker_Surv")) {
       param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
     }
  }

  
  # Phase 1: estimate SR params
  map <- list(logSgen=factor(rep(NA, data$N_Stks)), B_0=factor(NA), B_1=factor(NA)) # Determine which parameters to fix
  if(Mod %in% c("SR_IndivRicker_Surv", "SR_IndivRicker_SurvCap", "SR_IndivRicker_NoSurv")){
    obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, map=map)
  } else if(Mod %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")){
    obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)
  }
  

  # Call optimization:
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e10, iter.max = 1e10)) # LW: increased eval.max and iter.max from 1e5 to 1e10; helped model converge
  
  # Parameter estimate after phase 1 optimization:
  pl <- obj$env$parList(opt$par) # Parameter estimates after phase 1
  
  # pull out SMSY values
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ] # FLAG: Currently this is giving negative SMSY values. Why? 
  # set initial Sgen param as a function of Smsy
  pl$logSgen <- log(0.3*SMSYs) # FLAG: Currently this is returning NaN values, taking log of negative SMSY values
  
  #write.csv(All_Ests, "2020-12-18_estimates_phase1_debugging.csv", row.names=FALSE)

  #Phase 2: get Sgen, SMSY
  map2 <- list(B_0=factor(NA), B_1=factor(NA))
  if(Mod %in% c("SR_IndivRicker_Surv", "SR_IndivRicker_SurvCap", "SR_IndivRicker_NoSurv")){
    obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, map=map2)
  } else if(Mod %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")){
    obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA", map=map2)
  }
  
  ## Create upper & lower bounds vectors that are same length and order as nlminb start vector
  upper<-unlist(obj$par)
  upper[1:length(upper)]<-Inf
  upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy
  upper<-unname(upper)
  
  lower<-unlist(obj$par)
  lower[1:length(lower)]<--Inf
  lower[names(lower) =="logSgen"] <- log(0.001) # constrain Sgen to be positive
  lower<-unname(lower)
  
 # Call optimization:
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5), 
                upper = upper, lower=lower )
  
  # Parameter estimate after phase 2 optimization:
  pl2 <- obj$env$parList(opt$par) 
  
  # Phase 3: fit logistic model; hold other estimates constant
  if(Mod %in% c("SR_IndivRicker_Surv", "SR_IndivRicker_SurvCap", "SR_IndivRicker_NoSurv")){
    obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE)
  } else if(Mod %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")){
    map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
    # logMuA and logSigmaA disappear after mask, so need to save here
    HyperParams <- data.frame(summary(sdreport(obj))) %>% 
      mutate(Param = row.names(summary(sdreport(obj)))) %>% 
      filter(Param %in% c("logMuA", "logSigmaA"))
    obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE, map=map3)
  }
  
  # Create upper & lower bounds vectors that are same length and order as nlminb start vector
  upper<-unlist(obj$par)
  upper[1:length(upper)]<-Inf
  upper[names(upper) =="logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy
  upper<-unname(upper)
  
  lower<-unlist(obj$par)
  lower[1:length(lower)]<--Inf
  lower[names(lower) =="logSgen"] <- log(0.001) # constrain Sgen to be positive
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

    # Create Table of outputs
    All_Ests <- data.frame(summary(sdreport(obj)))
    All_Ests$Param <- row.names(All_Ests)
    if(Mod %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")){
      All_Ests <- full_join(All_Ests, HyperParams)
    }

    # put together readable data frame of values
    All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
    All_Ests$Mod <- Mod
    All_Ests$CU_ID[!(All_Ests$Param %in% c("logMuA", "logSigmaA", "B_0", "B_1", "Agg_LRP", "gamma", "logsigma", "prod"))] <- rep(0:(N_Stocks-1)) 
    All_Ests <- left_join(All_Ests, unique(SRDat[, c("CU_ID", "CU_Name")]))
    # don't want logged param values
    # need to unlog A,B,Sigma -- took A off this list for now
    All_Ests$Estimate[All_Ests$Param %in% c( "logB", "logSigma")] <- exp(All_Ests$Estimate[All_Ests$Param %in% c( "logB", "logSigma")] )
    All_Ests$Param[All_Ests$Param == "logSigma"] <- "sigma"
    All_Ests$Param[All_Ests$Param == "logB"] <- "B"
    #All_Ests$Param[All_Ests$Param == "logA"] <- "A" 
    All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% 
      mutate(Std..Error = Std..Error/Scale) # LW deleted a duplicate of this row. B gets divided by scale to unscale
    All_Ests[All_Ests$Param %in% c("Sgen", "SMSY", "SRep", "cap", "Agg_LRP"), ] <-  All_Ests %>% 
      filter(Param %in% c("Sgen", "SMSY", "SRep","cap","Agg_LRP")) %>% 
           mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
    All_Ests$Param[All_Ests$Param == "SMSY"] <- "Smsy" # LW: moved this below the multiply by scale line, because SMSY wasn't getting scaled
    Preds <- All_Ests %>% filter(Param == "Logit_Preds")
    Preds_Rec <- All_Ests %>% filter(Param == "Rec_Preds")
    All_Ests <- All_Ests %>% filter(!(Param %in% c( "logSgen", "Logit_Preds", "Rec_Preds"))) 
    
    out <- list()
    out$All_Ests <- All_Ests
    
    # also return agg abundance and num CUs over Sgen
    # depending on whether bernoulli or prop, grab correct N
    if(Bern_Logistic == T){
      N_CUs <- obj$report()$All_Above_BM
    } else {
      N_CUs <- obj$report()$N_Above_BM/N_Stocks
    }
    
    Logistic_Data <- data.frame(Mod = Mod, yr = Agg_Abund$yr, 
                                yy = N_CUs, xx = Agg_Abund$Agg_Esc)

    out$Logistic_Data <- Logistic_Data
    
    Logistic_Fits <- data.frame(xx = data$Pred_Abund*Scale, fit = inv_logit(Preds$Estimate),
                                lwr = inv_logit(Preds$Estimate - 1.96*Preds$Std..Error),
                                upr = inv_logit(Preds$Estimate + 1.96*Preds$Std..Error))
    
    out$Preds <- Logistic_Fits
    
    out$LRP <- data.frame(fit = All_Ests %>% filter(Param == "Agg_LRP") %>% pull(Estimate), 
                          lwr = All_Ests %>% filter(Param == "Agg_LRP") %>% mutate(xx =Estimate - 1.96*Std..Error) %>% pull(xx),
                          upr = All_Ests %>% filter(Param == "Agg_LRP") %>% mutate(xx =Estimate + 1.96*Std..Error) %>% pull(xx))
    
    out
    
}



# ==============================================================================================================

# Run TMB estimation for logistic regression fit only ========================================================
Run_LRP <- function(EscDat, Mod, useBern_Logistic, 
                           useGenMean, genYrs, p, TMB_Inputs) {

  # will define this here, since used so much
  Scale <- TMB_Inputs$Scale
  
  # create data vector to give TBM logistic regression
  data <- list()
  param <- list()
  
  N_Stocks <- length(unique(EscDat$CU_Name))
  data$N_Stks <- N_Stocks
  data$Bern_Logistic <- as.numeric(useBern_Logistic)
  
  # Calculate Generation Mean of Aggregate Abundance
  Agg_Abund <- EscDat %>% group_by(yr) %>% summarise(Agg_Esc = sum(Escp)) %>%
    mutate(Gen_Mean = rollapply(Agg_Esc, genYrs, gm_mean, fill = NA, align="right"))
  
  # Give model year for which will fit logistic model
  # -- only include years in which we have obs for all CUs?
  Num_CUs_Over_Time <- EscDat %>% group_by(yr) %>% summarise(n=length(unique((CU_Name))))
  Mod_Yrs <- Num_CUs_Over_Time$yr[Num_CUs_Over_Time$n == max(Num_CUs_Over_Time$n)]
  
  # -- and only include years with geometric mean of aggregagte abundance
  Gen_Mean_Yrs<-Agg_Abund$yr[!is.na(Agg_Abund$Gen_Mean)]
  Mod_Yrs<-Mod_Yrs[Mod_Yrs %in% Gen_Mean_Yrs]
  Logistic_Dat <- EscDat %>% filter(yr %in% Mod_Yrs)
  Agg_Abund <- Agg_Abund %>% filter(yr %in% Mod_Yrs)

   # need year as index
  # Logistic_Dat$yr_num <- group_indices(as.data.frame(Logistic_Dat), yr) - 1
  # LW - replaced line above with below. Reason: ... argument of group_indices() deprecated
  Logistic_Dat$yr_num <- group_by(as.data.frame(Logistic_Dat), yr) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
  
  if (Mod == "ThreshAbund_Subpop1000_ST") data$LM_CU_Status <- Logistic_Dat$HalfGrThresh
  if (Mod == "ThreshAbund_Subpop1000_LT") data$LM_CU_Status <- Logistic_Dat$AllGrThresh
  data$LM_Agg_Abund <- Agg_Abund$Gen_Mean / Scale
  # switch to this if not wanting to fit to geometric mean of aggregate escapement
  #data$LM_Agg_Abund <- Agg_Abund$Agg_Esc / Scale
  # switch to this if wanting to match Arbeider et al. approach of calculating geomean at sub-population level before aggregating (not; will need to add dum as a model input first)
  #dum2 <- dum2 %>% filter(yr %in% Mod_Yrs)
  #data$LM_Agg_Abund <- dum2$AggEscp.gm / Scale
  data$LM_yr <- Logistic_Dat$yr_num
  data$LM_stk <- Logistic_Dat$CU_ID
  # range of agg abund to predict from
  data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
  data$p <- p
 
  param$B_0 <- 2
  param$B_1 <- 0.1

  # Run TMB code
  obj <- MakeADFun(data, param, DLL="ThreshAbund_Subpop1000", silent=TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
  
  # Create Table of outputs
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  
  # Put together readable data frame of values
  All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
  All_Ests$Mod <- Mod
  All_Ests[All_Ests$Param == "Agg_LRP", ] <-  All_Ests %>% filter(Param == "Agg_LRP") %>% 
    mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
  Preds <- All_Ests %>% filter(Param == "Logit_Preds")
  All_Ests <- All_Ests %>% filter(!(Param == "Logit_Preds")) 

  out <- list()
  out$All_Ests <- All_Ests
  
  # also return agg abundance and num CUs over Sgen
  # depending on whether bernoulli or prop, grab correct N
  if(useBern_Logistic == T){
    N_CUs <- obj$report()$All_Above_BM
  } else {
    N_CUs <- obj$report()$N_Above_BM/N_Stocks
  }
  
  Logistic_Data <- data.frame(Mod = Mod, yr = Mod_Yrs, 
                              yy = N_CUs, xx = data$LM_Agg_Abund*Scale)
  
  out$Logistic_Data <- Logistic_Data
  
  Logistic_Fits <- data.frame(xx = data$Pred_Abund*Scale, fit = inv_logit(Preds$Estimate),
                              lwr = inv_logit(Preds$Estimate - 1.96*Preds$Std..Error),
                              upr = inv_logit(Preds$Estimate + 1.96*Preds$Std..Error))
  
  out$Preds <- Logistic_Fits
  
  out$LRP <- data.frame(fit = All_Ests %>% filter(Param == "Agg_LRP") %>% pull(Estimate), 
                        lwr = All_Ests %>% filter(Param == "Agg_LRP") %>% mutate(xx =Estimate - 1.96*Std..Error) %>% pull(xx),
                        upr = All_Ests %>% filter(Param == "Agg_LRP") %>% mutate(xx =Estimate + 1.96*Std..Error) %>% pull(xx))
  
  out
 
}


Run_ProjRicker_LRP<-function(SRDat, EscDat, BMmodel, LRPmodel, useGenMean, genYrs, p, TMB_Inputs, outDir) {
  
  # Run MCMC SR analysis to parameterize projections using available data
  
  Mod <- paste(BMmodel,"noLRP",sep="_")
  
  # Set-up for call to TMB
  Scale <- TMB_Inputs$Scale
  
  data <- list()
  data$Bayes <- 1 # What is this for?
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
  param$logB <- log(1/( (SRDat %>% group_by(CU_ID) %>% summarise(x=quantile(Spawners, 0.8)))$x/Scale) )
  param$logSigma <- rep(-2, N_Stocks)
  param$logSgen <-  log((SRDat %>% group_by(CU_Name) %>%  summarise(x=quantile(Spawners, 0.5)))$x/Scale)
  param$B_0 <- 2
  param$B_1 <- 0.1
  data$Tau_dist <- TMB_Inputs$Tau_dist
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
  #data$logMuA_mean <- TMB_Inputs$logMuA_mean
  #data$logMuA_sig <- TMB_Inputs$logMuA_sig
  #data$Tau_A_dist <- TMB_Inputs$Tau_A_dist

  # range of agg abund to predict from
  #data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund), length.out = 100)
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
  obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)

  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
  pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1


  # -- pull out SMSY values
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]

  pl$logSgen <- log(0.3*SMSYs)


  # Phase 2 get Sgen, SMSY etc. =================
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA")


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

  # Calculate correlation matrix in MPD recruitment residuals ========================
  resids<-as_tibble(data.frame(stock=data$stk,year=data$yr, resid=obj$report()$R_Resid))

  cor_matrix <- resids %>% pivot_wider(names_from = stock, names_prefix="stock",values_from = resid)  %>%
    select(-year) %>% cor()

  projDir <- paste(outDir,"DataOut/ProjectPars", sep="/")
  if (file.exists(projDir) == FALSE){
    dir.create(projDir)
  }

  write.table(cor_matrix, paste(projDir,"cohoCorrMat.csv",sep="/"),row.names=F, col.names=F, sep=",")

  # # Fit mcmc with STAN to get parameter estimates for projections ===============================
  fitmcmc <- tmbstan(obj, chains=3, iter=10000, init=opt$par,
                     control = list(adapt_delta = 0.95),upper=upper, lower=lower)

  ## Can get ESS and Rhat from rstan::monitor
  mon <- monitor(fitmcmc)
  max(mon$Rhat)
  min(mon$Tail_ESS)

  # pull out posterior vals
  post<-as.matrix(fitmcmc)
  post<-as_tibble(post)
  post<-post %>% add_column(iteration=as.numeric(row.names(post)))

  post_long_alpha<-post %>% select(starts_with("logA"), iteration) %>% pivot_longer(starts_with("logA"),names_to="stock", values_to="logA")
  post_long_alpha$stock<-rep(1:5,length=nrow(post_long_alpha))

  post_long_beta<-post %>% select(starts_with("logB"), iteration, gamma) %>% pivot_longer(starts_with("logB"),names_to="stock", values_to="logB")
  post_long_beta$stock<-rep(1:5,length=nrow(post_long_beta))

  post_long_sigma<-post %>% select(starts_with("logSigma") & !starts_with("logSigmaA"), iteration, gamma) %>%
    pivot_longer(starts_with("logSigma"),names_to="stock", values_to="logSigma")
  post_long_sigma$stock<-rep(1:5,length=nrow(post_long_sigma))

  post_long <- post_long_alpha %>%select(stk=stock, alpha=logA) %>% add_column(beta = exp(post_long_beta$logB)/Scale, sigma=exp(post_long_sigma$logSigma), gamma = post_long_beta$gamma)

  write.csv(post_long, paste(projDir,"cohoRickerSurv_mcmc.csv", sep="/"), row.names=F)

   
 #  # Run projections ===============================
 
    
  ## Check if necessary packages are available and install if necessary
  listOfPackages <- c("here", "parallel", "doParallel", "foreach",
                      "tidyverse", "tictoc", "samSim")
  newPackages <- listOfPackages[!(listOfPackages %in%
                                    installed.packages()[ , "Package"])]
  if (length(newPackages)) {
    install.packages(newPackages)
  }
  lapply(listOfPackages, require, character.only = TRUE)
 #  
 #  ## Load relevant input data (will move later to runFraserCoho.r, and input to function)
 #  
   # Simulation run parameters describing different scenarios
   simPar <- read.csv(paste(outDir, "DataIn/ProjPars/cohoSimPar.csv", sep="/"),
                      stringsAsFactors = F)
  # CU-specific parameters
  cuPar <- read.csv(paste(outDir, "DataIn/ProjPars/cohoCUPars_rickerSurv.csv", sep="/"),
                    stringsAsFactors=F)

  # Stock-recruit and catch data that are used to populate the simulation priming
  # period
  srDat <- read.csv(paste(outDir, "DataIn/ProjPars/cohoRecDatTrim.csv", sep="/"),
                    stringsAsFactors=F)

  # Posterior values of  CU-specific stock recruitment parameters for Ricker model; when available, passed and used to calculate alpha, beta and
  # sigma parameters rather than passing point values
  ricPars <- read.csv(paste(projDir,"cohoRickerSurv_mcmc.csv", sep="/"),
                      stringsAsFactors=F)

  #ricPars<-as.data.frame(post_long)

  dimnames(cor_matrix)=NULL
  corMatrix <- cor_matrix

  scenNames <- unique(simPar$scenario)
  dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),
                                                  sep = "_"))

  
# KH: To Add simPar loop ***************************  
# Loop over all scenarios [rows] in simPar ================= 
  genericRecoverySim(simPar[1, ], cuPar=cuPar, srDat=srDat,
                     variableCU=FALSE, ricPars=ricPars, cuCustomCorrMat = corMatrix,
                     dirName="test.co", nTrials=1000, makeSubDirs=FALSE, random=FALSE)
  genericRecoverySim(simPar[2, ], cuPar=cuPar, srDat=srDat,
                     variableCU=FALSE, ricPars=ricPars, cuCustomCorrMat = corMatrix,
                     dirName="test.co", nTrials=1000, makeSubDirs=FALSE, random=FALSE)
  genericRecoverySim(simPar[3, ], cuPar=cuPar, srDat=srDat,
                     variableCU=FALSE, ricPars=ricPars, cuCustomCorrMat = corMatrix,
                     dirName="test.co", nTrials=1000, makeSubDirs=FALSE, random=FALSE)
  genericRecoverySim(simPar[4, ], cuPar=cuPar, srDat=srDat,
                     variableCU=FALSE, ricPars=ricPars, cuCustomCorrMat = corMatrix,
                     dirName="test.co", nTrials=1000, makeSubDirs=FALSE, random=FALSE)
  genericRecoverySim(simPar[5, ], cuPar=cuPar, srDat=srDat,
                     variableCU=FALSE, ricPars=ricPars, cuCustomCorrMat = corMatrix,
                     dirName="test.co", nTrials=1000, makeSubDirs=FALSE, random=FALSE)


  # Read-in projection outputs to create spawner abundance plots
  
  for (i in 1:nrow(simPar)) {
    filename<-paste(simPar[i, "nameOM"],simPar[i, "nameMP"],"CUspwnDat.csv",sep="_" )
    dat.i<-read.csv(here("outputs", "SimData", "test.co", filename))
    dat.i$expRate<-simPar[i, "canER"] + simPar[i, "usER"] 
    if (i == 1) projSpwnDat<-dat.i
    if (i > 1) projSpwnDat<-rbind(projSpwnDat,dat.i)
  }
  
     makeSpawnerPlot<- function(i, plotDat, CUNames) {
        
      plotDat<-projSpwnDat %>% filter(CU==i) %>% group_by(year, expRate) %>% 
        summarise(medSpawn=median(spawners), lwr=quantile(spawners,0.10),upr=quantile(spawners,0.90))
  
      p <- ggplot(data=plotDat, mapping=aes(x=year,y=medSpawn, colour=factor(expRate))) +
       geom_ribbon(data=plotDat, aes(ymin = lwr, ymax = upr, x=year, fill=factor(expRate)), alpha=0.2) +
       geom_line(mapping=aes(x=year, y=medSpawn)) +
       geom_line(data=plotDat %>% filter(year < 18), col="black", size=1) +
       ggtitle(CUNames[i]) +
       xlab("Year") + ylab("Spawners") +
       theme_classic()  
      
     }
     
     ps<-lapply(1:N_Stocks, makeSpawnerPlot, CUNames=unique(SRDat$CU_Name))
     
     pdf(paste(cohoDir,"/Figures/", "IndivSRModel_projections.pdf", sep=""), 
         width=9, height=6)
     do.call(grid.arrange,  ps)
     dev.off()
     

  # Read in projection outputs to create input lists for logistic regression
  
  for (i in 1:nrow(simPar)) {
    filename<-paste(simPar[i, "nameOM"],simPar[i, "nameMP"],"lrpDat.csv",sep="_" )
    dat.i<-read.csv(paste("C:/github/SalmonLRP_RetroEval/outputs/simData/test.co", filename, sep="/"))
    dat.i<-dat.i %>% filter(year > max(SRDat$yr_num)+4)
    dat.i$expRate<-simPar[i, "canER"] + simPar[i, "usER"] 
    if (i == 1) projDat<-dat.i
    if (i > 1) projDat<-rbind(projDat,dat.i)
  }
 
     # Test alternative method of calculating CIs based on distribution of sAg
     projDat.pp<-projDat %>% group_by(ppnCUsLowerBM) %>% 
       summarise(LRP.50=median(sAg), LRP.95=quantile(sAg,0.95),LRP.05=quantile(sAg,0.05))
    # write.csv(projDat.pp, paste(cohoDir,"/Figures/", "Distbased_projLRP.csv", sep=""))
     
  out <- list()
  out$Proj<- data.frame(AggSpawners=projDat[,"sAg"], ppnCUs=projDat[, "ppnCUsLowerBM"])
  out$LRP <- as.data.frame(projDat.pp %>% filter(ppnCUsLowerBM==p) %>% select(fit=LRP.50, lwr=LRP.05, upr=LRP.95))
  
  out
  
}


