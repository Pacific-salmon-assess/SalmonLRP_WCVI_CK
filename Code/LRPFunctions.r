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
  # Logistic_Dat$yr_num <- group_indices(Logistic_Dat, yr) - 1
  # LW - replaced line above with below. Reason: ... argument of group_indices() deprecated
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
  
  # add data and parameters specific to Ricker model with survival co-variate
  if(Mod %in% c("SR_HierRicker_Surv", "SR_IndivRicker_Surv", "SR_HierRicker_SurvCap","SR_IndivRicker_SurvCap")){
    data$P_3 <- SRDat$Age_3_Recruits/SRDat$Recruits
    data$logSurv_3 <- log(SRDat$STAS_Age_3)
    data$logSurv_4 <- log(SRDat$STAS_Age_4)
    muSurv <- SRDat %>% group_by(CU_ID) %>% # get average smolt to adult survival 
       summarise(muSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))
    data$muLSurv <- log(muSurv$muSurv)
    data$Tau_dist <- TMB_Inputs$Tau_dist
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
  if(Mod %in% c("SR_IndivRicker_Surv", "SR_IndivRicker_SurvCap")){
    obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, map=map)
  } else if(Mod %in% c("SR_HierRicker_Surv", "SR_HierRicker_SurvCap")){
    obj <- MakeADFun(data, param, DLL=Mod, silent=TRUE, random = "logA", map=map)
  }

  # Call optimization:
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
  
  # Parameter estimate after phase 1 optimization:
  pl <- obj$env$parList(opt$par) # Parameter estimates after phase 1
  
  
  # pull out SMSY values
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
  # set initial Sgen param as a function of Smsy
  pl$logSgen <- log(0.3*SMSYs)
  
  #Phase 2: get Sgen, SMSY
  map2 <- list(B_0=factor(NA), B_1=factor(NA))
  if(Mod %in% c("SR_IndivRicker_Surv", "SR_IndivRicker_SurvCap")){
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
  if(Mod %in% c("SR_IndivRicker_Surv", "SR_IndivRicker_SurvCap")){
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
    #All_Ests$Param[All_Ests$Param == "logA"] <- "A"
    All_Ests$Param[All_Ests$Param == "logB"] <- "B"
    All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
    All_Ests$Param[All_Ests$Param == "logSigma"] <- "sigma"
    All_Ests$Param[All_Ests$Param == "SMSY"] <- "Smsy"
    All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
    All_Ests[All_Ests$Param %in% c("Sgen", "SMSY", "SRep", "cap", "Agg_LRP"), ] <-  All_Ests %>% filter(Param %in% c("Sgen", "SMSY", "SRep","cap","Agg_LRP")) %>% 
           mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
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








#============================================================================================

# Fit a logistic model to a binomial response variable that indicates whether all CUs were above their CU-level lower benchmarks
Run_LRP_Logistic_inR <- function(Dat, LBenchmark_byCU, genYrs, useBern_Logistic, useGenMean, p = 0.95) {
  
  Dat<-left_join(Dat, LBenchmark_byCU, by="CU_Name")
  
  # Add generational running average
  # this isn't grouped by stock though!
  #Dat$Gen_Mean <- rollapply(Dat$Escp, genYrs, gm_mean, fill=NA, align="right")
  
  # Specify CU-level status metric to be compared to Sgen
  if(useGenMean == T){
    Dat <- Dat %>% group_by(CU) %>% mutate(Gen_Mean = rollapply(Escp, genYrs, gm_mean, fill = NA, align="right"))  %>%
               filter(is.na(Dat$Gen_Mean) == FALSE)
    Status <-  "Gen_Mean" # compare generational average escapement calculated as a geometric mean
  } else {
    Status <- "Escp" # compare annual escapement
  }
  
  # Calculate proportion above LBM using generational averages (note: CUs with no LBM will be deleted from proportion calc)
  ModDat <- Dat %>% group_by(yr) %>% filter(is.na(LBM)==FALSE) %>% summarise(CUs_Above_LBM = sum(!!as.name(Status) > LBM),
                                                                n_CUs = n(), Prop_Above_LBM = sum(!!as.name(Status) > LBM)/n(),
                                                                Agg_Escp= sum(Escp))
  
  # And add generational geometric mean of aggregate escaoement
  ModDat$Gen_Agg_Escp <- rollapply(ModDat$Agg_Escp, genYrs, gm_mean, fill=NA, align="right")

  # Calculate proportion above LBM
  if(useBern_Logistic == T){   
     ModDat$yy <- ModDat$Prop_Above_LBM
     ModDat$yy[ModDat$Prop_Above_LBM < 1] <- 0
  } else {
     ModDat$yy <- ModDat$Prop_Above_LBM
  }
  
  if(useGenMean == T){
    ModDat$xx <- ModDat$Gen_Agg_Escp
  } else {
    ModDat$xx <- ModDat$Agg_Escp
  }
  # Fit model
  ModDat <-ModDat %>% drop_na()
  if(useBern_Logistic == T){
       Fit_Mod = glm(yy ~ xx, family = binomial, data = ModDat)
  } else {
       Success_Failure <- as.matrix(data.frame(ModDat$CUs_Above_LBM, ModDat$n_CUs - ModDat$CUs_Above_LBM))
       Fit_Mod = glm( Success_Failure ~ xx, family = binomial, data = ModDat)
  }
  
  #Get LRP (calculated as aggregate Escap number associated with 95% probability of all > LRP)
  # Note: this is done by re-arraranging the standard logisitic equation to solve for LRP when p(x) = p
  LRP <- (log(p/(1-p)) - Fit_Mod$coefficients[[1]])/ Fit_Mod$coefficients[[2]]
  # use MASS function to get "dose" 
  Dose <- dose.p(Fit_Mod, p=p)
      
  #  - Make x vector to predict with this model, for plotting
  xx <- data.frame(xx = seq(0, max(ModDat$xx*1.25), by=(max(ModDat$xx*1.25)/1000)))
  
  # - Create model predictions that include standard error fit
  preds <- predict(Fit_Mod, newdata = xx, type = "link", se.fit = TRUE)
      
  # Create a confidence interval (lwr, upr) on the link scale as the fitted value plus or minus 1.96 times the standard error
  critval <- 1.96 ## approx 95% CI
  upr <- preds$fit + (critval * preds$se.fit)
  lwr <- preds$fit - (critval * preds$se.fit)
  fit <- preds$fit
      
  # Transform confidene interval using the inverse of the link function to map the fitted values and the upper and lower limits of the interval back on to the response scale
  fit2 <- Fit_Mod$family$linkinv(fit)
  upr2 <- Fit_Mod$family$linkinv(upr)
  lwr2 <- Fit_Mod$family$linkinv(lwr)
      
  # Load predicted data on response scale into dataframe
  preddata<-xx
  preddata$fit<-fit2
  preddata$lwr <- lwr2 
  preddata$upr <- upr2 
      
  ### Calculate confidence intervals for LRP 
  LRP_lwr <- LRP - critval*attr(Dose, "SE")[[1]]
  LRP_upr <- LRP + critval*attr(Dose, "SE")[[1]]
  # These seem much tighter than other method?
  

  
    list.out<-list()
    list.out$Logistic_Data <- ModDat
    list.out$model <- Fit_Mod
    list.out$Preds <- preddata
    list.out$LRP<-data.frame(fit = LRP, lwr = LRP_lwr, upr = LRP_upr)


  list.out
}
