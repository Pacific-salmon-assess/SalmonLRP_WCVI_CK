# 


# Run TMB RIcker and LRP estimation, either Hier Ricker or regular ========================================================
Run_Ricker_LRP_UniformTau <- function(SRDat, EscDat, BM_Mod, Bern_Logistic, 
                                      useGenMean, genYrs, p,
                                      TMB_Inputs) {
  
  # Specify Mod (i.e., TMB model name) specific to each benchmark model type
  Mod <- BM_Mod
  
  
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
  summary(sdreport(obj))
  
  # pull out SMSY values
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
                upper = upper, lower=lower )
  
  pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
  summary(sdreport(obj))
  
  # Phase 3 fit logistic model
  # Hold other estimates constant
  map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
  # logMuA and logSigmaA disappear after mask, so need to save here
  HyperParams <- data.frame(summary(sdreport(obj))) %>% 
    mutate(Param = row.names(summary(sdreport(obj)))) %>% 
    filter(Param %in% c("logMuA", "logSigmaA"))
  obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE, map=map3)
  
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
  
  # Create Table of outputs
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  All_Ests <- full_join(All_Ests, HyperParams)
  
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
  All_Ests$Param[All_Ests$Param == "logSigma"] <- "sigma"
  All_Ests$Param[All_Ests$Param == "SMSY"] <- "Smsy"
  All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
  All_Ests[All_Ests$Param %in% c("Sgen", "SMSY", "Agg_LRP"), ] <-  All_Ests %>% filter(Param %in% c("Sgen", "SMSY", "Agg_LRP")) %>% 
    mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
  Preds <- All_Ests %>% filter(Param == "Logit_Preds")
  All_Ests <- All_Ests %>% filter(!(Param %in% c( "logSgen", "Logit_Preds"))) 
  
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


# Run TMB RIcker and LRP estimation, either Hier Ricker or regular ========================================================
RunStan_Ricker_LRP_UniformTau <- function(SRDat, EscDat, BM_Mod, Bern_Logistic, 
                                      useGenMean, genYrs, p,
                                      TMB_Inputs) {
  
  # Specify Mod (i.e., TMB model name) specific to each benchmark model type
  Mod <- BM_Mod
  
  
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
  upper<-obj$par
  upper[1:length(upper)]<-Inf
  upper[names(upper) == "logSigma"]<-TMB_Inputs$logSigma_upper
  upper[names(upper) == "logSigmaA"]<-TMB_Inputs$logSigmaA_upper
  
  lower<-obj$par
  lower[1:length(lower)]<--Inf
  lower[names(lower) == "logSigmaA"]<-TMB_Inputs$logSigmaA_lower
  lower[names(lower) == "logSigma"]<-TMB_Inputs$logSigma_lower
  
  # now call nlminb with bounds logSigma and logSigmaA
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5), lower=unname(lower), upper=unname(upper))
  
  
  
  pl <- obj$env$parList(opt$par) # Parameter estimate after phase 1
  summary(sdreport(obj))
  
  # pull out SMSY values
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  SMSYs <- All_Ests[grepl("SMSY", All_Ests$Param), "Estimate" ]
  
  pl$logSgen <- log(0.3*SMSYs)
  
  #Phase 2 get Sgen, SMSY etc.
  map2 <- list(B_0=factor(NA), B_1=factor(NA))
  
  obj <- MakeADFun(data, pl, DLL=Mod, silent=TRUE, random = "logA", map=map2)
  
  upper<-obj$par
  upper[1:length(upper)]<-Inf
  upper[names(upper) == "logSigma"]<-TMB_Inputs$logSigma_upper
  upper[names(upper) == "logSigmaA"]<-TMB_Inputs$logSigmaA_upper
  
  upper[names(upper) == "logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
  
  lower<-obj$par
  lower[1:length(lower)]<--Inf
  lower[names(lower) == "logSigma"]<-TMB_Inputs$logSigma_lower
  lower[names(lower) == "logSigmaA"]<-TMB_Inputs$logSigmaA_lower
  
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
                upper = unname(upper), lower=unname(lower) )
  
  pl2 <- obj$env$parList(opt$par) # Parameter estimate after phase 2
  summary(sdreport(obj))
  
  # Phase 3 fit logistic model
  # Hold other estimates constant
  map3 = list(logMuA = as.factor(NA), logSigmaA = as.factor(NA))
  # logMuA and logSigmaA disappear after mask, so need to save here
  HyperParams <- data.frame(summary(sdreport(obj))) %>% 
    mutate(Param = row.names(summary(sdreport(obj)))) %>% 
    filter(Param %in% c("logMuA", "logSigmaA"))
  obj <- MakeADFun(data, pl2, DLL=Mod, silent=TRUE, map=map3)
  
  opt <- tryCatch(
    {nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5),
            upper = unname(upper), lower=unname(lower))},
    error=function(cond) {
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      return(NA)
    }
  )   
  
  
  
  # Phase 4 do it all ....
  
  pl3<-obj$env$parList(opt$par)
  
  obj <- MakeADFun(data, pl3, DLL=Mod, silent=TRUE)
  
  upper<-obj$par
  upper[1:length(upper)]<-Inf
  upper[names(upper) == "logSigma"]<-TMB_Inputs$logSigma_upper
  upper[names(upper) == "logSigmaA"]<-TMB_Inputs$logSigmaA_upper
  
  upper[names(upper) == "logSgen"] <- log(SMSYs) # constrain Sgen to be less than Smsy (To do: confirm with Brooke)
  
  lower<-obj$par
  lower[1:length(lower)]<--Inf
  lower[names(lower) == "logSigma"]<-TMB_Inputs$logSigma_lower
  lower[names(lower) == "logSigmaA"]<-TMB_Inputs$logSigmaA_lower

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
  
  
  fitmcmc <- tmbstan(obj, chains=3, iter=100000, init=list(opt$par), 
                     control = list(adapt_delta = 0.95))
  
  ### The following is based on: https://github.com/kaskr/tmbstan#examples
  
  ## Can get ESS and Rhat from rstan::monitor
  mon <- monitor(fitmcmc)
  max(mon$Rhat)
  min(mon$Tail_ESS)
  
  ## Pairs plot of the fixed effects
 # pairs(fitmcmc, pars=names(obj$par))
  ## Trace plot
 # traceplot(fitmcmc, pars=names(obj$par), inc_warmup=TRUE)

  # pull out posterior vals
  post<-as.matrix(fitmcmc)
  
  # pull out logistic fit values
  AllAbove_Fit_Med <- obj$report(post[1,-ncol(All_Ests)])$All_Above_BM 
  NAbove_Fit_Med <- obj$report(post[1,-ncol(All_Ests)])$N_Above_BM
  
  AllAbove_Fit <- matrix(NA, nrow=nrow(post), ncol = length(AllAbove_Fit_Med))
  NAbove_Fit <- matrix(NA, nrow=nrow(post), ncol = length(NAbove_Fit_Med))
  
  for(i in 1:nrow(post)){
    r <- obj$report(post[i,-ncol(post)])
    AllAbove_Fit[i,] <- r$All_Above_BM
    NAbove_Fit[i,] <- r$N_Above_BM
  }
  # ============================================
  
  # Create Table of outputs
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  All_Ests <- full_join(All_Ests, HyperParams)
  
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
  All_Ests$Param[All_Ests$Param == "logSigma"] <- "sigma"
  All_Ests$Param[All_Ests$Param == "SMSY"] <- "Smsy"
  All_Ests[All_Ests$Param == "B",] <- All_Ests %>% filter(Param == "B") %>% mutate(Estimate = Estimate/Scale) %>% mutate(Std..Error = Std..Error/Scale)
  All_Ests[All_Ests$Param %in% c("Sgen", "SMSY", "Agg_LRP"), ] <-  All_Ests %>% filter(Param %in% c("Sgen", "SMSY", "Agg_LRP")) %>% 
    mutate(Estimate = Estimate*Scale) %>% mutate(Std..Error = Std..Error*Scale)
  Preds <- All_Ests %>% filter(Param == "Logit_Preds")
  All_Ests <- All_Ests %>% filter(!(Param %in% c( "logSgen", "Logit_Preds"))) 
  
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


## Fit with unifom priors ======================================================================


#TMB_Inputs_Uniform <- list(Scale = 1000, logA_Start = 1, logMuA_mean = 1, logMuA_sig = 2,
#                           logSigma_lower =log(0.000001), logSigma_upper = log(10),
#                           logSigmaA_lower =log(0.000001), logSigmaA_upper = log(10),
#                           gamma_mean = 0, gamma_sig = 100, S_dep = 1000, Sgen_sig = 1)


#Uniform_Mod <-Run_Ricker_LRP_UniformTau(SRDat = CoSRDat, EscDat = CoEscpDat, BM_Mod = "SR_HierRicker_Surv_UniformTau",
#                                        Bern_Logistic = F, useGenMean = F, genYrs=3, p=0.95, TMB_Inputs = TMB_Inputs_Uniform)


#Uniform_Mod_MCMC <-RunStan_Ricker_LRP_UniformTau(SRDat = CoSRDat, EscDat = CoEscpDat, BM_Mod = "SR_HierRicker_Surv_UniformTau",
#                                        Bern_Logistic = F, useGenMean = F, genYrs=3, p=0.95, TMB_Inputs = TMB_Inputs_Uniform)

#setwd(priorDir)



