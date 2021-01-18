

run_ScenarioProj <- function(SRDat, EscDat, BMmodel, scenarioName, useGenMean, genYrs, TMB_Inputs, outDir,
                        runMCMC, nMCMC, nProj, ERScalar=NULL, sigER, recCorScalar) {
  
  scenInputDir <- paste(outDir, "SamSimInputs", scenarioName, sep="/")
  scenOutputDir <- paste(outDir, "SamSimOutputs", scenarioName, sep="/")
    
    if (file.exists(scenInputDir) == FALSE){
      dir.create(scenInputDir)
    }
    if (file.exists(scenOutputDir) == FALSE){
      dir.create(scenOutputDir)
    }
  
  # Run MPD fit to parameterize samSim projections ==========
  mpdOut<-get_MPD_Fit(SRDat, BMmodel, TMB_Inputs, outDir)
  # save correlation matrix
  corMatrix<-mpdOut$corMatrix
  corMatrix<-corMatrix * recCorScalar
  corMatrix[col(corMatrix)==row(corMatrix)] <- 1
  write.table(corMatrix, paste(scenInputDir,"cohoCorrMat.csv",sep="/"),row.names=F, col.names=F, sep=",")
  
  # Run MCMC fit to parameterize samSim projections =======
  if (runMCMC == T) {
    mcmcOut<-get_MCMC_Fit(scenarioName, obj= mpdOut$obj, init=mpdOut$par, 
                          upper=mpdOut$upper, lower=mpdOut$lower, nMCMC=nMCMC, Scale = TMB_Inputs$Scale)
    mcmcOut<-as.data.frame(mcmcOut)
    write.csv(mcmcOut, paste(scenInputDir,"cohoRickerSurv_mcmc.csv", sep="/"), row.names=F)
  }
  if (runMCMC == F) {
    mpdOut<-read.csv(paste(scenInputDir,"cohoRickerSurv_mcmc.csv", sep="/"))
    
  }
 
  # Create recruitment input file in required samSim format ========================== 
  cohoRecDatTrim<-data.frame(stk=(SRDat$CU_ID + 1), 
                             yr=SRDat$BroodYear,
                             ets=SRDat$Spawners,
                             totalSpwn=SRDat$Spawners,
                             rec2=rep(0,length(SRDat$BroodYear)),
                             rec3=SRDat$Age_3_Recruits,
                             rec4=SRDat$Age_4_Recruits,
                             rec5=rep(0,length(SRDat$BroodYear)),
                             rec6=rep(0,length(SRDat$BroodYear)))
  write.csv(cohoRecDatTrim, paste(scenInputDir,"cohoRecDatTrim.csv", sep="/"), row.names=F)
 
  
  # Read-in CU pars file and re-write with updated scenario pars =====================
  CUpars<-read.csv(paste(outDir, "SamSimInputs/cohoCUPars.csv",sep="/"))
  # -- specify ER scenario
  CUpars$sigER <- rep(sigER,length(unique(SRDat$CU_ID)))
  # -- fill-in MPD fits
  CUpars$alpha <- mpdOut$All_Ests[grepl("logA", mpdOut$All_Ests$Param), "Estimate" ]
  CUpars$beta0 <- exp(mpdOut$All_Ests[grepl("logB", mpdOut$All_Ests$Param), "Estimate" ])/TMB_Inputs$Scale
  CUpars$sigma <- exp(mpdOut$All_Ests[grepl("logSigma", mpdOut$All_Ests$Param), "Estimate" ])
  CUpars$coef1 <- rep(mpdOut$All_Ests[grepl("gamma", mpdOut$All_Ests$Param), "Estimate" ],length(unique(SRDat$CU_ID)))
  # -- add initial marine survival covariate based on recent average
  means<-SRDat %>% group_by(CU_ID) %>% filter(BroodYear > 1999) %>% summarise(coVariate=mean(STAS_Age_3))
  CUpars$covarInit <- means$coVariate
  # -- specify distribution for marine survival for annual sampling
  means_log<-SRDat %>% filter(BroodYear > 1999) %>% group_by(CU_ID) %>% summarise(coVariate=mean(log(STAS_Age_3)))
  CUpars$mu_logCovar1 <- means_log$coVariate
  sigs_log<-SRDat %>% filter(BroodYear > 1999) %>% group_by(CU_ID) %>% summarise(coVariate=sd(log(STAS_Age_3)))
  CUpars$sig_logCovar1 <- sigs_log$coVariate
  CUpars$min_logCovar <- rep(-5.914504,length(unique(SRDat$CU_ID)))
  CUpars$max_logCovar <- rep(-2.701571,length(unique(SRDat$CU_ID)))
  # -- add mean recruitment, by age
  CUpars$meanRec2 <- rep(0,length(unique(SRDat$CU_ID)))
  rec3<-SRDat %>% group_by(CU_ID) %>% summarize(rec3=Age_3_Recruits/Recruits)
  CUpars$meanRec3 <- mean(rec3$rec3)
  rec4<-SRDat %>% group_by(CU_ID) %>% summarize(rec4=Age_4_Recruits/Recruits)
  CUpars$meanRec4 <- mean(rec4$rec4)
  # -- add median annual recruitment and quantiles
  quantRec<- SRDat %>% group_by(CU_ID) %>% summarize(med=median(Recruits), lowQ=quantile(Recruits,0.25),highQ=quantile(Recruits,0.75))
  CUpars$medianRec <- quantRec$med
  CUpars$lowQRec <- quantRec$lowQ
  CUpars$highQRec <- quantRec$highQ
  
  write.csv(CUpars, paste(scenInputDir,"cohoCUPars.csv", sep="/"), row.names=F)
  
  # Read-in sim par file and re-write with updated scenario pars =====================
  simPars<-read.csv(paste(outDir, "SamSimInputs/cohoSimPar.csv",sep="/"))
  simPars$nameOM<-rep(scenarioName,nrow(simPars))
  simPars$scenario<-paste(simPars$nameOM,simPars$nameMP,sep="_")
  
  write.csv(simPars, paste(scenInputDir,"cohoSimPars.csv", sep="/"), row.names=F)
  
  
  ## Run projections =================================================================================
  
  
  ## Check if necessary packages are available and install if necessary
  listOfPackages <- c("here", "parallel", "doParallel", "foreach",
                      "tidyverse", "tictoc", "samSim")
  newPackages <- listOfPackages[!(listOfPackages %in%
                                    installed.packages()[ , "Package"])]
  if (length(newPackages)) {
    install.packages(newPackages)
  }
  lapply(listOfPackages, require, character.only = TRUE)
  
  dimnames(corMatrix)=NULL
  
  scenNames <- unique(simPars$scenario)
  dirNames <- scenNames#sapply(scenNames, function(x) paste(x, unique(simPars$species),sep = "_"))
  


  
  for (i in 1:nrow(simPars)) {
    
  genericRecoverySim(simPars[i, ], cuPar=CUpars, srDat=cohoRecDatTrim,
          variableCU=FALSE, ricPars=mcmcOut, cuCustomCorrMat = corMatrix,
          dirName=dirNames[[i]], nTrials=nProj, makeSubDirs=FALSE, random=FALSE, outDir=outDir)
  
  }
   
  
  # Read-in projection outputs
  for (i in 1:nrow(simPars)) {
    
    filename<-paste(simPars[i, "nameOM"],simPars[i, "nameMP"],"CUspwnDat.csv",sep="_" )
    filename2<-paste(simPars[i, "nameOM"],simPars[i, "nameMP"],"lrpDat.csv",sep="_" )
    datCUSp.i<-read.csv(here(outDir,"SamSimOutputs", "simData", dirNames[[i]], filename))
    datCUSp.i$expRate<-simPars[i, "canER"] + simPars[i, "usER"] 
    datLRP.i<-read.csv(here(outDir,"SamSimOutputs", "simData", dirNames[[i]], filename2))
    datLRP.i$expRate<-simPars[i, "canER"] + simPars[i, "usER"] 
    
    if (i == 1) {
      projSpwnDat<-datCUSp.i
      projLRPDat<-datLRP.i
    }
    
    if (i > 1) {
      projSpwnDat<-rbind(projSpwnDat,datCUSp.i)
      projLRPDat<-rbind(projLRPDat,datLRP.i)
    }
    
  }
  
  write.csv(projSpwnDat,paste(here(outDir,"SamSimOutputs", "simData"), paste("projSpwnDat_",simPars$nameOM[1],".csv",sep=""),sep="/"))
  write.csv(projLRPDat,paste(here(outDir,"SamSimOutputs", "simData"), paste("projLRPDat_",simPars$nameOM[1],".csv",sep=""),sep="/"))
}



get_MPD_Fit<-function (SRDat, BMmodel, TMB_Inputs, outDir) {
  
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
  
  mpdOut<-list()
  mpdOut$obj<-obj
  mpdOut$par<-opt$par
  mpdOut$upper<-upper
  mpdOut$lower<-lower
  mpdOut$All_Ests<-All_Ests
  mpdOut$corMatrix<-cor_matrix
  
 mpdOut
  
}


get_MCMC_Fit<-function (scenarioName, obj, init, upper, lower, nMCMC, Scale) {
 
# # Fit mcmc with STAN to get parameter estimates for projections ===============================
fitmcmc <- tmbstan(obj, chains=3, iter=nMCMC, init=init,
                   control = list(adapt_delta = 0.99),upper=upper, lower=lower)

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

post_long



}
  