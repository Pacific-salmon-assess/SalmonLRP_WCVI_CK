  # Testing out MA's basic model from IFC RPA
# Brooke Davis
# Feb 27/2020

library(R2jags)
library(TMB)
library(dplyr)

# read in data
codeDir<-getwd()
setwd('..')
tmpDir<-getwd()
cohoDir<-paste(tmpDir,"/IFCohoStudy",sep="")
setwd(cohoDir)
CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")
setwd(codeDir)

# first fit exact code using JAGS

Niter=75000
Nthin=50
burnin=25000
  
  
# prepare data to feed into model
# number of sites
sites <- unique(Data$CU_Name)
Ncu <- length(sites)
pAge3 <- Data$Age_3_Recruits/Data$Recruits
Sp <- Data$Spawners
LRS <- log(Data$Recruits/Data$Spawners)

# need muLSurv to get "prod" going to assume this is average weighted?
muLSurv = Data %>% group_by(CU_ID) %>% 
            summarise(muLSurv = mean(STAS_Age_3*(Age_3_Recruits/Recruits) + STAS_Age_4*(Age_4_Recruits/Recruits)))


# create data input for JAGS

Dlist <- list ("Ncu" = Ncu, "Nrecs" = dim(Data)[1], "CUid" = Data$CU_ID+1, "pAge3" = pAge3, 
                "Sp" = Data$Spawners, "LSurvAge3" = log(Data$STAS_Age_3), "LSurvAge4" = log(Data$STAS_Age_4),
               "LRS" = LRS, "muLSurv" = log(muLSurv$muLSurv ))
jags.inits <- function(){
  list("alpha" = rep(2, Ncu), "beta" = rep(0.001, Ncu), "mu_alpha" = 3, "tau_alpha" = 1, "gamma" = 0.5)
}
simoutJAGS <- jags(Dlist, model.file = "JAGS_Code/MA_Ricker.txt",
                   parameters.to.save=c("alpha", "beta", "mu_alpha", "tau_alpha", "gamma", "prod"), n.chains = 3, 
                   n.iter = Niter, n.thin=Nthin, n.burnin=burnin, inits = jags.inits)
  
  #save output for later
 fname <- "MA_Mod"
    saveRDS(simoutJAGS, file=paste(cohoDir, "/DataOut/", fname,".rds", sep=""))
  
  
  # change to bugs output to easily extract values
  simoutBUGS <- simoutJAGS$BUGSoutput
  
  #extract posterior medians and concatenate to vectors
  # all vectors will be length 7
  A <- simoutBUGS$median$alpha
  B <- simoutBUGS$median$beta
  Prod <- simoutBUGS$median$prod
  
  # Calc BM's and CI's around them
  SMSYUp <- SGen <- LowMin <- LowMax <- UppMin <- UppMax <- NULL
  
  Posts<-list()
  for(i in 1:Ncu){
    BMdat <- data.frame(A = simoutBUGS$sims.list$prod[,i], B=simoutBUGS$sims.list$beta[,i])
    Posts[[i]] <- CalcBenchmarks(data=BMdat, threshSMSY=1, giveName=T, nameA="A", nameB="B")
    #return median, quantiles of upper and lower
    SMSYUp[i] <- median(Posts[[i]]$SMSYUp, na.rm=T)
    SGen[i] <- median(Posts[[i]]$SGen, na.rm=T)
    LowCI <- quantile(Posts[[i]]$SGen, c(0.025, 0.975), na.rm=T)
    LowMin[i] <- LowCI[1]
    LowMax[i] <- LowCI[2]
    UppCI <- quantile(Posts[[i]]$SMSYUp, c(0.025, 0.975), na.rm=T)
    UppMin[i] <- UppCI[1]
    UppMax[i] <- UppCI[2]
  }
  
  #only return final year of data and parameter estimates
  AllDat <- data.frame(CU=sites, RickerA=A, RickerB= B, 
                       Mod = rep(fname, Ncu), 
                       "SMSYUp"=SMSYUp, "SGen"=SGen, "LowMin"=LowMin, "LowMax"=LowMax, "UppMin"=UppMin, 
                       "UppMax"=UppMax)
  
  # Now compare these values to 
  

CalcBenchmarks <- function( data , giveName=F, nameA, nameB, threshSMSY) {
  if(giveName==F){
    A <- data$RickerA
    B <- data$RickerB
  } else{
    A <- data[, which(names(data)==nameA)]
    B <- data[, which(names(data)==nameB)]
  }
  # Initialize an empty vectors
  SRep <-   SGen <- negLL <- SGen <- SMSY <-  NULL
  # Loop over rows in data
  for( i in 1:nrow(data) ) {
    if( !is.na(A[i]) & !is.na(B[i]) & A[i]>0){
      # Calculate spawners at replacement/equilibrium (Ricker 1975)
      SRep[i] <- log(A[i]) / B[i]
      # Calculate spawners at MSY (approximation; Hilborn and Walters 1992)
      SMSY[i] <- SRep[i] * ( 0.5 - 0.07*log(A[i]) )  
    } else SRep[i] <- SMSY[i] <- SGen[i] <- NA
    # Only proceed if A and B are not NA
    if( !is.na(A[i]) & !is.na(B[i]) & SMSY[i] > 0 & !is.na(SMSY[i])) {
      # Calculate SGen by minimizing the objective function, ObjectiveSGen
      opt <- optimize( f=ObjectiveSGen1, interval=c(0, SMSY[i]), SpMSY=SMSY[i], 
                       alpha=A[i], beta=B[i] )
      # Get SGen from the optimized output (i.e., minimum neg log-like)
      SGen[i] <- opt$minimum
      # Get the minimum negative log-likelihood
      negLL[i] <- opt$objective
    } else { # End if not NA(s), otherwise
      # Set SGen, SMSYUp, negLL to NA
      SGen[i] <- SMSY[i] <- negLL[i] <- NA
    }# End if NA(s)
  }  # End i loop over rows
  # Add benchmarks to data frame
  # delete return( list(SGen=SGen, SMSYUp=SMSY*threshSMSY, SRep=SRep, negLL=negLL) )
  data<-cbind( data, SGen=SGen, SMSYUp=SMSY*threshSMSY, SRep=SRep, negLL=negLL )
  data
}  # End CalcBenchmarks function

ObjectiveSGen1 <- function( S, SpMSY, alpha, beta ) {
  # Recruits to get to SMSY (Holt et al. 2009)
  logR <- log(alpha) + log(S) - beta * S 
  # Residual difference between ln(SMSY) and ln(R)
  delta <- log( SpMSY ) - logR  
  # Calculate the negative log-likelihood
  negLL <- -1 * dnorm( x=delta, mean=0, sd=1, log=TRUE )
  # Return the negative log-likelihood (this is the value to minimize)
  return( negLL )
}  # End ObjectiveSGen function

#======================================================================


