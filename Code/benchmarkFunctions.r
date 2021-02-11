# This file will contain fucntions that calculate percentile benchmarks, independent SR-based benchmarks
#  (i.e., for cases where Sgen estimates are not incorporated into), etc



# GetRickerParams() -  Calculate Ricker params and benchmarks (Smsy, Sgen) using Base R
library(gsl)


#Calculate Ricker params using Base R ============================================================
GetRickerParams <- function(SRDat){
  
  
  CUList<-unique(SRDat$CU_Name)
  
  # Loop over CUs and calculate individual Sgen parameters ========= (Question: should this be done with some sort of apply fxn?)
  for (ss in 1:length(CUList)) {
    
    Dat <- SRDat %>% filter(CU_Name == CUList[ss])
    
    # Calculate CU-level benchmarks ==================
    Dat$RPS <- Dat$Recruits/Dat$Spawners
    Dat$logRPS <-log(Dat$RPS)
    fit <- lm( logRPS ~ Spawners, data=Dat )
    # Get the Ricker alpha and beta parameters and fit statistics (and R2, pValue)
    A <- exp( as.numeric(coef(fit)[1]) )  # Intercept, log(Ricker A)
    B <- -1 * as.numeric( coef(fit)[2] )  # Slope, - Ricker b
  
    # SRep<- log(A) / B
    # # Calculate spawners at MSY (approximation; Hilborn and Walters 1992)
    # SMSY <- SRep * ( 0.5 - 0.07*log(A) )
    
    # Calculate using Scheurell 2016 approach
    SMSY <- (1 - gsl::lambert_W0(exp(1 - log(A)))) / B
  
  # solve for Sgen
    ObjectiveSGen <- function( S, SpMSY, alpha, beta ) {
      # Recruits to get to SMSY (Holt et al. 2009)
      R <- alpha * S * exp( -beta * S )
      # Residual difference between ln(SMSY) and ln(R)
      delta <- log( SpMSY ) - log( R )
      # Calculate the negative log-likelihood
      negLL <- -1 * dnorm( x=delta, mean=0, sd=1, log=TRUE )
      # Return the negative log-likelihood (this is the value to minimize)
      return( negLL )
    }  # End ObjectiveSGen function
  
    opt <- optimize( f=ObjectiveSGen, interval=c(0, SMSY), SpMSY=SMSY,
                   alpha=A, beta=B )
    
    # Get SGen from the optimized output (i.e., minimum neg log-like)
    SGen <- opt$minimum
    # And sigma
    sigma <- summary(fit)[[6]]

    ## Save CU-level benchmark estimates
    if (ss == 1) {
      output.df<-data.frame(CUList[ss], Dat$CU_ID[ss], A, B, SMSY, SGen, sigma)
      names(output.df)<-c("CU_Name", "CU_ID", "est_A", "est_B", "est_Smsy", "est_Sgen", "est_sigma")
    } else {
      newRow.df<-data.frame(CUList[ss], Dat$CU_ID[ss], A, B, SMSY, SGen, sigma)
      names(newRow.df)<-c("CU_Name", "CU_ID","est_A", "est_B", "est_Smsy", "est_Sgen", "est_sigma")
      output.df<-rbind(output.df,newRow.df)
    }
  
  } # End of CU loop
    
  output.df
} # end get Ricker param. function


# Function to get percentile benchmarks from historic time series of spawner abundance (e.g., for Inside South Coast Chum)
get_percentile_benchmarks <- function( ) {
  
}