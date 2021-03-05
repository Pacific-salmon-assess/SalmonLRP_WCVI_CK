#-------------------------------------------------------------------------------
# Functions to run diagnostics on logistic regression

# Function 1: LRdiagnostic, performs diagnostics using the entire time-series 
# of data

# Function 2: LOO_LRdiagnostic, performs leave-one-out cross-validation by 
# iteratively re-running LRP estimation of subsets of data, leaving one data 
# point out each time, and estimating classification accuracy on removed data point
#-------------------------------------------------------------------------------

# -------------------------------------------------------------------------------
# Assumptions of logistic regression
# 1. The relationship between aggregate abundance and logit(proportions) is
#  linear
# 2. The observed proportions are independent (not autocorrelated)
# 3. There are no influential values (outliers)
# 4. There is no multicollinearity among predictors (NA, only 1 predictor)
# -------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Function 1: LRdiagnostic, performs diagnostics using the entire time-series 
# of data

# Steps:
# 1. Estimate Pearson resiudals and deviance residuals (Assumption 3).
#   Are deviance residuals >2?
# 2a. Plot residuals against fitted values (Assumption 1).
#   Is there a trend in residuals over fitted values?
# 2b. Plot autocorrelation among residuals Assumption 2).
#   Are residuals autocorrelated?
# 3. Evaluate statistical significance of model coefficients.
# 4. Evaluate Pearson Chi-squared statistic (goodness of fit).
#   Is there statistical evidence for lack of fit?
# 5. Evaluate Deviance G-squared statistic (goodness of fit)
#   Is there statistical evidence for lack of fit?
# 6. Evaluate quasi-Rsquared.
#   What is the ratio of the fit to the null model?
# 7. Evaluate Wald test
#   Is deviance significantly reduced with predictor, Aggregate abundance?
# 8. Evaluate hit rate.
#   What is the classification accuracy of the LRP based on the logistic
#   regression?
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Arguments:
# All_Ests = Dataframe containing parameters of logistic regression, B_0 and B_1
#   with p values (All_Ests <- data.frame(summary(sdreport(obj), p.value=TRUE)))
# AggAbund = Vector of scaled observed aggregate abundances, with NAs removed
# obsPpnAboveBM = Vectors of observed ppn of CUs above their lower benchmarks
#   with NAs removed
# p = the proportion of CUs that must be > their benchmark when defining LRP
# nLL = negLogLikehood = ans from TMB file, outputted with REPORT(ans); in TMB
#   and then called in R with obj$report()$ans
# dir = name of directory where plots should be saved followed by /
# plotname = filename for residual plots

# # Carrie's inputs for testing
# zz <- Get.LRP(remove.EnhStocks = TRUE)
# All_Ests <- zz$out$All_Ests
# AggAbundRaw <- zz$SMU_Esc
# digits <- count.dig(AggAbundRaw)
# ScaleSMU <- min(10^(digits -1 ), na.rm=T)
# AggAbund <- AggAbundRaw/ScaleSMU
# AggAbund <- AggAbund[!is.na(AggAbund)]
# obsPpnAboveBM <- zz$out$Logistic_Data$yy
# p <- zz$LRPppn
# nLL <- zz$nLL
# dir <- "DataOut/"
# plotname <- "WCVI_ResidPlots"
# 
# input<- list(All_Ests=All_Ests, AggAbund=AggAbund,
#              obsPpnAboveBM=obsPpnAboveBM, p=p, nLL=nLL,dir="", plotname="test")
# save(input, file="DataIn/Input_LRdiagnostics.rda")
# 
# LRdiagOut <- LRdiagnostics(All_Ests=All_Ests, AggAbund=AggAbund, 
#                             obsPpnAboveBM=obsPpnAboveBM, p=p, nLL=nLL, dir=dir, 
#                             plotname=plotname)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Returns:
# - png of plots of residuals and their autocorrelation (step 2) stored in 
#   directory, dir
# - list of outputs from steps 1, 3-8
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Function:

LRdiagnostics <- function(All_Ests, AggAbund, obsPpnAboveBM, p, nLL, dir, 
                          plotname){
  
  #-------------------------------------------------------------------------------
  # Source functions and libraries
  #source("R/helperFunctions.r")
  library(patchwork)
  #-------------------------------------------------------------------------------
  
  
  #-------------------------------------------------------------------------------
  # Step 1. Estimate Pearson resiudals and deviance residuals (Assumption 3). 
  #   Are deviance residuals >2?
  #-------------------------------------------------------------------------------
  
  # Get observed and predicted ppn of CUs above their lower benchmark
  B_0 <- All_Ests %>% filter(Param=="B_0") %>% pull(Estimate)
  B_1 <- All_Ests %>% filter(Param=="B_1") %>% pull(Estimate)
  predPpnAboveBM <- inv_logit(B_0 + B_1*AggAbund)
  
  
  # Pearson residuals: Eq3.15 https://data.princeton.edu/wws509/notes/c3s8
  # setting n=1 (number of trials at each observation of x)
  
  PearResid <- ( obsPpnAboveBM - predPpnAboveBM ) / sqrt( predPpnAboveBM * 
                                                            (1 - predPpnAboveBM) ) 
  
  # Deviance residual: Eq3.16 https://data.princeton.edu/wws509/notes/c3s8
  # setting n=1 (number of trials at each observation of x)
  # To avoid NANs, use ifelse statement in the eqn, as suggested here:
  # https://www.datascienceblog.net/post/machine-learning/interpreting_
  # generalized_linear_models/#:~:text=For%20type%20%3D%20%22pearson%22%2
  # 0%2C,%CB%86f(xi)
  
  binom.resid <- function(y, mu) {
    y * log( ifelse(y== 0, 1, y/mu)) + (1-y) * log( 
      ifelse(y==1 ,1,(1-y)/(1-mu) ) )  
  }
  
  
  DevResid <- sign(obsPpnAboveBM - predPpnAboveBM ) * 
    sqrt( 2 * binom.resid(y=obsPpnAboveBM, mu=predPpnAboveBM) ) 
  
  # Observations with a deviance residual in excess of two may indicate 
  # lack of fit. (https://data.princeton.edu/wws509/notes/c3s8)
  
  ## Testing. Residuals match output from R using glm objects
  # ModDat <- data.frame(xx=data$LM_Agg_Abund, yy=SMUlogisticData$ppn)
  # Fit_Mod <- glm( yy ~ xx , family = quasibinomial, data=ModDat)
  # B_0 <- Fit_Mod$coef[1]
  # B_1 <- Fit_Mod$coef[2]
  # obsPpnAboveBM <- ModDat$yy
  # predPpnAboveBM <- inv_logit(B_0 + B_1*ModDat$xx)
  # residuals(Fit_Mod, type="pearson")
  # residuals(Fit_Mod, type="deviance")
  
  #-------------------------------------------------------------------------------
  # Step 2:
  # 2a. Plot residuals against fitted values (Assumption 1). 
  #   Is there a trend in residuals over fitted values?
  # 2b. Plot autocorrelation among residuals Assumption 2). 
  #   Are residuals autocorrelated?
  #-------------------------------------------------------------------------------
  
  
  # Put data for diagnostics in a dataframe for plotting
  diagData <- data.frame(predPppnAboveBM = predPpnAboveBM, 
                         PearResid = PearResid, DevResid=DevResid)
  
  p1 <- ggplot(diagData, aes(predPpnAboveBM, PearResid)) +
    geom_point(size=3) + 
    #geom_smooth(method="lm", formula=y~x) + 
    geom_smooth(method="loess", formula=y~x, span=1) + 
    geom_hline(linetype="dashed", yintercept=0) +
    xlab("Predicted proportions") + ylab("Pearson's Residuals") +
    ggtitle("Pearson's Residuals") +
    theme_classic() + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(size = 20)
    ) 
  
  p2 <- ggplot(diagData, aes(predPpnAboveBM, DevResid)) +
    geom_point(size=3) + 
    #geom_smooth(method="lm", formula=y~x) + 
    geom_smooth(method="loess", formula=y~x, span=1) + 
    geom_hline(linetype="dashed", yintercept=0) +
    xlab("Predicted proportions") + ylab("Deviance Residuals") +
    ggtitle("Deviance Residuals") +
    theme_classic() + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(size = 20)
    ) 
  
  # See ggplot.cor function in "helperFunctions.r"
  p3 <- ggplot.corr(data=PearResid, title="Pearsons's residuals") 
  p4 <- ggplot.corr(data=DevResid, title="Deviance residuals") 
  
  ggsave(filename=paste0(dir, "/",plotname, ".png"), p1+p2+p3+p4)
  #ggsave(filename=paste0(dir, plotname, "1.png"), p2)
  #ggsave(filename=paste0(dir, plotname, "2.png"), p4)
  
  
  #-------------------------------------------------------------------------------
  # step 3. 
  # Evaluate statistical significance of model coefficients.
  #-------------------------------------------------------------------------------
  
  signTable <-  All_Ests %>% filter(Param %in% c("B_0", "B_1")) %>% 
    rename(P.value=Pr...z.2..) %>% dplyr::select(Param, Estimate, 
                                                 Std..Error, z.value, 
                                                 P.value)
  # Look at table here:
  # knitr::kable(signTable)
  
  #-------------------------------------------------------------------------------
  # Step 4. 
  # Evaluate Pearson Chi-square statistic (goodness of fit). 
  #   Is there statistical evidence for lack of fit?
  #-------------------------------------------------------------------------------
  
  # Evaluate goodness of fit by comparing the residual deviance to a Chi-square
  # distribution  
  
  # Sum of squared Pearson's residuals
  Pearson <- sum(PearResid^2)
  
  # Statistical test of the goodness of fit 
  # Section 6.5.6. https://bookdown.org/roback/bookdown-bysh/ch-logreg.html
  # Section 1.4 of https://www.flutterbys.com.au/stats/tut/tut10.5a.html
  p.PearChiSq <- 1 - pchisq(q=Pearson, df=length(PearResid)-2)
  #values < 0.05 indicate statistically significant evidence for lack of fit
  
  
  #-------------------------------------------------------------------------------
  # Step 5. 
  # Evaluate Deviance G-squared statistic (goodness of fit)
  #   Is there statistical evidence for lack of fit?
  #-------------------------------------------------------------------------------
  
  # Deviance statistic (sum of deviance residuals)
  # The deviance is a key concept in logistic regression. It measures the
  # deviance of the fitted logistic model with respect to a perfect model (the 
  # saturated model) https://bookdown.org/egarpor/SSS2-UC3M/logreg-deviance.html
  Deviance <- sum(DevResid^2)
  p.DevChiSq <- 1-pchisq(q=Deviance, df=length(DevResid)-2)
  #values < 0.05 indicate statistically significant evidence for lack of fit
  # See section 6.5.6: https://bookdown.org/roback/bookdown-bysh/ch-logreg.html
  #-------------------------------------------------------------------------------
  # Step 6. 
  # Evaluate quasi-Rsquared.
  #   What is the ratio of the fit to the null model?
  #-------------------------------------------------------------------------------
  
  # What is ratio of fit to the null model? This is a measure of the
  # strength of the relationship
  NullDev <- deviance(glm( obsPpnAboveBM ~ 1 , family = 
                             quasibinomial))
  quasiR2 <- 1 - Deviance/NullDev
  # This is not the percentage of variance explained by the logistic model, 
  # but rather a ratio indicating how close is the fit to being perfect 
  # or the worst. It is not related to any correlation coefficient
  # https://bookdown.org/egarpor/SSS2-UC3M/logreg-deviance.html
  # However, shen we forced the y-intercept near zero,  the model actually fits 
  # worse than the null model and quasiR2 can be negative
  
  
  #-------------------------------------------------------------------------------
  # Step 7. 
  # Evaluate Wald test
  #   Is the predictor, 'Aggregate Abundance', a significant predictor of the ppn
  #   of CUs > their lower benchmark?
  #-------------------------------------------------------------------------------
  
  # The Wald test evaluates the significance of a predictor based on difference 
  # in Deviances between the specificed model and a reduced model without the 
  # predictor
  # See Section 4.7 https://bookdown.org/egarpor/SSS2-UC3M/logreg-deviance.html
  p.Wald <- signif( pchisq(q=(Deviance-NullDev), df=1), digits=2)
  
  
  #-------------------------------------------------------------------------------
  # Step 8. 
  # Evaluate hit rate from a confusion matrix
  #   What is the classification accuracy of the LRP based on the logistic 
  #   regression?
  #-------------------------------------------------------------------------------
  
  
  # In which years did the model predict aggregate abundances >LRP?
  yHat <- predPpnAboveBM > p
  # In which years were observed aggregate abundances >LRP?
  y <- obsPpnAboveBM > p
  
  # Confusion Matrix
  confMat <- table(y, yHat)
  
  # What is the accuracy in classifying observed aggregate abundances?
  # Hit ratio = ratio of correct classification
  hitRatio <- sum(diag(confMat))/sum(confMat)
  hitRatio <- round(hitRatio, digits=2)
  
  # The hit matrix will be always biased towards unrealistically good
  # classification rates when computed with the same sample used for
  # fitting the logistic model. An approach based on data-splitting/
  # cross-validation is therefore needed to estimate the hit matrix 
  # unbiasedly  https://bookdown.org/egarpor/PM-UC3M/glm-modsel.html
  
  #-------------------------------------------------------------------------------
  # Output
  out <- list()
  out$DevResid <- DevResid
  out$signTable <- signTable
  out$p.PearChiSq <- p.PearChiSq 
  out$p.DevChiSq <- p.DevChiSq 
  out$quasiR2 <- quasiR2
  out$p.Wald <- p.Wald
  out$confMat <- confMat
  out$hitRatio <- hitRatio
  
  return(out)
  
  #-------------------------------------------------------------------------------
  
  
} # End of Function 1: LRdiagnostics()


#-------------------------------------------------------------------------------
# Function 2: LOO_LRdiagnostic, performs leave-one-out cross-validation by 
# iteratively re-running LRP estimation of subsets of data, leaving one data 
# point out each time, and estimating classification accuracy on removed data point

# Steps:
# 1. Run logistic regression iteratively removing one year in time-series 
# each time
# 2. Get observed ppn of CUs above their lower benchmark for the year that
# was held out (iteratively)
# 3. Get predicted ppn of CUs above their lower benchmark for the year that
# was held out (iteratively)
# 4. Evaluate confusion matrix hit rate on the resulting time-series of 
# predictied and observed ppns.
#   What is the out-of-sample classification accuracy?

# Note, this function will have to be adapted to case studies based no the 
# logistic regression model used for those stocks. I revised my Get.LRP() 
# function by adding an argument LOO, where LOO = numeric for leave-one-out
# cross validation of the logistic regression. This number is the index of 
# the time-series of ppn of CUs and aggregate abundances that are removed
# prior to implementing the logistic regression in TMB. It is set to NA as 
# a default (no values removed). Note, the outputted time-series ('out') 
# contain all the data, but parameter estimates are derived from 
# time-series without LOO index value. Within the Get.LRP function, I added 
# the following code when setting up data list for TMB:
# if(!is.na(LOO)) {
#  data$LM_Agg_Abund <- data$LM_Agg_Abund[-LOO]
#  data$N_Above_BM <- data$N_Above_BM[-LOO]
# }



#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Arguments:
# remove.EnhStocks = logical specifying if enhanced stocks are removed 
  # (PNI<0.5)
# n = length of time-series over with to iteravely remove a data point in the 
  # logistic regression


# Note, LOO_LRdiagnists() calls the function Get.LRP() which generates the 
# following outputs required for calculation of this diagnostic:
# out$All_Ests = Dataframe containing parameters of logistic regression, B_0 and B_1
  # with p values (All_Ests <- data.frame(summary(sdreport(obj), p.value=TRUE)))
# SMU_Esc = Vector of scaled observed aggregate abundances, with NAs removed
# out$Logistic_Data$yy = Vectors of observed ppn of CUs above their lower benchmarks
  # with NAs removed
# LRPppn = the proportion of CUs that must be > their benchmark when defining LRP


#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Returns:
# - Hit Ratio from confusion matrix, based on leave-one-out predictions
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------
# Function 2:

LOO_LRdiagnostics <- function(remove.EnhStocks=TRUE, n=18){
  
  
  # Step 1: estimate logistic regression iteratively, removing a single year 
  # each time
  
  predPpnAboveBM <- NA
  
  for (i in 1:n){
    # Estimate logistic regression using function Get.LRP
    zz <- Get.LRP(remove.EnhStocks = TRUE, LOO=i)
    All_Ests <- zz$out$All_Ests

    if(i==1){ # These remain constant over iterations
      # Step 2: Get observed time-series of aggregate raw abundances that includes all
      # data and then scale to units near 1-10 
      AggAbundRaw <- zz$out$Logistic_Data$xx
      digits <- count.dig(AggAbundRaw)
      ScaleSMU <- min(10^(digits -1 ), na.rm=T)
      AggAbund <- AggAbundRaw/ScaleSMU
      # Get time-series of observed ppns of CUs> benchamark, including all
      # data
      obsPpnAboveBM <- zz$out$Logistic_Data$yy
      # Get threshold p value (ppn of CUs>benchmark) used to estimate LRP
      p <- zz$LRPppn
      #dir <- "DataOut/"
    }
 
   # Step 3: Get predicted ppn of CUs above their lower benchmark for the year 
    # that was held out
    B_0 <- All_Ests %>% filter(Param=="B_0") %>% pull(Estimate)
    B_1 <- All_Ests %>% filter(Param=="B_1") %>% pull(Estimate)
    #predPpnAboveBM <- inv_logit(B_0 + B_1*AggAbund)
    predPpnAboveBM[i] <- inv_logit(B_0 + B_1*AggAbund[i])
  } # End of for i in 1:18
  
  # Step 4: Calculate Hit Ratio
  
  # In which years did the model predict aggregate abundances >LRP?
  yHat <- predPpnAboveBM > p
  # In which years were observed aggregate abundances >LRP?
  y <- obsPpnAboveBM > p
  
  # Confusion Matrix
  confMat <- table(y, yHat)
  
  # What is the accuracy in classifying observed aggregate abundances?
  # Hit ratio = ratio of correct classification
  hitRatio <- sum(diag(confMat))/sum(confMat)
  hitRatio <- round(hitRatio, digits=2)
  
  return(hitRatio=hitRatio)
}# End of Function 2: LOO_LRdiagnostics()




#NB the logit function from Brooke's helperFunction code needs an extra set of
# ()s. I don't think it was actually implemented in the code...?

logit <- function(x){
  log(x/(1-x))
}
