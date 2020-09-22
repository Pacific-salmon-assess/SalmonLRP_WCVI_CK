# Trying to estimate SMSY using Lambert's W instead of Hilborn and Walters eqn.

library(TMB)
library(dplyr)

# use data from retro evaluation
codeDir<-getwd()
setwd('..')
tmpDir<-getwd()
cohoDir<-paste(tmpDir,"/IFCohoStudy",sep="")
setwd(cohoDir)
Ests <- read.csv("DataOut/Bin_095_Esc/annualRetro__SRparsByCU.csv")
setwd(codeDir)

dyn.unload(dynlib("TMB_Files/Lambert_W"))
compile("TMB_Files/Lambert_W.cpp")

dyn.load(dynlib("TMB_Files/Lambert_W"))

# pull out alpha and beta pair

A <- Ests %>% filter(retroYr == 2018 ) %>% pull(est_A) 
b <- Ests %>% filter(retroYr == 2018 ) %>% pull(est_B)
SMSY <- Ests %>% filter(retroYr == 2018) %>% pull(est_Smsy)

data <- list()
data$a <- log(A)
data$b <- b
param <- list()
param$dummy <- 0

obj <- MakeADFun(data, param, DLL="Lambert_W", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr)
summary(sdreport(obj))

# Check that they are actually better
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))

SMSY_Ws <- All_Ests %>% filter(Param == "SMSY_W") %>% pull(Estimate)
SMSY_Hs <- All_Ests %>% filter(Param == "SMSY_Hil") %>% pull(Estimate)

# calculate Yield for each
 Y_Ws <- SMSY_Ws * exp(log(A) - b*SMSY_Ws) - SMSY_Ws
 Y_Hs <- SMSY_Hs * exp(log(A) - b*SMSY_Hs) - SMSY_Hs
 # Yield slightly better using W method! (less than one fish though)

