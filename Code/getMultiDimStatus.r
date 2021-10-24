
library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(here)
library(zoo)



getMultiDimStatus<-function(statusDat,outDir) {
 
  rapidStatus<-NA
  
  for (i in 1:nrow(statusDat)) {
    
    if (statusDat$Gen_Mean[i] < 1500) {
      rapidStatus[i] <- "Red" }
    else {
    
      if (statusDat$Gen_Mean[i] > 10000) {
        rapidStatus[i]<-ifelse(statusDat$Gen_Mean[i] < statusDat$Sgen[i], "Red", ifelse((statusDat$Gen_Mean[i]*1.1) > statusDat$Smsy0.8[i],"Green","Amber"))
      } else {
        rapidStatus[i]<-ifelse(statusDat$Gen_Mean[i] < statusDat$Sgen[i], "Red", "Amber")
      }
    }
    
  }
  
  statusOut<-statusDat %>% add_column(rapidStatus=rapidStatus)
  
  write.csv(statusOut, paste(outDir, "DataOut/multiDimStatusEsts.csv", sep="/"), row.names=F)
  
  return(statusOut)
  
}