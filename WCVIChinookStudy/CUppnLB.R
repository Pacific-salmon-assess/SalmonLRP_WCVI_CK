# Code to derive data file for whether individual CUs are above LB
# from samSimOutputs
# 2 April 2024

# setwd('..')
rootDir<-getwd()
codeDir<-paste(rootDir,"/Code",sep="")
wcviCKDir<-paste(rootDir,"/WCVIChinookStudy",sep="")
outDir <- wcviCKDir
OM <- "baseER"
SA <- "_AllExMH"
simPars<-read.csv(paste(outDir, "SamSimInputs/SimPars.csv",sep="/"))
# CU_bench <- readRDS("C:/github/SalmonLRP_RetroEval/WCVIChinookStudy/SamSimOutputs/simData/baseER_CoreInd/baseER_CoreInd_baseER_CUaboveLB.RData")
CU_bench <- readRDS(paste(outDir, "/SamSimOutputs/simData/", OM, SA, "/", OM, SA ,"_", OM, "_CUaboveLB.RData", sep=""))
# Get sAg from lrpDat.csv?
sAgDat <- read.csv(paste(outDir, "/SamSimOutputs/simData/", OM, SA, "/", OM, SA ,"_", OM, "_lrpDat.csv", sep=""))



# Need to concatenate CULRP data, to have similar strucure at lrpDat.csv, columns: year, interation, sAg, ppnCUsLowerBM
# I.e make long data frame

nYears <- dim(CU_bench)[1]
year <- data.frame(year=1:nYears)

# need to pivot_longer first, putting array elements iterations) in rows (currently years)
# Change array to list of dataframes then add column of years and iteration to each list element,
# and the concatenate list elements as additional rows in a df
nTrials <- dim(CU_bench)[3]
CU_bench_list <- list()
for (z in 1:nTrials){
  dum <- data.frame(CU_bench[,,z])
  dum <- dum %>% tibble::add_column(year)
  dum <- dum %>% tibble::add_column(data.frame(iteration=z))
  CU_bench_list[[z]] <- dum
}

for (z in 1:nTrials){
  if(z == 1) CU_bench_long <- as.data.frame(CU_bench_list[[z]])
  if(z > 1) CU_bench_long <- CU_bench_long %>%
      tibble::add_row( as.data.frame(CU_bench_list[[z]]) )
}


# Tests:
# Then select column sAg making sure nYears and nTrials are the same as CU_bench
if(length(unique(sAgDat$year)) != length(pull((unique(CU_bench_long["year"])), year))) {
  print("Error. CU_bench has differnt dimensions than SAg")
}
if(length(unique(sAgDat$iteration)) != length(pull((unique(CU_bench_long["iteration"])), iteration))) {
  print("Error. CU_bench has differnt dimensions than SAg")
}

# Join sAG and CU_bench_long
CU_bench_long <- CU_bench_long %>%
  left_join(sAgDat, by=join_by(year, iteration))
CU_bench_long <- CU_bench_long %>%
  select(!ppnCUsLowerBM)


# read to a csv file
filename <- paste("projCUBenchDat_",OM, SA,".csv",sep="")
write.csv(CU_bench_long, file=paste(outDir, "/SamSimOutputs/simData/", filename, sep="") )

# This format would be the same as "_projLRPDat_baseER_CoreInd.csv"
# where columns are wethere each CU 1:K was above lower benchmark

