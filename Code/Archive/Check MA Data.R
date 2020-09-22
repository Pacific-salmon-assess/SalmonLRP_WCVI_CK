codeDir<-getwd()
setwd('..')
tmpDir<-getwd()
cohoDir<-paste(tmpDir,"/IFCohoStudy",sep="")
setwd(cohoDir)
CoEscpDat <- read.csv("DataIn/IFCoho_escpByCU.csv")
CoSRDat <- read.csv("DataIn/IFCoho_SRbyCU.csv")
MA_Esc <- read.table("DataIn/SpRecByCalYr_Arbeider18.dat", header = T)
MA_SR <- read.table("DataIn/SR_Arbeider18_Type1.dat", header = T)
setwd(codeDir)


 # check that data is the same

Esc_Join <- left_join(CoEscpDat, MA_Esc, by = c("CU_Name" = "CUnm", "ReturnYear" = "Year"))
Esc_Join %>% filter(is.na(Sp) == F & Escapement != Sp)

# there are some other years htat are different -- ask Michael why is his files don't go back to 1984

# check S-R data
#start with just SR
SR_Join <- left_join(CoSRDat[, c(1,3,4,5)], MA_SR[, c(1, 3, 4, 5)], by = c("CU_Name" = "CUnm", "BroodYear" = "Byr"))
SR_Join %>% filter(Spawners != Sp) # all match
SR_Join %>% filter(Recruits != Rec_total) # don't match

# Look if esc in MA's SR df match those in the spawner time series
MA_Join <- left_join(MA_Esc, MA_SR[c(1,3,4)], by = c("CUnm", "Year" = "Byr"))
MA_Join %>% filter(Sp.x != Sp.y)

