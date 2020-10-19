
# load packages
library(dplyr)
library(tidyr)

# Read in Return data from P.V.W received 3/7/2016 =========================
WildRetWide <- read.csv("DataIn/WildReturnsPVW_2013.csv", check.names=F)

# Need to get PVW return data into long form
# First need to "collapse" areas to CU level
WildByCU <- WildRetWide  %>% group_by(CU_Name)   %>%  summarise_at(vars("2013": "1953"),sum,na.rm = T)
WildByCU$CU_Name <- as.character(WildByCU$CU_Name)



# Now into long form
WildRetLong <- WildByCU %>% gather("Year", "Return", 2:62)


# Read in Return data from P.V.W received in 2020 =======================
RetWide <- read.csv("DataIn/WildReturnsPVW_2018.csv", check.names=F)

# Need to get PVW return data into long form
# First need to "collapse" areas to CU level
ByCU <- RetWide  %>% group_by(CU_Name)  %>%  summarise_at(vars("2018": "1954"),sum,na.rm = T)
ByCU$CU_Name <- as.character(ByCU$CU_Name)

# Now into long form
RetLong <- ByCU %>% gather("Year", "Return", 2:66)

CUNames<-unique(RetLong$CU_Name)

# Loop over CUs to plot

par(mfrow=c(3,3),mar=c(2,2,1,1))

for (i in 1:length(CUNames)) {
  
  dat.2013<-WildRetLong %>% filter(CU_Name == CUNames[i])
  dat.2018<-RetLong %>% filter(CU_Name == CUNames[i])
  
  y.max<-max(c(dat.2013$Return,dat.2018$Return))
  
  plot(dat.2013$Year,dat.2013$Return,typ="l",col="black", xlim=c(1953,2018),ylim=c(0,y.max), main=CUNames[i])
    points(dat.2013$Year,dat.2013$Return, pch=19)
    points(dat.2018$Year,dat.2018$Return, col="coral")
    lines(dat.2018$Year,dat.2018$Return, col="coral", lty=1)
  
}

plot(0, xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')
legend(x=0, legend=c("2013 data", "2018 data"), pch=c(19,1), lty=1, col=c("black", "coral")) 


