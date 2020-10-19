
# Functions in this file were taken from "Functions.R" (from "Retrospective Analysis BD" folder) 
#    written by B. Davis for Holt et al. 2018 Chum CSAS paper



##################
###  Infill   ###
#################

# Takes data frame, variable to group by, and infills, also takes unique identifier for sites (rabcode, here)
# Must have column of "Escape"
# Can do site average as geometric mean or arithmatic mean according to avg

Infill <- function( data, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME", avg="GeoMean", EscCol="Escape"){
  
  #Make sure column with Esc is called "Escape"
  names(data)[which(names(data)==EscCol)] <- "Escape"
  
  # Calculate geometric mean of each unit within group
  if(avg=="GeoMean"){
    #SiteAvg <- data %>% group_by_(.dots=c(unit, groupby, Uid)) %>% mutate(EscAvg = GeoMean(Escape))
    #KH - changed to:
    SiteAvg <- data %>% group_by(.dots=c(unit, groupby, Uid)) %>% mutate(EscAvg = GeoMean(Escape))
  }
  if(avg=="Amean"){
    #SiteAvg <- data %>% group_by_(.dots=c(unit, groupby, Uid)) %>% mutate(EscAvg = mean(Escape, na.rm=T))
    #KH - changed to:
    SiteAvg <- data %>% group_by(.dots=c(unit, groupby, Uid)) %>% mutate(EscAvg = mean(Escape, na.rm=T))
  }
  
  #need Year to be character string for group_by_ and .dots notation to work
  Year="Year"
  # proportion of the group that this unit average makes up
  #GroupAvg <- SiteAvg %>% group_by_(.dots=c(groupby, Year)) %>% mutate(GroupSum = sum(EscAvg))
  # KH - change to:
  GroupAvg <- SiteAvg %>% group_by(.dots=c(groupby, Year)) %>% mutate(GroupSum = sum(EscAvg))
  GroupAvg$Props <- GroupAvg$EscAvg/GroupAvg$GroupSum
  
  #go through and add up all props for sites with counts that aren't NA
  GroupAvg$Present <- ifelse(is.na(GroupAvg$Escape), 0,1)
  GroupAvg$PresentProps <- GroupAvg$Present*GroupAvg$Props
  
  # go through each group and calculate 1/(present proportions) for that year
  #AllDat <- GroupAvg %>% group_by_(.dots=c(groupby, Year)) %>% mutate(F=1/sum(PresentProps))
  # KH - changed to:
  AllDat <- GroupAvg %>% group_by(.dots=c(groupby, Year)) %>% mutate(F=1/sum(PresentProps))
  
  # Now for each group total is F*(total Escapement)
  #AllDat2 <- AllDat %>% group_by_(.dots=c(groupby, Year)) %>% mutate(GroupCount=sum(Escape, na.rm=T)) 
  #AllDat3 <- AllDat2 %>% group_by_(.dots=c(groupby, Year)) %>% mutate(GroupEsc=F*sum(Escape, na.rm=T))
  # KH - changed to:
  AllDat2 <- AllDat %>% group_by(.dots=c(groupby, Year)) %>% mutate(GroupCount=sum(Escape, na.rm=T)) 
  AllDat3 <- AllDat2 %>% group_by(.dots=c(groupby, Year)) %>% mutate(GroupEsc=F*sum(Escape, na.rm=T))
  # calculate estimated contributions for missing data
  #AllDat4 <- AllDat3 %>% group_by_(.dots=c(groupby, Year)) %>% mutate(ContrEsc=GroupEsc*Props)
  # KH - changed to:
  AllDat4 <- AllDat3 %>% group_by(.dots=c(groupby, Year)) %>% mutate(ContrEsc=GroupEsc*Props)
  #fill in NA's with ContrEsc
  AllDat4$SiteEsc <- ifelse(is.na(AllDat4$Escape), AllDat4$ContrEsc, AllDat4$Escape)
  
  # Now summarise so just see each group, Year, will use mean arbitrarily since all the same
  #ByGroup <- AllDat4 %>% group_by_(.dots=c(groupby, Year))  %>% summarise(GroupEsc=mean(GroupEsc))
  # KH changed to:
  ByGroup <- AllDat4 %>% group_by(.dots=c(groupby, Year))  %>% summarise(GroupEsc=mean(GroupEsc))
  
  #get rid of NaN's
  ByGroup$GroupEsc[which(is.nan(ByGroup$GroupEsc)==T)] <- NA
  
  # Return with data by site, and also by group
  list( AllDat4, ByGroup)
  
}

#################
## %notin%    ##
###############

#Usefull operator works like !%in%, which doesn't work

`%notin%` <- function(x,y) { !(x %in% y) }



################
###  GeoMean  ##
###############

#Takes vector and gives geometric mean, ignoring NA's

GeoMean <- function(x){
  xx <- x[which(is.na(x)==F)]
  exp(mean(log(xx)))
}
