
# Functions in this file were taken from "Functions.R" (from "Retrospective Analysis BD" folder) 
#    written by B. Davis for Holt et al. 2018 Chum CSAS paper



##################
###  Infill   ###
#################

# Takes data frame, variable to group by, and infills, also takes unique identifier for sites (rabcode, here)
# Must have column of "Escape"
# Can do site average as geometric mean or arithmetic mean according to avg

Infill <- function( data, groupby=c("CU_Name"), Uid = c("Rabcode", "GroupName"), unit="NME", avg="GeoMean", EscCol="Escape"){
  
  #Make sure column with Esc is called "Escape"
  names(data)[which(names(data)==EscCol)] <- "Escape"
  
  # Calculate geometric mean of each unit (stream) within group, over all years, and add column. The average escapement for each steam, over all years. 
  if(avg=="GeoMean"){
    SiteAvg <- data %>% group_by(.dots=c(unit, groupby, Uid)) %>% mutate(EscAvg = GeoMean(Escape)) # groups by NME (stream), CU, and Uid (combo of Rabcode and Group name), adds column with geometric mean of escapement across all years
  }
  if(avg=="Amean"){
    SiteAvg <- data %>% group_by(.dots=c(unit, groupby, Uid)) %>% mutate(EscAvg = mean(Escape, na.rm=TRUE))
  }
  
  Year="Year" #need Year to be character string for group_by and .dots notation to work
  
  # proportion of the group that this unit average makes up
  GroupAvg <- SiteAvg %>% group_by(.dots=c(groupby, Year)) %>% mutate(GroupSum = sum(EscAvg)) # groups data frame from above by CU and year, and sums the average escapement. Sum of average stream escapement for each CU.
  GroupAvg$Props <- GroupAvg$EscAvg/GroupAvg$GroupSum # Proportion of average escapement for each stream from its group (here, CU). The contribution of each stream to its CU. Or, average stream escapement / sum of average stream escapement across CU
  
  #go through and add up all props for sites with counts that aren't NA
  GroupAvg$Present <- ifelse(is.na(GroupAvg$Escape), 0,1) # does stream have escapement number?, if yes, 1, if not, 0
  GroupAvg$PresentProps <- GroupAvg$Present*GroupAvg$Props # multiply proportion of each stream by whether it is present or not
  
  # go through each group and calculate 1/(present proportions) for that year. Value = 1 if all streams were counted, increases with number of uncounted streams
  AllDat <- GroupAvg %>% group_by(.dots=c(groupby, Year)) %>% mutate(inverse_prop=1/sum(PresentProps)) # Note that if there are no steams counted within a CU in a year, this becomes infinity
  
  # Now for each group total is F*(total Escapement)
  AllDat2 <- AllDat %>% group_by(.dots=c(groupby, Year)) %>% mutate(GroupCount=sum(Escape, na.rm=TRUE)) # add column with sum of escapement by group (CU) and year
  AllDat3 <- AllDat2 %>% group_by(.dots=c(groupby, Year)) %>% mutate(GroupEsc=inverse_prop*sum(Escape, na.rm=TRUE)) # multiply factor (1/PresentProps) by group:year escapement total. Expands CU escapement by the missing proportion of streams. If all streams were counted, it is the same as the sum of escapement. 
  # calculate estimated contributions for missing data by multiplying the expanded group escapement by the stream-specific proportion
  AllDat4 <- AllDat3 %>% group_by(.dots=c(groupby, Year)) %>% mutate(ContrEsc=GroupEsc*Props)
  #fill in NA's with ContrEsc
  AllDat4$SiteEsc <- ifelse(is.na(AllDat4$Escape), AllDat4$ContrEsc, AllDat4$Escape) # new column. if stream was counted, use actual value. If it wasn't counted, use estimated value based on proportions
  
  # Now summarise so just see each group, Year, will use mean arbitrarily since all the same
  # and add a column that is actual counts per CU
  ByGroup <- AllDat4 %>% group_by(.dots=c(groupby, Year))  %>% summarise(GroupEsc=mean(GroupEsc), SumRawEsc = sum(Escape, na.rm=TRUE) )
  
  #get rid of NaN's
  ByGroup$GroupEsc[which(is.nan(ByGroup$GroupEsc)==TRUE)] <- NA
  
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
  xx <- x[which(is.na(x)==FALSE)]
  exp(mean(log(xx)))
}
