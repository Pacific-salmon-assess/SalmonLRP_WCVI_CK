# Inverse logit and logit funcitons can come in handy =====================================================================
inv_logit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  log(x/1-x)
}



gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

stdErr <- function(x) sd(x)/sqrt(length(x))

# Crop data: Filter to only include years with data in all years between startYr and endYr for MU =======================

cropData<-function(Dat,startEndYrs, MU_ID) {
  # Extract start and end years for that MU  
  startYr<-startEndYrs[as.character(startEndYrs$MU_ID)==MU_ID,"LRPStartYr"]
  endYr<-startEndYrs[as.character(startEndYrs$MU_ID)==MU_ID,"LRPEndYr"]
  
  # Filter to only include years with data in all years between startYr and endYr for MU
  Dat.MU <- Dat%>% filter(MU == MU_ID & yr >= startYr & yr <= endYr)
  required_nYrs<-length(startYr:endYr)
  nYrs_byCU <- Dat.MU %>% group_by(CU)  %>% summarize(nYrs = length(na.omit(ets))) 
  
  Dat.MU<-left_join(Dat.MU, nYrs_byCU,by = "CU") %>% filter(nYrs == required_nYrs)
  
  Dat.MU
 
}