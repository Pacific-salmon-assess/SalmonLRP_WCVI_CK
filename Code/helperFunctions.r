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



ggplot.corr <- function(data, lag.max = 10, ci = 0.95, title="") {
  # Adapted from https://rh8liuqy.github.io/ACF_PACF_by_ggplot2.html
  
  list.acf <- acf(data, lag.max = lag.max, type = "correlation", plot = FALSE)
  N <- as.numeric(list.acf$n.used)
  df1 <- data.frame(lag = list.acf$lag, acf = list.acf$acf)
  df1$lag.acf <- dplyr::lag(df1$acf, default = 0)
  df1$lag.acf[2] <- 0
  df1$lag.acf.cumsum <- cumsum((df1$lag.acf)^2)
  df1$acfstd <- sqrt(1/N * (1 + 2 * df1$lag.acf.cumsum))
  df1$acfstd[1] <- 0
  df1 <- df1 %>% dplyr::select(lag, acf, acfstd)
  
  plot.acf <- ggplot(data = df1, aes( x = lag, y = acf)) +
    geom_col(fill = "#4373B6", width = 0.7) +
    geom_hline(yintercept = qnorm((1+ci)/2)/sqrt(N), 
               colour = "sandybrown",
               linetype = "dashed") + 
    geom_hline(yintercept = - qnorm((1+ci)/2)/sqrt(N), 
               colour = "sandybrown",
               linetype = "dashed") + 
    scale_x_continuous(breaks = seq(0,max(df1$lag),6)) +
    scale_y_continuous(name = element_blank(), 
                       limits = c(min(df1$acf, - qnorm((1+ci)/2)/sqrt(N)) ,1)) +
    ggtitle(paste("ACF for", title)) +
    theme_bw() + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(size = 14))
  
  
  return(plot.acf)
  
}
