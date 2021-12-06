
# **** Plotting functions for looking at data availability ***********
# plot_Num_CUs_Over_Time()      - Plot number of CUs over time for all MUs
# plot_CU_DataObs_Over_Time()   - Plot available CU estimates over time, by MU
# plot_CU_Escp_Over_Time()
# plot_CU_EscpRatio_Over_Time()



# *** Plotting functions that show retrospective results *********

# plotAnnualRetro() 
# plotLogistic()


library(cowplot)


# =================================================================================
# Functions for plotting data availability
# ================================================================================

# Plot number of CUs over time for all MUs
plot_Num_CUs_Over_Time<-function(Dat,Dir, plotName) {
  
  pdf(paste(Dir,"/Figures/",plotName,".pdf",sep=""))
  
  # How many stocks per CU?
  Max_CUs <- Dat %>% group_by(MU) %>% summarise(n=length(unique((CU))))
  
  # Num_CUs_Over_Time <- Dat %>%  filter(is.na(Escp) == F) %>% group_by(MU, yr) %>% summarise(n=length(unique((CU))))
  # LW - changed above line to: (is.na(Escp) was giving unexpected results with chum data. Maybe a filter/numeric object thing?)
  Num_CUs_Over_Time <- Dat %>%  filter(!is.na(Escp)) %>% group_by(MU, yr) %>% summarise(n=length(unique((CU))))
  
  MUs <- unique(Num_CUs_Over_Time$MU)
  for(mm in 1:length(MUs)){
    Dat <- Num_CUs_Over_Time %>% filter(MU == MUs[mm])
    if(mm==1){
      plot(Dat$yr, Dat$n, col = mm, type="l", ylim = c(1, max(Max_CUs$n)), lwd=2, 
           xlab = "Year", ylab = "Number of CUs", main = "CUs per MU Over Time")
    } else {
      lines(Dat$yr, Dat$n, col = mm, lwd = 2)
    }
    legend("topleft", lty=1, col = 1:length(MUs), legend = MUs, lwd=2)
  }
  
  dev.off()
}


# Plot available CU estimates over time, by MU
plot_CU_DataObs_Over_Time <-function(Dat, Dir, plotName) {
  
  theme_update(plot.title = element_text(hjust = 0.5))
  
  pdf(paste(Dir,"/Figures/",plotName,".pdf",sep=""))
  
  MUs <- unique(Dat$MU)
  
  p <- list()
  
  for(mm in 1:length(MUs)){
    Dat.MU <- Dat %>% filter(MU == MUs[mm])
    binDat <- Dat.MU %>% filter(MU == MUs[mm]) %>% group_by(CU, yr) %>% summarise(value=ifelse(Escp >= 0, 1, 0))
    p[[mm]]<-qplot(data=na.omit(binDat),x=yr,y=CU) + 
      labs(title=MUs[mm], x = "Year", y = "CU") + 
      theme(plot.title = element_text(size=14, face="bold",hjust=0.5))
  }
  
  do.call(grid.arrange,p)
  
  dev.off()
  
}

# Plot available CU estimates over time, by CU
plot_CU_Escp_Over_Time <-function(Dat, Dir, plotName, samePlot = T, withSgen=NULL, SgenFileName=NULL, addGenMean=TRUE) {
  
  theme_update(plot.title = element_text(hjust = 0.5))
  
  p <- list()
  
  if(samePlot == T){
    
    png(paste(Dir,"/Figures/",plotName,".png",sep=""), width=380, height=225)
    
    Dat$CU_Name[Dat$CU_Name == "Fraser_Canyon"] <- "Fraser Canyon"
    Dat$CU_Name[Dat$CU_Name == "Middle_Fraser"] <- "Middle Fraser"
    Dat$CU_Name[Dat$CU_Name == "North_Thompson"] <- "North Thompson"
    Dat$CU_Name[Dat$CU_Name == "South_Thompson"] <- "South Thompson"
    Dat$CU_Name[Dat$CU_Name == "Lower_Thompson"] <- "Lower Thompson"

      MUs <- unique(Dat$MU)
      for(mm in 1:length(MUs)){
        Dat.MU <- Dat %>% filter(MU == MUs[mm])
        escpDat <- Dat.MU %>% filter(MU == MUs[mm]) %>% group_by(CU_Name, yr) %>% summarise(Escp=sum(Escp))
        p[[mm]]<-ggplot(data=escpDat,aes(x=yr,y=Escp)) + 
              geom_line(aes(colour=CU_Name),size=1.2 ) + 
              labs(x = "Year", y = "Escapement") + 
              theme(plot.title = element_text(size=14, face="bold",hjust=0.5)) + 
              theme_classic()
        
      }
      
    
      do.call(grid.arrange, args = c(p, nrow=1))
        
  } else {
    
    png(paste(Dir,"/Figures/",plotName,".png",sep=""),width=780, height=510)
    
    
    CUs <- unique(Dat$CU_Name)
    
    if (withSgen == T) {
      # Get Sgen estimates for IM model
       ests<-read.csv(paste(Dir,"/DataOut/",SgenFileName,".csv", sep=""))
       Sgen_base<-ests[ests$Param=="Sgen","Estimate"]
       SRep_base<-ests[ests$Param=="SRep","Estimate"]
       cap.ests<-read.csv(paste(Dir,"/DataOut/",SgenFileName,"_priorCap.csv",sep=""))
       Sgen_cap<-cap.ests[cap.ests$Param=="Sgen","Estimate"]
       SRep_cap<-cap.ests[cap.ests$Param=="SRep","Estimate"]
      }
    
    for(ss in 1:length(CUs)){
    
      CU_label<-str_replace(CUs[ss], "_", " ")
      
       escpDat <- Dat %>% filter(CU_Name == CUs[ss]) %>%
         mutate(Gen_Mean = rollapply(Escp, 3, gm_mean, fill = NA, align="right"))
      
       p[[ss]]<-ggplot(data=escpDat,aes(x=yr,y=Escp)) +
        geom_line(col="darkgrey", size=1) +
        labs(title=CU_label, x = "Year", y = "Escapement") + 
        theme_classic() + ylim(0,28000) +
        theme(plot.title = element_text(size=16, face="bold",hjust=0.5),
              axis.text = element_text(size=12), axis.title =element_text(size=15)) 
    
       if (addGenMean ==TRUE) 
       {
         
         p[[ss]]<- p[[ss]] + geom_line(data=escpDat, aes(x=yr,y=Gen_Mean), size=1.4)
         
       }
       
       
       if (withSgen ==TRUE) 
       {
         
         p[[ss]]<- p[[ss]] + 
           geom_hline(yintercept=Sgen_base[[ss]], linetype="solid", color = "red") +
           geom_hline(yintercept=Sgen_cap[[ss]], linetype="dashed", color = "red") #+
         #geom_hline(yintercept=SRep_base[[ss]], linetype="solid", color = "blue") +
         #geom_hline(yintercept=SRep_cap[[ss]], linetype="dashed", color = "blue")
         
       }
       

    }
    
    do.call(grid.arrange, args = c(p, nrow=2))
    
  }

  
  
  
  
  dev.off()
  
}




plot_CU_Escp_withStatus<-function(Dat, Dir, plotName, SgenFileName=NULL) {

  
  statusDat<- Dat %>% group_by(CU_Name) %>%
    mutate(Gen_Mean = rollapply(Escp, 3, gm_mean, fill = NA, align="right"))  
    
  # Get Sgen estimates for IM model
  ests<-read.csv(paste(Dir,"/DataOut/",SgenFileName,".csv", sep=""))
  ests<-as_tibble(ests)
  ests <- ests %>% filter(Param == "Sgen") %>% select(Estimate, CU_Name)
  ests <- rename(ests, BenchmarkBase = Estimate)
  
  cap.ests<-read.csv(paste(Dir,"/DataOut/",SgenFileName,"_priorCap.csv",sep=""))  
  cap.ests<-as_tibble(cap.ests)
  cap.ests <- cap.ests %>% filter(Param == "Sgen") %>% select(Estimate, CU_Name)
  cap.ests <- rename(cap.ests, BenchmarkCap = Estimate)
  
  statusDat<-left_join(x=statusDat, y= ests, by="CU_Name")
  statusDat<-left_join(x=statusDat, y= cap.ests, by="CU_Name")
  
  # Rename to make prettier for plotting
  statusDat$CU_Name[Dat$CU_Name == "Fraser_Canyon"] <- "Fraser Canyon"
  statusDat$CU_Name[Dat$CU_Name == "Middle_Fraser"] <- "Middle Fraser"
  statusDat$CU_Name[Dat$CU_Name == "North_Thompson"] <- "N. Thompson"
  statusDat$CU_Name[Dat$CU_Name == "South_Thompson"] <- "S. Thompson"
  statusDat$CU_Name[Dat$CU_Name == "Lower_Thompson"] <- "L. Thompson"
  

  statusBase <- statusDat %>% summarise(StatusBase=ifelse(Gen_Mean>BenchmarkBase, "Above", "Below"))
  statusBase<-statusBase$StatusBase
  
  statusCap <- statusDat %>% summarise(StatusCap=ifelse(Gen_Mean>BenchmarkCap, "Above", "Below"))
  statusCap<-statusCap$StatusCap
  
  statusDat<-statusDat %>% add_column(statusBase=statusBase, statusCap = statusCap)
  
  g<-ggplot(statusDat) +
    geom_path(aes(y=Escp, x=yr), colour='black', alpha=0.5) +
    geom_point(aes(y=Gen_Mean, x=yr, colour=statusBase)) +
    geom_hline(aes(yintercept=BenchmarkBase), colour="orange") +
    geom_hline(aes(yintercept=BenchmarkCap), colour="orange", linetype=2) +
    ylab("Escapement") + xlab("Year") +
    scale_colour_manual(guide = NULL, breaks = c("Above", "Below"), values=c("gray40", "red")) +
    facet_wrap(~interaction(CU_Name), scales = "free_y") +
    theme_classic()
  
  # Save plot
  ggsave(paste(Dir,"/Figures/",plotName,".png",sep=""), plot = g,
         width = 10, height = 5, units = "in") 
  
  
}





plot_Subpop_Escp_withStatus <- function(Dat, Dir, plotName) {
  
  
  escpDat <- Dat %>% group_by(Subpop_Name, CU_Name, yr) %>% summarise(Escp=sum(Escp)) %>% 
    mutate(Gen_Mean = rollapply(Escp, 3, gm_mean, fill = NA, align="right")) 
  
  statusDat <- escpDat %>%
    add_column(Status=ifelse(escpDat$Gen_Mean > 1000, "Above", "Below"))
 
  
  # Change some subpop names to make it prettier for plotting
  
  statusDat$Subpop_Name[statusDat$Subpop_Name == "Lower Middle Fraser"] <- "Lower Middle Fraser"
  statusDat$Subpop_Name[statusDat$Subpop_Name == "Lower North Thompson"] <- "Lower North Thomp."
  statusDat$Subpop_Name[statusDat$Subpop_Name == "Middle and Lower Shuswap Rivers"] <- "Mid. and L.Shuswap Rivers"
  statusDat$Subpop_Name[statusDat$Subpop_Name == "Shuswap Lake Tributaries"] <- "Shuswap Lake Tribs."
  statusDat$Subpop_Name[statusDat$Subpop_Name == "Middle North Thompson"] <- "Middle North Thomp."
  statusDat$Subpop_Name[statusDat$Subpop_Name == "Upper North Thompson"] <- "Upper North Thomp."
  
  statusDat$CU_Name[statusDat$CU_Name == "South_Thompson"] <- "South Thompson"
  statusDat$CU_Name[statusDat$CU_Name == "Fraser_Canyon"] <- "Fraser Canyon"
  statusDat$CU_Name[statusDat$CU_Name == "Middle_Fraser"] <- "Middle Fraser"
  statusDat$CU_Name[statusDat$CU_Name == "North_Thompson"] <- "North Thompson"
  statusDat$CU_Name[statusDat$CU_Name == "Lower_Thompson"] <- "Lower Thompson"
  
  #g<-ggplot(statusDat) +
  #  geom_path(aes(y=Escp, x=yr), colour='black', alpha=0.5) +
  #  geom_point(aes(y=Gen_Mean, x=yr, colour=Status)) +
  #  geom_hline(aes(yintercept=1000), colour="orange") +
  #  ylab("Escapement") + xlab("Year") +
  #  scale_colour_manual(guide = NULL, breaks = c("Above", "Below"), values=c("gray40", "red")) +
  #  facet_wrap(~interaction(Subpop_Name), scales = "free_y") +
  #  theme_classic()

  cus<-unique(statusDat$CU_Name)
  pp<-list()
  for(n in seq_along(cus)){

    statusDato <-statusDat %>% filter(CU_Name==cus[n]) 
    g<-ggplot(statusDato) +
      geom_path(aes(y=Escp, x=yr), colour='black', alpha=0.5) +
      geom_point(aes(y=Gen_Mean, x=yr, colour=Status)) +
      geom_hline(aes(yintercept=1000), colour="orange") +
      ylab(paste(cus[n],"\n Escapement")) + xlab("Year") +
      scale_colour_manual(guide = NULL, breaks = c("Above", "Below"), values=c("gray40", "red")) +
      facet_wrap(~interaction(Subpop_Name), scales = "free_y",ncol=1) +
      theme_classic()
    pp[[n]]<-g
  }



tres<- plot_grid(pp[[1]],pp[[4]],ncol=2,align="h")
dois<- plot_grid(pp[[3]],pp[[5]],NULL,ncol=2, rel_heights = c(.66, .33),
  align="h")
um <- plot_grid(pp[[2]],NULL,NULL,ncol=1,align="h",
  rel_heights = c(.36, .32,.32))

gf<-plot_grid(tres, dois, um, nrow=1,rel_widths = c(.33, .33,.165))



  # Save plot
  ggsave(paste(Dir,"/Figures/",plotName,".png",sep=""), plot = gf,
         width = 12, height = 5, units = "in") 
    
  
}






# Plot available CU estimates over time, by CU
plot_Subpop_Escp_Over_Time <-function(Dat, Dir, plotName, samePlot = T, withThresh=TRUE, addGenMean=TRUE) {
  
  theme_update(plot.title = element_text(hjust = 0.5))
  
  png(paste(Dir,"/Figures/",plotName,".png",sep=""),width=570, height=480)
  
  p <- list()
 
  if(samePlot == T){
    #escpDat <- Dat %>% group_by(Subpop_Name, yr) %>% summarise(value=sum(Escp))
    escpDat <- Dat %>% filter(Subpop_Name == Subs[ss]) %>% group_by(yr) %>% summarise(Escp=sum(Escp)) %>% 
      mutate(Gen_Mean = rollapply(Escp, 3, gm_mean, fill = NA, align="right"))
    p<-ggplot(data=escpDat,aes(x=yr,y=Escp)) + 
      geom_line(aes(colour=Subpop_Name),size=1.1 ) + 
      labs( x = "Year", y = "Escapement") + 
      theme(plot.title = element_text(size=12, face="bold",hjust=0.5)) + 
      theme_classic()
    print(p)
  } else {
    Subs <- unique(Dat$Subpop_Name)
    for(ss in 1:length(Subs)){
      #escpDat <- Dat %>% filter(Subpop_Name == Subs[ss]) %>% group_by(yr) %>% summarise(Escp=sum(Escp))
      escpDat <- Dat %>% filter(Subpop_Name == Subs[ss]) %>% group_by(yr) %>% summarise(Escp=sum(Escp)) %>% 
        mutate(Gen_Mean = rollapply(Escp, 3, gm_mean, fill = NA, align="right"))
      
      if (Subs[ss] == "Middle and Lower Shuswap Rivers") Subs[ss] <- "Middle + Lower Shuswap" 
      
        p[[ss]]<-ggplot(data=escpDat,aes(x=yr,y=Escp)) +
        geom_line(colour="darkgrey") +
        labs(title=Subs[ss], x = "Year", y = "Escapement") + 
        theme(plot.title = element_text(size=8, face="bold",hjust=0.5)) + 
        theme_classic()
    
        
      if (addGenMean == T) {
        p[[ss]]<-p[[ss]] + geom_line(data=escpDat, aes(x=yr,y=Gen_Mean), size=1.4)
        
      }
        
        
      if (withThresh==T) {
        p[[ss]]<-p[[ss]] + geom_hline(yintercept=1000,colour="red")
      }
      
    }
    do.call( grid.arrange, args = c(p, ncol=3))
  }
  
  
  dev.off()
  
}


plot_CU_Escp_withMultiStatus<-function(statusEsts,outDir,plotName) {
  
  # Rename to make prettier for plotting
  statusEsts$CU_Name[statusEsts$CU_Name == "Fraser_Canyon"] <- "Fraser Canyon"
  statusEsts$CU_Name[statusEsts$CU_Name == "Middle_Fraser"] <- "Middle Fraser"
  statusEsts$CU_Name[statusEsts$CU_Name == "North_Thompson"] <- "N. Thompson"
  statusEsts$CU_Name[statusEsts$CU_Name == "South_Thompson"] <- "S. Thompson"
  statusEsts$CU_Name[statusEsts$CU_Name == "Lower_Thompson"] <- "L. Thompson"
  
 # LRPStatus <- statusEsts %>% summarise(LRPStatus=ifelse(rapidStatus=="Red", "Below", "Above"))
#  status<-as.data.frame(LRPStatus)$LRPStatus
  
  status<-ifelse(statusEsts$rapidStatus == "Red", "Below", "Above")
  
  statusEsts<-statusEsts %>% add_column(LRPStatus=status)
  
  g<-ggplot(statusEsts) +
    geom_path(aes(y=Escp, x=yr), colour='black', alpha=0.5) +
    geom_point(aes(y=Gen_Mean, x=yr, colour=LRPStatus)) +
    geom_hline(aes(yintercept=Sgen), colour="orange") +
    ylab("Escapement") + xlab("Year") +
    scale_colour_manual(guide = NULL, breaks = c("Above", "Below"), values=c("gray40", "red")) +
    facet_wrap(~interaction(CU_Name), scales = "free_y") +
    theme_classic()
  

  # Save plot
  ggsave(paste(outDir,"/Figures/",plotName,".png",sep=""), plot = g,
         width = 10, height = 5, units = "in") 
  

  
}



# =================================================================================
# Functions for plotting retrospective analysis results
# ================================================================================

plotAnnualRetro<-function(Dat, Name ,outDir, useGenMean = F, genYrs = 3) {
  
  # Calculate aggregate escapement (note: CUs with no LBM but with escapement data will be included)
  Dat.AggEscp <-  Dat %>% group_by(yr) %>%  summarise(Agg_Escp = sum(Escp))
  
  if(useGenMean == T){
    Dat.AggEscp$Agg_Escp <- rollapply(Dat.AggEscp$Agg_Escp, genYrs, gm_mean, fill = NA, align="right")
    Ylab <- "Gen. Mean of Agg. Escapement"
  } else {
    Ylab <- "Agg. Escapement"
  }
  
  # Extract retrospective LRP estimates
  LRPs<-read.csv(paste(outDir,"/DataOut/AnnualRetrospective/", Name, "/annualRetro_LRPs.csv", sep=""))
  
  # Create plot
  annual_LRP_plot <- ggplot(data=LRPs, mapping = aes(x=retroYear, y=LRP)) +
    geom_line(aes(col = "LRP"), size = 2) +
    geom_line(data=Dat.AggEscp, mapping=aes(x=yr,y=Agg_Escp), size = 1.2) +
    xlab("Year") + ylab(Ylab) +
    scale_color_manual(values = 'steelblue2') +
    labs(color = '') +
    #coord_cartesian(ylim=c(0,2500000)) + # LW: optional for making all plots have same y limits, for comparisons.
    theme_classic() + 
    theme(legend.position = "none")

  #if(length(LRPs$LRP_lwr[is.na(LRPs$LRP_lwr)=="TRUE"])==0) { # LW: removed this if statement to plot confidence intervals, even though some will have NA years
    annual_LRP_plot <- annual_LRP_plot + 
      geom_ribbon(aes(ymin = LRP_lwr, ymax = LRP_upr, x=retroYear), fill = "steelblue2", alpha = 0.5)
  #}
  
  # Save plot
  ggsave(paste(outDir,"/Figures/AnnualRetrospective/", Name, "/Annual_LRP.pdf",sep=""), plot = annual_LRP_plot,
         width = 3.5, height = 3, units = "in")    
  
}



#=================================================================================
# Compare multiple retrospective LRPs from different probability thresholds on one plot
# --- separate plot for each method 
#=============================================================================


plotAnnualRetro_CompareProbs<-function(Dat, Names, L_Names = NA, pList, outDir, useGenMean, genYrs) {
  
  if(is.na(L_Names[1])){
    L_Names <- Names
  }
  
  # Calculate aggregate escapement (note: CUs with no LBM but with escapement data will be included)
  Dat.AggEscp <-  Dat %>% group_by(yr) %>%  summarise(Agg_Escp = sum(Escp))
  
  if(useGenMean == T){
    Dat.AggEscp$Agg_Escp <- rollapply(Dat.AggEscp$Agg_Escp, genYrs, gm_mean, fill = NA, align="right")
    Ylab <- "Abundance"
  } else {
    Ylab <- "Agg. Escapement"
  }
  
  # Extract retrospective LRP estimates
  for (i in 1:length(Names)) {
    for (pp in 1:length(pList)) {
      if (i ==1 & pp ==1) {
        LRPs <- read.csv(paste(outDir,"/DataOut/AnnualRetrospective/", Names[i],"_", pList[pp], "/annualRetro_LRPs.csv", sep=""))
        LRPs$Mod <- L_Names[1]
        LRPs$pVal<-pList[pp]
        } else {
        LRP_New <- read.csv(paste(outDir,"/DataOut/AnnualRetrospective/", Names[i],"_", pList[pp], "/annualRetro_LRPs.csv", sep=""))
        LRP_New$Mod <- L_Names[i]
        LRP_New$pVal<-pList[pp]
        LRPs <- rbind(LRPs, LRP_New)
      }
    } # end of p loop 
  } # end of Names loop
  
  LRPs$pVal<-as.factor(LRPs$pVal)

 # Embedded function to plot one model type:
    plotByModel<-function(i, Dat, L_Names) {
    LRPs<-Dat %>% filter(Mod == L_Names[i])
    # Create plot
    annual_LRP_plot <- ggplot() +
      geom_line(data=Dat.AggEscp, mapping=aes(x=yr,y=Agg_Escp), size = 1.2) +
      geom_ribbon(data = LRPs, aes(x=retroYear, y=LRP, colour = pVal, ymin = LRP_lwr, ymax = LRP_upr, fill = pVal), alpha = 0.5) +
      geom_line(data = LRPs,aes(x=retroYear, y=LRP, colour = pVal), size = 1.2) +
        xlab("Year") + ylab(Ylab) +
      theme_classic() +
      ggtitle(L_Names[i]) +
      theme(legend.position = "none")
    annual_LRP_plot
  }
  
  # Make multi-panel plot
    pLRPs<-lapply(1:length(L_Names), plotByModel, Dat = LRPs, L_Names)
    pMulti<-do.call(grid.arrange,  pLRPs)
    
  # Save multi=panel plot
    ggsave(paste(outDir,"/Figures/LRP_compareRetro.pdf",sep=""), plot = pMulti) 
  
}





#=================================================================================
# Compare multiple retrospective LRPs from different methods on one plot
# --- separate plot for each probability threshold
#=============================================================================


plotAnnualRetro_CompareMethods<-function(Dat, Names, L_Names = NA, pList, outDir, useGenMean, genYrs) {
  
  if(is.na(L_Names[1])){
    L_Names <- Names
  }
  
  # Calculate aggregate escapement (note: CUs with no LBM but with escapement data will be included)
  Dat.AggEscp <-  Dat %>% group_by(yr) %>%  summarise(Agg_Escp = sum(Escp))
  
  if(useGenMean == T){
    Dat.AggEscp$Agg_Escp <- rollapply(Dat.AggEscp$Agg_Escp, genYrs, gm_mean, fill = NA, align="right")
    Ylab <- "Abundance"
  } else {
    Ylab <- "Agg. Escapement"
  }
  
  # Extract retrospective LRP estimates
  for (i in 1:length(Names)) {
    for (pp in 1:length(pList)) {
      if (i ==1 & pp ==1) {
        LRPs <- read.csv(paste(outDir,"/DataOut/AnnualRetrospective/", Names[i],"_", pList[pp], "/annualRetro_LRPs.csv", sep=""))
        LRPs$Model <- L_Names[1]
        LRPs$pVal<-pList[pp]
      } else {
        LRP_New <- read.csv(paste(outDir,"/DataOut/AnnualRetrospective/", Names[i],"_", pList[pp], "/annualRetro_LRPs.csv", sep=""))
        LRP_New$Model <- L_Names[i]
        LRP_New$pVal<-pList[pp]
        LRPs <- rbind(LRPs, LRP_New)
      }
    } # end of p loop 
  } # end of Names loop
  
  LRPs$pVal<-as.factor(LRPs$pVal)
  
  
  # Embedded function to plot one model type:
  plotByProb<-function(i, Dat, pList) {
    LRPs<-Dat %>% filter(pVal == pList[i])
    # Create plot
    annual_LRP_plot <- ggplot() +
      geom_line(data=Dat.AggEscp, mapping=aes(x=yr,y=Agg_Escp), size = 1.2) +
      #geom_ribbon(data = LRPs, aes(x=retroYear, y=LRP, colour = Mod, ymin = LRP_lwr, ymax = LRP_upr, fill = Mod), alpha = 0.3) +
      geom_ribbon(data = LRPs, aes(x=retroYear, y=LRP, ymin = LRP_lwr, ymax = LRP_upr, fill = Model), alpha = 0.3) +
      geom_line(data = LRPs,aes(x=retroYear, y=LRP, colour = Model), size = 1.2) +
      xlab("Year") + ylab(Ylab) +
      theme_classic() +
      #ggtitle(paste("Probability = ", pList[i], "%", Sep="")) +
      theme(legend.position = c(0.2, 0.2))
   
    annual_LRP_plot
  }
  
  # Make multi-panel plot
 # pLRPs<-lapply(1:length(L_Names), plotByModel, Dat = LRPs, L_Names)
   pLRPs<-lapply(1:length(pList), plotByProb, Dat = LRPs, pList=pList)
   
   pMulti<-do.call(grid.arrange,  pLRPs)
  
  # Save multi=panel plot
  ggsave(paste(outDir,"/Figures/coho_LRP_compareRetro.png",sep=""), plot = pMulti, 
         units="in", width=5.5, height=4.4) 
  
}



#==================================================================
# Plot logistic model from 
#==============================================================

plotLogistic <- function(Data, Preds, LRP, useGenMean = F, plotName, outDir, p=0.95, useBern_Logistic = F){

  # y label different depending on model
  if(useBern_Logistic){
    Ylab = "Pr(All CUs > Lower Benchmark)"
  } else {
    Ylab = "Prop. CUs > Lower Benchmark"
  }
  
  if(useGenMean){
      Xlab = "Gen. Mean Agg. Escapement"
    }else {
      Xlab = "Aggregate Escapement"
  }
  
  Preds<-Preds[is.na(Preds$fit)==FALSE,]
  
  # Create plot
  annual_LRP_plot <- ggplot(data=Preds, mapping=aes(x=xx,y=fit)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, x=xx), fill = "grey85") +
    geom_line(mapping=aes(x=xx, y=fit), col="red", size=1, na.rm=TRUE) +
    geom_line(mapping=aes(x=xx, y=upr), col="grey85", na.rm=TRUE) +
    geom_line(mapping=aes(x=xx, y=lwr), col="grey85", na.rm=TRUE) +
    geom_point(data=Data, aes(x = xx, y = yy)) +
    #geom_segment(aes(y=0.01, yend=0.01, x=2200, xend=230873), colour="dodgerblue") + #LW: add in likelihood penalty SD
    #geom_point(aes(y=0.01, x=mean(c(2200, 230873))), colour="dodgerblue") + #LW: add in likelihood penalty mean
    xlab(Xlab) + ylab(Ylab) +
    coord_cartesian(ylim=c(0,1)) +
    #  ggtitle(paste("Retrospective Year", max(Dat.LRP$yr), sep=" ")) +
    theme_classic()
  
  
  
  # If LRP$lwr is <0, assume can't fit LRP
  if(LRP$lwr > 0 & is.na(LRP$lwr) == FALSE) {
    annual_LRP_plot <- annual_LRP_plot + 
      geom_vline(xintercept=LRP$fit, color="orange", size = 1) +
      geom_hline(yintercept= p, linetype="dotted", color="black", size = 0.5) + 
      geom_vline(xintercept = LRP$lwr, linetype = "dashed", color = "orange", size = 1) +
      geom_vline(xintercept = LRP$upr, linetype = "dashed", color = "orange", size = 1) 
  }
  
  # Save plot
  ggsave(paste(outDir, "/",plotName,".pdf",sep=""), plot = annual_LRP_plot,
         width = 4, height = 3, units = "in")      
  
}  



#===================================================================================
# Plot status by year 
#===================================================================================


plotStatusBarsCoho_byYear<-function(LRP_estYr, retroYears, Dir, genYrs,AggEscp,EscpDat,  modelFitList, 
                                    multiDimList, projLRPList = NULL, mutliDimList, ps_Prop,
                                    WSP_estYr=NULL, WSP_AboveLRP=NULL, outDir, fName) {

  
  Status_DF <- data.frame(LRP_estYr = numeric(), retroYear=numeric(), Name = character(), AboveLRP = character())
  
  # Extract / calculate data for plot ==========================================================

  # Step 1: Get output summaries for logistic-regression model LRPs"
    # # --- read in LRPs for each model option
    
  for (mm in 1:length(modelFitList)) {
    
    retroResults <- read.csv(paste(Dir,"/DataOut/AnnualRetrospective/", modelFitList[mm], "/annualRetro_LRPs.csv", sep=""))
    
    
    # Labels if only using IM models:
    if (retroResults$BMmodel[1] == "SR_IndivRicker_Surv") name1<-"Abund: Logistic:Sgen-Ricker"
    if (retroResults$BMmodel[1] == "SR_IndivRicker_SurvCap") name1<-"Abund: Logistic:Sgen-priorCap"
    if (retroResults$BMmodel[1] == "ThreshAbund_Subpop1000_ST") name1<-"Abund: Logistic:IFCRT"
    
    # If including probability in label:
    # name3<-strsplit(modelFitList[mm], "_")[[1]][2]
    # retroResults$Name <- paste(name1,name3,sep="_")

    # If no probability in label:
    retroResults$Name<-name1
    
      if(mm == 1){
         LRP_DF <- retroResults
      } else {
         New_Rows <- retroResults
         LRP_DF <- rbind(LRP_DF, New_Rows)
      }
    
    # Loop over years and summarize status as above (True) or below (False) the LRP for each year
    for (yy in 1:length(retroYears)) { 
      # Check if LRP is valid (Error bounds don't cross 0)
      LRPs <- retroResults %>% filter(retroYear == LRP_estYr) %>% select(LRP, LRP_lwr, LRP_upr)
      # Now compare to that year's Agg Escapement
      if(LRPs$LRP_lwr < 0 | LRPs$LRP == Inf) { 
        Status <- NA 
      } else {
        Curr_Esc <- AggEscp %>% filter(yr == retroYears[yy]) %>% pull(Gen_Mean)
        Status <-  Curr_Esc > LRPs$LRP
      }
      New_Row <- data.frame(LRP_estYr = LRP_estYr, retroYear = retroYears[yy], Name = retroResults$Name[1], AboveLRP=Status)
      Status_DF <- rbind(Status_DF, New_Row)
      
    } # end of year (yy) loop
    
    
  } # end of model (mm) loop

  
  # Step 2: Add projected LRPs (Optional)
  if (!is.null(projLRPList)) {
    
    projResults <- as_tibble(read.csv(paste(Dir,"/DataOut/ProjectedLRPs/projectedLRPs.csv", sep="")))
    
    for (j in 1:length(projLRPList)) {
      
      tmp<-strsplit(projLRPList[j],"_")
      
      if (tmp[[1]][2] == "priorCap") {
        OM.j<-paste(tmp[[1]][1],tmp[[1]][2], sep="_")
        probThresh.j<-as.numeric(tmp[[1]][3])/100
      } else {
        OM.j<-tmp[[1]][1]
        probThresh.j<-as.numeric(tmp[[1]][2])/100
      }
     
      dum<-projResults %>% filter(OM==OM.j,ProbThresh==probThresh.j)
      LRP<-dum$LRP
      
      # if (OM.j == "IM.Base") name1<-"AggAb_Proj: Sgen_IM"
      # if (OM.j == "HM.Base") name1<-"AggAb_Proj: Sgen_HM"
      # if (OM.j == "IMCap.Base") name1<-"AggAb_Proj: Sgen_IM.HiSrep"
      # if (OM.j == "HMCap.Base") name1<-"AggAb_Proj: Sgen_HM.HiSrep"
      
      
      # Labels if only using IM models:
      if (OM.j == "Ricker") name1<-"Abund: Proj:Sgen-Ricker"
      if (OM.j == "Ricker_priorCap") name1<-"Abund: Proj:Sgen-priorCap"
      
      # If includuding probability in label:
      # name2<-probThresh.j * 100
      # Name <- paste(name1,name2,sep=" ")
      
      # if not including probability:
      Name <- name1
      
      # Loop over years and summarize status as above (True) or below (False) the LRP for each year
      for (yy in 1:length(retroYears)) { 
    
        Curr_Esc <- AggEscp %>% filter(yr == retroYears[yy]) %>% pull(Gen_Mean)
        Status <-  Curr_Esc > LRP
        New_Row <- data.frame(LRP_estYr = LRP_estYr, retroYear = retroYears[yy], Name = Name, AboveLRP=Status)
        Status_DF <- rbind(Status_DF, New_Row)
               
      } # end of year (yy) loop
      
    } # end of projection model loop
    
  } # end of !is.null(projLRPList)
  

  # Step 3: Assess status for proportion-based LRP using Sgen (optional)
  
   # --- for each year, extract Sgens and calc proportion above Sgen

  # Make a list of only modelFits that are based on SR model fits to get Sgen:
  # -- For coho, I can do this by excluding the distributional models based on absolute abundance thresholds
  
  SRmodelList<-modelFitList[-grep("SPopAbundThresh",modelFitList)]
  
  for (mm in 1:length(SRmodelList)) {

    SRmodName<-strsplit(SRmodelList[mm],"_")[[1]][1]
    

     if(SRmodName =="Bern.IndivRickerSurv") {
       propName<-"Prop: Sgen-Ricker"
       fitName<-"AllEsts_Indiv_Ricker_Surv"
     }
     if(SRmodName =="Bern.IndivRickerSurvCap") {
       propName<-"Prop: Sgen-priorCap"
       fitName<-"AllEsts_Indiv_Ricker_Surv_priorCap"
     }
     
    CU_Params <- read.csv(paste(Dir, "/DataOut/ModelFits/",fitName,".csv",sep=""))
    CU_Params <- CU_Params %>% filter(Param == "Sgen")

    # -- add generational means to CU-level escapements; we will compare CU benchmarks to these
    EscpDat.mm <- EscpDat %>% group_by(CU) %>% mutate(Gen_Mean = rollapply(Escp, genYrs, gm_mean, fill = NA, align="right"))  %>%
      filter(is.na(Gen_Mean) == F)

    # -- join together Escp data with Sgens
    CU_Status <- left_join(EscpDat.mm[ , c("CU_Name", "yr", "Escp", "Gen_Mean")],
                         CU_Params[, c( "CU_Name", "Estimate")],
                         by = c("CU_Name" = "CU_Name"))
    colnames(CU_Status)[colnames(CU_Status)=="Estimate"] <- "Sgen"

    # --- for each year, get proportion of stocks above Sgen
    NCUs <- length(unique(CU_Status$CU_Name))
    CU_Status_Summ <- CU_Status %>%  group_by(yr) %>% filter(yr %in% retroYears) %>%
    summarise(Prop = sum(Gen_Mean > Sgen)/NCUs)

    # --- for each year, add to Status_DF using p thresholds, Ps
    #    --- note: should make these an input variable in the future
    Ps <- ps_Prop
    for(pp in 1:length(Ps)){
      # Would use this is wanted to show the proportion of CUs above Sgen
      #Name <- paste(propName,Ps[pp]*100)
      Name <- propName
      Status <- CU_Status_Summ$Prop >= Ps[pp]
      New_Rows <- data.frame(LRP_estYr, retroYear = CU_Status_Summ$yr, Name , AboveLRP=Status)
      Status_DF <- rbind(Status_DF, New_Rows)
    }
  }


  # Step 4: Add mutlidimensional results (optional)
 
  for (j in 1:length(multiDimList)) {
    
    tmp<-strsplit(multiDimList[j],"_")
    
    if (tmp[[1]][2] == "priorCap") {
      OM.j<-paste(tmp[[1]][1],tmp[[1]][2], sep="_")
      probThresh.j<-as.numeric(tmp[[1]][3])/100
    } else {
      OM.j<-tmp[[1]][1]
      probThresh.j<-as.numeric(tmp[[1]][2])/100
    }
    
  
    if (OM.j == "Ricker") name1<-"Prop: Multidim-Ricker"
    if (OM.j == "Ricker_priorCap") name1<-"Prop: Multidim-priorCap"
      
    multiDimResults<-as_tibble(read.csv(paste(Dir,"/DataOut/multiDimStatusEsts_",OM.j,".csv", sep="")))
    
    # Loop over years and summarize status as above (True) or below (False) the LRP for each year
    for (yy in 1:length(retroYears)) { 
      
      multiDimDat<-multiDimResults %>% filter(yr == retroYears[yy])
      
      nRed<-nrow(multiDimDat[multiDimDat$rapidStatus=="Red",])
      Status<-ifelse(nRed>0,FALSE,TRUE)
      New_Row <- data.frame(LRP_estYr = LRP_estYr, retroYear = retroYears[yy], Name = name1, AboveLRP=Status)
      Status_DF <- rbind(Status_DF, New_Row)
      
    } # end of year (yy) loop
    
    
  } # end of multidimMethods
  
  # Step 5: Add row to Status_DF for 2014 status assessment (Optional)
  if (!is.null(WSP_estYr)) {
    New_Row <- data.frame(LRP_estYr,retroYear = WSP_estYr, Name = "Prop: WSP-2014", AboveLRP = WSP_AboveLRP)
    Status_DF <- rbind(Status_DF, New_Row)
    Status_DF <- arrange(Status_DF, Name)
  }

  
  # Make Plot =============================================================
  
  methods <- unique(Status_DF$Name)
  
  # Hack to re-order methods for plotting ===============
  methods<-c("Prop: WSP-2014", "","Prop: Multidim-Ricker", "Prop: Sgen-Ricker", "Abund: Logistic:Sgen-Ricker", "Abund: Proj:Sgen-Ricker", "",
             "Prop: Multidim-priorCap", "Prop: Sgen-priorCap","Abund: Logistic:Sgen-priorCap", "Abund: Proj:Sgen-priorCap", "",
             "Abund: Logistic:IFCRT") 
             
  # --- set-up pdf to save to
  #pdf(paste(outDir,"/Figures/", fName, ".pdf", sep=""), width=8.5, height=6.5)
  png(paste(outDir,"/Figures/", fName, ".png", sep=""), width=700, height=580)
  par( oma=c(3,10,5,3), mar=c(3,3,3,3), lend=2, xpd=T)
  
  # ---- specify colouts 
  #  cols <- c( "#FF0900", "#08AB0B") # (red, green)
  cols <- c( "#FF0900", "grey80") # (red, grey)
  
  
  #---- set xlims for all sites
  Xlow <- min(retroYears)
  Xhigh <- max(retroYears)
  
  # --- standardize aggregate escapement values
  AggEscp <- AggEscp %>% filter(yr %in% retroYears)
  AggEscp$StdGen_Mean<- AggEscp$Gen_Mean/mean(AggEscp$Gen_Mean)
  #AggEscp$StdAgg_Escp <- AggEscp$Agg_Escp/mean(AggEscp$Agg_Escp, na.rm=T) # - do not need this at present bc only comparing to Gen_Mean
  
  # --- plot margin, also use for jittering
  low <- min(AggEscp$StdGen_Mean, na.rm=T) 
  high <- max(AggEscp$StdGen_Mean, na.rm=T)
  
  # --- create empty plotting region
  # plot(1, type="n", xlab="", ylab="", xlim=c(Xlow, Xhigh), ylim=c(low-(high-low)/8*length(methods), high), axes=F,
  #      main = paste("LRP Estimation Year = ",LRP_estYr, sep=""))
  plot(1, type="n", xlab="", ylab="", xlim=c(Xlow, Xhigh), ylim=c(low-(high-low)/8*length(methods), high), axes=F)
  
  # --- add generational mean escapement
  lines(x=AggEscp$yr, y=AggEscp$StdGen_Mean, lwd=1.5)
  
  # --- loop over methods and ....
  for(mm in 1:length(methods)){
    #subset data by method
    Mdat <- Status_DF %>% filter(Name == methods[mm] & is.na(AboveLRP)==F)
    
    if(dim(Mdat)[1] > 0){
      #set y location for that method
      #set increment
      inc <- (high-low)/8
      y <- low-inc*(mm)
      #do each plot segment by segment with appropriate colors
      for(k in 1:dim(Mdat)[1]){
        segments(x0=Mdat$retroYear[k]-0.5, x1=Mdat$retroYear[k]+0.5, y0=y, y1=y, col=cols[Mdat$AboveLRP[k]+1], lwd=4)
      }
      #label the method
      text(x=(Xlow-((Xhigh-Xlow)/2)), y=y, labels=methods[mm], 
           xpd=NA, pos=4, cex=1)
    }
  }
  
  
  # Add vertical lines to plots to help distinguish years
  addVertLines_minor<-function(x,low,high,n) {
    segments(x0=x, y0=low-(high-low)/8*n, x1=x, y1=high, lty=3,col="grey85")
  }
  
  addVertLines_major<-function(x,low,high,n) {
    segments(x0=x, y0=low-(high-low)/8*n, x1=x, y1=high, lty=2)
  }
  

  sapply(seq(1999.5,2018.5,by=1), addVertLines_minor,low=low,high=high,n=length(methods))
  sapply(seq(1999.5,2015.5,by=5), addVertLines_major,low=low,high=high,n=length(methods))
  
  #sapply(seq(1998,2018,by=1), addVertLines_minor,low=low,high=high,n=length(methods))
  #sapply(seq(2000,2015,by=5), addVertLines_major,low=low,high=high,n=length(methods))
 
  # add x-axis label
  axis(1)
  mtext(side = 1, "Year", outer = T, cex=1.2)
  text(x=(Xlow-((Xhigh-Xlow)/2)),y=0.6,labels=bquote(underline("Method")), xpd=NA, pos=4, cex=1.2)
  
  # legend( "topright", bty="n", lty=c(1,1,1,1), lwd=c(2,2,4,2) , 
  #         col=c("black", "grey", "#08AB0B", Tcols[1]),
  #         legend=c( "Gen. Mean Esc.", "Escapement", 
  #                   "Status", "Status 95% CI"),
  #         cex = 0.7) 
  
  dev.off()
  
}

#---------------------------------------------------------------------
# Make Status Bar plots for WCVI Chinook
#---------------------------------------------------------------------

plotStatusBarsChinook_byYear<-function(LRP_estYr, retroYears,  genYrs,
                                       AggEscp,EscpDat, pLRP, 
                                       ps_Prop, WSP_estYr=NULL, 
                                       WSP_AboveLRP=NULL, 
                                       outDir, fName) {
  
  Status_DF <- data.frame(LRP_estYr = numeric(), retroYear=numeric(), Name = character(), AboveLRP = character())
  
  # Extract / calculate data for plot ==========================================================
  
  
  # Step 1: Get projection-basedLRPs 

    
   projResults <- as_tibble( read.csv( 
    paste(outDir,
          "/DataOut/ProjectedLRPs/ProjectedLRPsbaseER_ALLp.csv", 
          sep= "" )))
    
  for (j in 1:length(pLRP)) {
    
    dum<-projResults %>% filter(ProbThresh==pLRP[j])
    LRP<-dum$LRP
      
    Name <- paste("ProjLPR_p",pLRP[j], sep="")
      
      
    # Loop over years and summarize status as above (True) or below (False) 
    # the LRP for each year. Note, the raw aggregate escapement is use instead
    # of geometric generational mean, as there are frequent NAs, so running
    # mean would be over 2,3, or 4 years depending on the calendar year. To
    # avoid this inconsistency, the raw aggregate abundances are used (as in 
    # WSP assessment of sBC CK)
    for (yy in 1:length(retroYears)) { 
      Curr_Esc <- AggEscp %>% filter(yr == retroYears[yy]) %>% pull(Agg_Escp)
      Status <-  Curr_Esc > LRP
      New_Row <- data.frame(LRP_estYr = LRP_estYr, retroYear = retroYears[yy], 
                            Name = Name, AboveLRP=Status)
      Status_DF <- rbind(Status_DF, New_Row)
      
    } # end of year (yy) loop
      
  } # end of projection model loop
    

  
  
  
  # Step 2: Assess status for data-based LRP options based on the observed 
  # proportion of CUs above LRP
  
  
  # Inlet-specific Sgen's are estimated in the Watershed-Area-Model with 
  # assumptions about productivity derived from life-stage model and expert 
  # opinion. See Repository "Watershed-Area-Model", WCVILRPs.r calculated using 
  # all data up to 2020
  
  Inlet_Names <- unique(EscpDat$Inlet_Name)
  CU_Params <- read.csv(paste (wcviCKDir, "/DataIn/wcviRPs_noEnh.csv", 
                               sep="")) %>% filter(Stock %in% Inlet_Names) %>% 
    rename("Inlet_Name"="Stock")
  
  # Use raw escapaement instead ofgenerational means, becuase of large number 
  # of gaps 
  
  # -- join together Escp data with Sgens
  Inlet_Status <- left_join(EscpDat[ , c("Inlet_Name", "yr", "Escp")],
                         CU_Params[, c( "Inlet_Name", "SGEN")],
                         by = c("Inlet_Name" = "Inlet_Name"))
    
  # --- for each year, get proportion of stocks above Sgen
  NInlets <- length(unique(Inlet_Status$Inlet_Name))
  Inlet_Status_Summ <- Inlet_Status %>%  group_by(yr) %>% filter(yr %in% retroYears) %>%
    summarise(Prop = sum(Escp > SGEN)/NInlets)
    
  
  
  # --- for each year, add to Status_DF using p thresholds, Ps
  #    --- note: should make these an input variable in the future
  Ps <- ps_Prop
  for(pp in 1:length(Ps)){
    # Would use this if wanted to show the proportion of CUs above Sgen
    #Name <- paste(propName,Ps[pp]*100)
    Status <- Inlet_Status_Summ$Prop >= Ps[pp]
    New_Rows <- data.frame(LRP_estYr, retroYear = Inlet_Status_Summ$yr,  Name = "Ppn1: Sgen", AboveLRP=Status)
    New_Rows2 <- data.frame(LRP_estYr, retroYear = Inlet_Status_Summ$yr,  Name = "Ppn2: Multi-D", AboveLRP=Status)
    Status_DF <- rbind(Status_DF, New_Rows, New_Rows2)
  }
 
  # Step 4: Add row to Status_DF for 2011 status assessment (Optional)
  if (!is.null(WSP_estYr)) {
    New_Row <- data.frame(LRP_estYr,retroYear = WSP_estYr, Name = "WSP (2014 only)", AboveLRP = WSP_AboveLRP)
    Status_DF <- rbind(Status_DF, New_Row)
    Status_DF <- arrange(Status_DF, Name)
  }
  
  
  
  # Make Plot =============================================================
  
  methods <- unique(Status_DF$Name)
  
  # # Hack to re-order methods for plotting ===============
  # methods<-c("AggAb_Hist: Sgen", "AggAb_Proj: Sgen", "Prop: Sgen",
  #            "AggAb_Hist: Sgen.HiSrep", "AggAb_Proj: Sgen.HiSrep", "Prop: Sgen.HiSrep",
  #            "AggAb_Hist: Dist")
  
  
  # --- set-up pdf to save to
  #pdf(paste(outDir,"/Figures/", fName, ".pdf", sep=""), width=8.5, height=6.5)
  png(paste(outDir,"/Figures/", fName, ".png", sep=""), width=700, height=580)
  par( oma=c(3,10,5,3), mar=c(3,3,3,3), lend=2, xpd=T)
  
  # ---- specify colouts 
  #  cols <- c( "#FF0900", "#08AB0B") # (red, green)
  cols <- c( "#FF0900", "grey80") # (red, grey)
  
  
  #---- set xlims for all sites
  Xlow <- min(retroYears)
  Xhigh <- max(retroYears)
  
  # --- standardize aggregate escapement values
  AggEscp <- AggEscp %>% filter(yr %in% retroYears)
  AggEscp$Std<- AggEscp$Agg_Escp/mean(AggEscp$Agg_Escp, na.rm=T)
  #AggEscp$StdAgg_Escp <- AggEscp$Agg_Escp/mean(AggEscp$Agg_Escp, na.rm=T) # - do not need this at present bc only comparing to Gen_Mean
  
  # --- plot margin, also use for jittering
  low <- min(AggEscp$Std, na.rm=T) 
  high <- max(AggEscp$Std, na.rm=T)
  
  # --- create empty plotting region
  plot(1, type="n", xlab="", ylab="", xlim=c(Xlow, Xhigh), ylim=c(low-(high-low)/8*length(methods), high), axes=F,
       main = paste("LRP Estimation Year = ",LRP_estYr, sep=""))
  
  # --- add generational mean escapement
  lines(x=AggEscp$yr, y=AggEscp$Std, lwd=1.5)
  
  # --- loop over methods and ....
  for(mm in 1:length(methods)){
    #subset data by method
    Mdat <- Status_DF %>% filter(Name == methods[mm] & is.na(AboveLRP)==F)
    
    if(dim(Mdat)[1] > 0){
      #set y location for that method
      #set increment
      inc <- (high-low)/8
      y <- low-inc*(mm)
      #do each plot segment by segment with appropriate colors
      for(k in 1:dim(Mdat)[1]){
        segments(x0=Mdat$retroYear[k]-0.5, x1=Mdat$retroYear[k]+0.5, y0=y, y1=y, 
                 col=cols[Mdat$AboveLRP[k]+1], lwd=4)
      }
      noAss <- !(retroYears %in%  Mdat$retroYear)
      for(p in 1:length(retroYears)){
        if(noAss[p]) points(x=retroYears[p], y=y, pch=4, col=grey(0.8), 
                            cex=0.5)
      }
      #label the method
      text(x=(Xlow-((Xhigh-Xlow)/2)), y=y, labels=methods[mm], 
           xpd=NA, pos=4, cex=0.8)
    }
  }
  
  
  # Add vertical lines to plots to help distinguish years
  addVertLines_minor<-function(x,low,high,n) {
    segments(x0=x, y0=low-(high-low)/8*n, x1=x, y1=high, lty=3,col="grey85")
  }
  
  addVertLines_major<-function(x,low,high,n) {
    segments(x0=x, y0=low-(high-low)/8*n, x1=x, y1=high, lty=2)
  }
  
  
  sapply(seq(1999.5,2018.5,by=1), addVertLines_minor,low=low,high=high,n=length(methods))
  sapply(seq(1999.5,2015.5,by=5), addVertLines_major,low=low,high=high,n=length(methods))
  
  #sapply(seq(1998,2018,by=1), addVertLines_minor,low=low,high=high,n=length(methods))
  #sapply(seq(2000,2015,by=5), addVertLines_major,low=low,high=high,n=length(methods))
  
  # add x-axis label
  axis(1)
  mtext(side = 1, "Year", outer = T)
  text(x=(Xlow-((Xhigh-Xlow)/2)),y=0.6,labels=bquote(underline("Method")), xpd=NA, pos=4)
  
  # legend( "topright", bty="n", lty=c(1,1,1,1), lwd=c(2,2,4,2) , 
  #         col=c("black", "grey", "#08AB0B", Tcols[1]),
  #         legend=c( "Gen. Mean Esc.", "Escapement", 
  #                   "Status", "Status 95% CI"),
  #         cex = 0.7) 
  
  dev.off()
  

}


#---------------------------------------------------------------------
# Make Status Bar plots for Inside South Coast Chum
#---------------------------------------------------------------------

  
plotStatusBarsChum_byYear<-function(Status_DF, AggEscp, fName) {
    LRP_estYr <- max(Status_DF$year)
    # Make Plot =============================================================
    
    methods <- sort(unique(Status_DF$data_name))

    # --- set-up pdf to save to
    png(paste("Figures/", fName, ".png", sep=""), width=8, height=6, units="in", res=300)
    par( oma=c(1,8,1,0), mar=c(3,3,3,3), lend=2, xpd=T)
    
    # ---- specify colouts 
    cols <- c( "#FF0900", "grey80") # (red, grey)
    
    #---- set xlims for all sites
    Xlow <- min(Status_DF$year)
    Xhigh <- max(Status_DF$year)
    
    # --- standardize aggregate escapement values
    AggEscp <- AggEscp %>% filter(yr %in% unique(Status_DF$year))
    AggEscp$StdGen_Mean<- AggEscp$Gen_Mean/mean(AggEscp$Gen_Mean, na.rm=TRUE)
    #AggEscp$StdAgg_Escp <- AggEscp$Agg_Escp/mean(AggEscp$Agg_Escp, na.rm=T) # - do not need this at present bc only comparing to Gen_Mean
    
    # --- plot margin, also use for jittering
    low <- min(AggEscp$StdGen_Mean, na.rm=TRUE) 
    high <- max(AggEscp$StdGen_Mean, na.rm=TRUE)
    
    # --- create empty plotting region
    plot(1, type="n", xlab="", ylab="", xlim=c(Xlow, Xhigh), ylim=c(low-(high-low)/8*length(methods), high), axes=FALSE,
         main = paste("LRP Estimation Year = ",LRP_estYr, sep=""))
    
    # --- add generational mean escapement
    lines(x=AggEscp$yr, y=AggEscp$StdGen_Mean, lwd=1.5)
    
    # --- loop over methods and ....
    for(mm in 1:length(methods)){
      #subset data by method
      Mdat <- Status_DF %>% filter(data_name == methods[mm] & is.na(AboveLRP)==FALSE)
      
      if(dim(Mdat)[1] > 0){
        #set y location for that method
        #set increment
        inc <- (high-low)/8
        y <- low-inc*(mm)
        #do each plot segment by segment with appropriate colors
        for(k in 1:dim(Mdat)[1]){
          segments(x0=Mdat$year[k]-0.5, x1=Mdat$year[k]+0.5, y0=y, y1=y, col=cols[Mdat$AboveLRP[k]+1], lwd=4)
        }
        #label the method
        text(x=(Xlow-((Xhigh-Xlow)/2)), y=y, labels=methods[mm], 
             xpd=NA, pos=4, cex=0.8)
      }
    }
    
    
    # Add vertical lines to plots to help distinguish years
    addVertLines_minor<-function(x,low,high,n) {
      segments(x0=x, y0=low-(high-low)/8*n, x1=x, y1=high, lty=3,col="grey85")
    }
    
    addVertLines_major<-function(x,low,high,n) {
      segments(x0=x, y0=low-(high-low)/8*n, x1=x, y1=high, lty=2)
    }
    
    
    # sapply(seq(1999.5,2018.5,by=1), addVertLines_minor,low=low,high=high,n=length(methods))
    # sapply(seq(1999.5,2015.5,by=5), addVertLines_major,low=low,high=high,n=length(methods))
    
    #sapply(seq(1998,2018,by=1), addVertLines_minor,low=low,high=high,n=length(methods))
    #sapply(seq(2000,2015,by=5), addVertLines_major,low=low,high=high,n=length(methods))
    
    # add x-axis label
    axis(1)
    mtext(side = 1, "Year", outer = TRUE)
    text(x=(Xlow-((Xhigh-Xlow)/2)),y=0.6,labels=bquote(underline("Method")), xpd=NA, pos=4)
    
    # legend( "topright", bty="n", lty=c(1,1,1,1), lwd=c(2,2,4,2) , 
    #         col=c("black", "grey", "#08AB0B", Tcols[1]),
    #         legend=c( "Gen. Mean Esc.", "Escapement", 
    #                   "Status", "Status 95% CI"),
    #         cex = 0.7) 
    
    dev.off()
    
  }





plotAggStatus_byNCUs <- function(yearList, nCUList, LRPmodel, BMmodel, p, Dir, inputPrefix, plotAveLine, ymax) {
  
  # Read in data for all NCU levels and combine
  for (cc in 1:length(nCUList)) {
    
    if (cc == 1) {
      nCUDat<-read.csv(paste(Dir,"/DataOut/nCUCombinations/",inputPrefix,"_",nCUList[cc],"CUs.csv",sep=""))
      } else {
      newDat<-read.csv(paste(Dir,"/DataOut/nCUCombinations/",inputPrefix,"_",nCUList[cc],"CUs.csv",sep=""))
      nCUDat<-rbind(nCUDat,newDat)
    }
 
  }
  
  makeYrPlot<-function (year, Dat, aveLine, ymax) {
    
    Dat<-Dat %>% filter(retroYear==year, BMmodel == BMmodel)

    maxNCU<-max(Dat$nCUs)
    AllCU<-Dat %>% filter(nCUs == maxNCU)
    AllCU_status_lwr<-AllCU$Gen_Mean / AllCU$LRP_lwr
    AllCU_status_upr<- AllCU$Gen_Mean /AllCU$LRP_upr
    
    # manually create a jitter for 3 and 4 CU cases
    nCUs.jit<-Dat$nCUs
    nCUs.jit[nCUs.jit == 4] <- seq(3.8,4.2,length.out=length(nCUs.jit[nCUs.jit == 4]))
    nCUs.jit[nCUs.jit == 3] <- seq(2.7,3.3,length.out=length(nCUs.jit[nCUs.jit == 3]))
    # append jittered nCUs to Dat
    Dat<-add_column(Dat,nCUs.jit)
    
    # flag any LRP status estimates with error estimates that didn't converge and/or make sense
    errorFlag<-rep(FALSE,length(nCUs.jit))
    errorFlag[Dat$LRP_lwr < 0] <- TRUE
    errorFlag[is.na(Dat$LRP_lwr) == TRUE] <- TRUE
    errorFlag[is.na(Dat$LRP) == TRUE] <- TRUE
    
    # calculate nCU-specific mean and SE
    Dat2 <- Dat[errorFlag == FALSE,]
    Mean<-Dat2 %>% group_by(nCUs) %>% summarize(meanStatus=gm_mean(status))
    #SE<- Dat2 %>% group_by(nCUs) %>% summarize(meanStatus=stdErr(status))
    
    Dat2$status_lwr<-Dat2$Gen_Mean/Dat2$LRP_upr
    Dat2$status_upr<-Dat2$Gen_Mean/Dat2$LRP_lwr
    
    ymaxFlag<-rep(FALSE,nrow(Dat2))
    ymaxFlag[Dat2$status_upr > ymax] <- TRUE
    Dat2$yTrunc<-ymaxFlag
    Dat2$status_upr[Dat2$status_upr > ymax] <- ymax
    
   
    g<-ggplot(Dat2, aes(x=nCUs.jit, y=status)) +
        scale_x_reverse(breaks=unique(Dat2$nCUs)) +
        geom_errorbar(aes(x=nCUs.jit, ymax = status_upr, ymin = status_lwr), width = 0, size = 0.8, colour="grey") +
        geom_point(size=2.5) +   
        labs(title=year, x = "Number of CUs", y = "Aggregate Status") +
        theme_classic() +
        theme(axis.text=element_text(size=16), axis.title=element_text(size=16), plot.title = element_text(hjust = 0.5, size=16, face="bold")) +
        ylim(0,ymax) + 
        # add additional points for cases that didn't converge
        geom_point(data = Dat[errorFlag==TRUE,],aes(x=nCUs.jit, y=rep(AllCU$status,length(Dat[errorFlag==TRUE,"nCUs.jit"]))),shape=4, size=2.5)
    
    
    # theme(axis.text=element_text(size=12),
    #       axis.title=element_text(size=14,face="bold"))
    
    if (length(Dat2$yTrunc[Dat2$yTrunc==TRUE]) > 0) {
        # add additional points to show when ylim truncated
        g<-g+geom_point(data=Dat2[Dat2$yTrunc==TRUE,], aes(x=nCUs.jit, y = rep(ymax,length(nCUs.jit))), shape=17, size=2)
    }
    
    if (aveLine == TRUE) {

      g <- g + geom_segment(data = Mean %>% filter(nCUs==3), aes(x=2.6, xend=3.4, y=AllCU$status,yend=AllCU$status),colour="red", linetype="dashed") +
        geom_segment(data = Mean %>% filter(nCUs==4), aes(x=4.4, xend=3.7, y=AllCU$status,yend=AllCU$status),colour="red", linetype="dashed")
     
    }
    
}
  
  LRP.mod<-LRPmodel
  
  ps <- lapply(yearList, makeYrPlot, Dat = nCUDat %>% filter(LRPmodel == LRP.mod), 
               aveLine=plotAveLine, ymax=ymax)

  outputDir<-paste(Dir,"/Figures/nCUCombinations",sep="")
  
  if (file.exists(outputDir) == FALSE){
    dir.create(outputDir)
  } 
  
 
  outName1<-str_split(inputPrefix, "[.]")[[1]][2]
  outName2<-paste(str_split(outName1, "_")[[1]][1],str_split(outName1, "_")[[1]][2],sep="-") 
  
   plotName <- paste("coho-StatusByNCUs-",outName2, sep="")
  
   png(paste(outputDir,"/", plotName, ".png", sep=""), width=580, height=700)
   
   do.call(grid.arrange,  ps)
  
   dev.off()
   
}





plotAggStatus_byNCUs_Compare<-function(estYear, nCUList, Names, labelNames, p, Dir, plotAveLine=TRUE)  {
  
  
  makeMethodPlot<-function (methodName, aveLine, year) {
    
    if (methodName == "Bern.IndivRickerSurv") labName<-"Sgen: HM"
    if (methodName == "Bern.HierRickerSurv") labName<-"Sgen: IM"
    if (methodName == "Bern.IndivRickerSurvCap") labName<-"Sgen: HM.HiSrep"
    if (methodName == "Bern.HierRickerSurvCap") labName<-"Sgen: IM.HiSrep"
    if (methodName == "Bern.SPopAbundThreshST") labName<-"Dist"
    
    # Read in data for all NCU levels and combine
    for (cc in 1:length(nCUList)) {
      
      if (cc == 1) {
        nCUDat<-read.csv(paste(Dir,"/DataOut/nCUCombinations/",methodName,"_",p*100,"_",nCUList[cc],"CUs.csv",sep=""))
      } else {
        newDat<-read.csv(paste(Dir,"/DataOut/nCUCombinations/",methodName,"_",p*100,"_",nCUList[cc],"CUs.csv",sep=""))
        nCUDat<-rbind(nCUDat,newDat)
      }
      
    }
    
    Dat<-nCUDat %>% filter(retroYear==year)
    maxNCU<-max(Dat$nCUs)
    AllCU<-Dat %>% filter(nCUs == maxNCU)
    AllCU_status_lwr<-AllCU$Gen_Mean / AllCU$LRP_lwr
    AllCU_status_upr<- AllCU$Gen_Mean /AllCU$LRP_upr
    
    # manually create a jitter for 3 and 4 CU cases
    nCUs.jit<-Dat$nCUs
    nCUs.jit[nCUs.jit == 4] <- seq(3.8,4.2,length.out=length(nCUs.jit[nCUs.jit == 4]))
    nCUs.jit[nCUs.jit == 3] <- seq(2.7,3.3,length.out=length(nCUs.jit[nCUs.jit == 3]))
    # append jittered nCUs to Dat
    Dat<-add_column(Dat,nCUs.jit)
    
    # flag any LRP status estimates with error estimates that didn't converge and/or make sense
    errorFlag<-rep(FALSE,length(nCUs.jit))
    errorFlag[Dat$LRP_lwr < 0] <- TRUE
    errorFlag[is.na(Dat$LRP_lwr) == TRUE] <- TRUE
    
    # calculate nCU-specific mean and SE
    Dat2 <- Dat[errorFlag == FALSE,]
    Mean<-Dat2 %>% group_by(nCUs) %>% summarize(meanStatus=gm_mean(status))
    SE<- Dat2 %>% group_by(nCUs) %>% summarize(meanStatus=stdErr(status))
    
    g<-ggplot(Dat, aes(x=nCUs.jit, y=status)) +
      scale_x_reverse(breaks=unique(Dat$nCUs)) +
      geom_errorbar(aes(x=nCUs.jit, ymax = Gen_Mean/LRP_lwr, ymin = Gen_Mean/LRP_upr), width = 0, colour="grey") +
      geom_point() +   
      labs(title=labName, x = "Number of CUs", y = "Aggregate Status") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0,4.5) + 
      geom_point(data = Dat[errorFlag==TRUE,],aes(x=nCUs.jit, y=status),shape=4, size=3) # add additional points for cases that didn't converge
    
    if (aveLine == TRUE & nrow(Mean)>1) {
      
      g <- g + geom_segment(data = Mean %>% filter(nCUs==3), aes(x=2.6, xend=3.4, y=meanStatus,yend=meanStatus),colour="red", linetype="dashed") +
        geom_segment(data = Mean %>% filter(nCUs==4), aes(x=4.4, xend=3.7, y=meanStatus,yend=meanStatus),colour="red", linetype="dashed")
      
    }
    
    print(g)
    
  }
  
  ps <- lapply(Names, makeMethodPlot, aveLine=plotAveLine, year=estYear)
  
  pdf(paste(Dir,"/Figures/compare_StatusByNCUs.pdf", sep=""), width=6.5, height=6.5)
  do.call(grid.arrange,  ps)
  dev.off()
  
}




plotLRP.CV_by_nCUs<-function(yearList, nCUList, LRPmodel, BMmodel, p, Dir, inputPrefix, fName) {
  
  # Read in data for all NCU levels and combine
  for (cc in 1:length(nCUList)) {
    
    if (cc == 1) {
      nCUDat<-read.csv(paste(Dir,"/DataOut/nCUCombinations/",inputPrefix,"_",p*100,"_",nCUList[cc],"CUs.csv",sep=""))
    } else {
      newDat<-read.csv(paste(Dir,"/DataOut/nCUCombinations/",inputPrefix,"_",p*100,"_",nCUList[cc],"CUs.csv",sep=""))
      nCUDat<-rbind(nCUDat,newDat)
    }
    
  }
  
  makeYrPlot<-function (year, Dat) {
    
    Dat<-Dat %>% filter(retroYear==year)
    
    maxNCU<-max(Dat$nCUs)
    #AllCU<-Dat %>% filter(nCUs == maxNCU)
    #AllCU_status_lwr<-AllCU$Gen_Mean / AllCU$LRP_lwr
    #AllCU_status_upr<- AllCU$Gen_Mean /AllCU$LRP_upr
    
    # manually create a jitter for 3 and 4 CU cases
    nCUs.jit<-Dat$nCUs
    nCUs.jit[nCUs.jit == 4] <- seq(3.8,4.2,length.out=length(nCUs.jit[nCUs.jit == 4]))
    nCUs.jit[nCUs.jit == 3] <- seq(2.7,3.3,length.out=length(nCUs.jit[nCUs.jit == 3]))
    # append jittered nCUs to Dat
    Dat<-add_column(Dat,nCUs.jit)
    
    # flag any LRP status estimates with error estimates that didn't converge and/or make sense. and then remove
    errorFlag<-rep(FALSE,length(nCUs.jit))
    errorFlag[Dat$LRP_lwr < 0] <- TRUE
    errorFlag[is.na(Dat$LRP_lwr) == TRUE] <- TRUE
    Dat <- Dat[errorFlag == FALSE,]

    g<-ggplot(Dat, aes(x=nCUs.jit, y=cv)) +
      scale_x_reverse(breaks=unique(nCUs)) +
      #geom_errorbar(aes(x=nCUs.jit, ymax = Gen_Mean/LRP_lwr, ymin = Gen_Mean/LRP_upr), width = 0, colour="grey") +
      geom_point() +   
      labs(title=year, x = "Number of CUs", y = "CV on LRP estimate") +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0,0.5)
  }
  
  LRP.mod<-LRPmodel
  BM.mod<-BMmodel
  
  ps <- lapply(yearList, makeYrPlot, Dat = nCUDat %>% filter(LRPmodel == LRP.mod, BMmodel == BM.mod) )
  
  if (LRPmodel == "BernLogistic") plotName <- paste("LRP_CVs_ByNCUs_Bern.",BM.mod,p*100, sep="")
  if (LRPmodel == "BinLogistic") plotName <- paste("LRP_CVs_ByNCUs_Bin.",BM.mod,p*100, sep="")
  png(paste(Dir,"/Figures/nCUCombinations/", plotName, ".png", sep=""))
  
  do.call(grid.arrange,  ps)
  
  dev.off()
  
  
}



plotLRPCompare<-function(LRP_estYr, modelFitList,  pList, outDir, fName) {
  
  # Step 1: Get output summaries for modelled LRPs (i.e., integrated Bern and Bin regressions):
  
  # # --- read in LRPs for each model option
  
  for (mm in 1:length(modelFitList)) {
    
    for (pp in 1:length(pList)) {
    
      retroResults <- read.csv(paste(outDir,"/DataOut/AnnualRetrospective/", modelFitList[mm], "_",pList[pp],"/annualRetro_LRPs.csv", sep=""))
  
      # Set-up name for labelling plot ====================
      if (retroResults$BMmodel[1] == "SR_HierRicker_Surv") name1<-"HM"
      if (retroResults$BMmodel[1] == "SR_IndivRicker_Surv") name1<-"IM"
      if (retroResults$BMmodel[1] == "SR_HierRicker_SurvCap") name1<-"HM.HiSrep"
      if (retroResults$BMmodel[1] == "SR_IndivRicker_SurvCap") name1<-"IM.HiSrep"
      if (retroResults$BMmodel[1] == "ThreshAbund_Subpop1000_ST") name1<-"Dist"
      
      
      retroResults$Name <- name1
    
      LRPs <- retroResults %>% filter(retroYear == LRP_estYr) %>% select(Name, LRP, LRP_lwr, LRP_upr)
      LRPs$p<-pList[pp]
      
      if(mm == 1 & pp == 1){
       LRP_DF <- LRPs
      } else {
       New_Rows <- LRPs
       LRP_DF <- rbind(LRP_DF, New_Rows)
      }
    }
  }

  # assign plot order:
  LRP_DF$Name<-factor(LRP_DF$Name, levels = c("IM", "HM", "IM.HiSrep", "HM.HiSrep", "Dist"))
  LRP_DF$p<-factor(LRP_DF$p, levels=pList)
  
 pdf(paste(outDir,"/Figures/", fName,".pdf", sep=""),width=5, height=4)
  dodge <- position_dodge(.4)
  p.lrp<- ggplot(data=LRP_DF, mapping=aes(x=Name, y=LRP)) +
    geom_errorbar(aes(x=Name,ymax=LRP_upr,ymin=LRP_lwr,col=p), width=0, position=dodge) +
    geom_point(mapping=aes(x=Name, y=LRP, col=p), size=2, position=dodge) +
    xlab("") + ylab("Aggregate LRP") +
    ylim(min(LRP_DF$LRP_lwr), max(LRP_DF$LRP_upr)) +
    theme_classic() +
    scale_color_grey()
   # scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  print(p.lrp)
  dev.off()
  
  
}



plotProjected<-function(Data, LRP, plotName, outDir, p) {
 
  g <- ggplot(data=Data, mapping=aes(x=as.factor(ppnCUs),y=AggSpawners, fill = as.factor(ppnCUs))) + geom_violin()
  
  # g <- g + stat_summary(
  #  fun.min = function(z) { quantile(z,0.05) },
  #  fun.max = function(z) { quantile(z,0.95) },
  #  fun = median, col="navy")
  
  p.ind<-which(sort(unique(Data$ppnCUs)) == p)
  
  g <- g  + theme_light() + scale_fill_manual(values=rep("grey80", length(unique(Data$ppnCUs)))) +
    labs(y="Aggregate Escapement",
         x="Prop. CUs > Lower Benchmark") + theme(legend.position = "none") +
          geom_hline(yintercept=LRP$fit, linetype="solid", color="orange", size = 0.5) +
          geom_vline(xintercept= p.ind, linetype="dashed", color="black", size = 0.3) +
          geom_hline(yintercept = LRP$lwr, linetype = "dotted", color = "orange", size = 0.5) +
          geom_hline(yintercept = LRP$upr, linetype = "dotted", color = "orange", size = 0.5) +
          ylim(0,100000)
  
  g <- g + coord_flip()
  
  # Save plot
  ggsave(paste(outDir, "/",plotName,".pdf",sep=""), plot = g,
         width = 4, height = 3, units = "in")      
  
}



plotPostHist<-function(x, post, parName, CUNames) {
  
 margPost<-post %>% filter(stk==x) %>% select(parName)
  
 round.int<-2
 
  if (parName == "adjProd") {
    margPost <- margPost$adjProd
    xlab<-"Adjusted Productivity"
  }
 
  if (parName == "Sgen") {
    margPost <- margPost$Sgen
    xlab<-"Sgen"
    round.int<-0
  }
  if (parName == "alpha") {
    margPost <- margPost$alpha
    xlab<-"Alpha"
  }
  if (parName == "sigma") {
    margPost <- margPost$sigma
    xlab<-"Sigma"
  }
  if (parName == "beta") {
    margPost <- margPost$beta
    round.int<-5
    xlab<-"Beta"
  }
  if (parName == "cap") {
    margPost <- margPost$cap
    xlab<-"SRep"
  }

  margPost<-as.data.frame(margPost)

  x.text<-quantile(margPost$margPost, probs=0.85)[[1]]
  y.text<-max(hist(margPost$margPost,plot=F)$counts)*0.80
  
  p<-ggplot(margPost, aes(x=margPost)) + geom_histogram() +
    labs(title=CUNames[x], x=xlab, y = "Count") +
    geom_vline(aes(xintercept = mean(margPost)),col='red',size=1) +
    geom_text(aes(label=round(mean(margPost),round.int),y=0,x=mean(margPost)*1.25),
              vjust=-1,col='red',size=5) +
    geom_vline(aes(xintercept = quantile(margPost,0.95)[[1]]),col='blue',size=1) +
    geom_text(aes(label=round(quantile(margPost,0.95)[[1]],round.int),y=200,x=quantile(margPost,0.95)[[1]]*1.25),
              vjust=-1,col='blue',size=5) +
    geom_vline(aes(xintercept = quantile(margPost,0.05)[[1]]),col='blue',size=1) +
    geom_text(aes(label=round(quantile(margPost,0.05)[[1]],round.int),y=200,
                  x=quantile(margPost,0.05)[[1]]*1.25), vjust=-1,col='blue',size=5)
  
    #annotate("text", x=x.text, y=13000, label= "boat")
  
}



plotAgeProp_byCU<-function(CUages, outDir, plotName) {

  
  
  CUages$CU_Names[CUages$CU_Names=="Fraser_Canyon"]<-"Fraser Canyon"
  CUages$CU_Names[CUages$CU_Names=="Lower_Thompson"]<-"Lower Thompson"
  CUages$CU_Names[CUages$CU_Names=="Middle_Fraser"]<-"Middle Fraser"
  CUages$CU_Names[CUages$CU_Names=="North_Thompson"]<-"North Thompson"
  CUages$CU_Names[CUages$CU_Names=="South_Thompson"]<-"South Thompson"

  
  
  #cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  colList<-c("#999999","#E69F00", "#56B4E9", "#009E73", "#D55E00")
  
  
  g<-ggplot(CUages,aes(x=Year, y=age3, col=CU_Names)) + geom_line(size=1.1) +
    theme_classic() +
    xlab("Year") + ylab("Proportion Age 3") + labs(color='Conservation Unit')  +
    theme(legend.position = c(0.2, 0.25)) +
    scale_color_manual(values=colList) +
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=12,face="bold"))
  
  ggsave(paste(outDir,"/",plotName,".png",sep=""), plot = g,
         width = 5.25, height = 4.25, units = "in") 
  
  
}