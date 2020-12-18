# ==================================================================

# Functions that loop over different subsets of data (either via annual retrospective analyses or bootstrap of years or CUs).
### For each subset of data tested, functions from LRPFunctions.r will be called and results will be saved.
#
# runAnnualRetro    -  Runs retrospective analysis using specified model type. Calls Run_LRP or Run_Ricker_LRP functions from LRPFunctions.r (latter if integratedModel=T) 
# runNCUsRetro      -  Runs retrospective analysis using specified model type. Calls runAnnualRetro function
#
# ==================================================================


runAnnualRetro<-function(EscpDat, SRDat, startYr, endYr, BroodYrLag, genYrs, p = 0.95,
                         BMmodel, LRPmodel=NULL, integratedModel=F, useGenMean=T, 
                         TMB_Inputs=NULL, outDir, RunName, bootstrapMode=F, plotLRP=T) {
  
  yearList<-startYr:endYr
  
  # Put together Run info
  RunInfo <- data.frame(BMmodel, LRPmodel, useGenMean, integratedModel)
  
  # Loop over years ===============================================
  for (yy in 1:length(yearList)) {
  
      # Only use SR data for brood years that have recruited by current year yy 
       # (note: most recent brood year is calculated by subtracting BroodYearLag (e.g. 2 years) from current year)
      if (!is.null(SRDat)) {
        Dat <- SRDat %>%  filter(BroodYear <= (yearList[yy])-BroodYrLag)
      }
      EscpDat.yy <- EscpDat %>% filter(yr <= yearList[yy]) 
  
      # get trigger ready so can be used with both basic and TMB
      if (LRPmodel == "BernLogistic") useBern_Logistic <- TRUE
      if (LRPmodel == "BinLogistic") useBern_Logistic <- FALSE
      
      # Model Not Integrated (i.e., seperate benchmark and LRP functions are being used) ===========================
      if(integratedModel == F){
        
        # 1) Call specified benchmark function:

        
        # Case 1: LRP model is based on projections
        if (LRPmodel == "Proj_BinLogistic" | LRPmodel == "Proj_BernLogistic" ) {
          
          if (LRPmodel == "Proj_BernLogistic") useBern_Logistic <- TRUE
          if (LRPmodel == "Proj_BinLogistic") useBern_Logistic <- FALSE
          
          Dat$yr_num <- group_by(Dat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
          Dat$CU_ID <- group_by(Dat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
          EscDat <- EscpDat.yy %>%  right_join(unique(Dat[,c("CU_ID", "CU_Name")]))
          
          LRP_Mod<-Run_ProjRicker_LRP(SRDat = Dat, EscDat = EscDat, BMmodel = BMmodel, LRPmodel = LRPmodel, 
                                      useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs, outDir=outDir)
          
        }
        
      
        # Case 2: BM Model == 1000 fish threshold at sub-population level
        
        if (BMmodel == "ThreshAbund_Subpop1000_ST" | BMmodel == "ThreshAbund_Subpop1000_LT") {
        
          # Calcualte geomertic means for subpopulation escapements (where, subpopulation is sum of tributaries) 
          EscpDat.cu <- EscpDat.yy %>% group_by(CU_Name, CU_ID, Subpop_Name, yr)  %>% summarise(SubpopEscp=sum(Escp)) %>%
            mutate(Gen_Mean = rollapply(SubpopEscp, genYrs, gm_mean, fill = NA, align="right"))
          # Add a column indicating whether geometric mean escapement iss > 1000 fish threshold for each subpopulation 
          Above1000<-ifelse(EscpDat.cu$Gen_Mean >= 1000, 1, 0)
          EscpDat.cu<-EscpDat.cu %>% add_column(Above1000)
          # Calculate the number of subpopulations that are above 1000 fish threshold in each CU
          LBM_status_byCU <- EscpDat.cu %>% group_by(CU_Name, CU_ID, yr) %>% 
            summarise(Escp=sum(SubpopEscp), N = length(unique(Subpop_Name)),N_grThresh=sum(Above1000))
          # Add a column indicating whether >= 50% of subpopulations were above 1000 fish threshold in each CU (1 = yes, 0 = no) (short-term recovery goal)
          HalfGrThresh<-ifelse(LBM_status_byCU$N_grThresh >=  ceiling(LBM_status_byCU$N/2),1,0)
          # Add a column indicating wheter all subpopulations were above 1000 fish threshold in each CU (long-term recovery goal) 
          AllGrThresh<-ifelse(LBM_status_byCU$N_grThresh == LBM_status_byCU$N,1,0)
          LBM_status_byCU <- LBM_status_byCU %>% add_column(HalfGrThresh) %>% add_column(AllGrThresh)
          
          # Added dum 2 as a test to confirm our results match M. Arbeider code when we take geometric mean at the sub-population level before aggregating escapements
          # --- June 5, 2020: Have confirmed that our results do match . Our results for long-term recovery target are the same as those 
          #         from the "IFC Conservation Goal Targets.Rmd" file Michael provided us.
          #dum<-EscpDat.yy %>% group_by(Subpop_Name,yr) %>% summarise(Escp.sp = sum(Escp)) %>% mutate(Gen_Mean = rollapply(Escp.sp, genYrs, gm_mean, fill = NA, align="right"))
          #dum2<- dum %>% group_by(yr) %>% summarise(AggEscp.gm=sum(Gen_Mean))
          


      # 2) Call specified LRP function:
        
        LRP_Mod<-Run_LRP(EscDat=LBM_status_byCU, Mod = BMmodel, useBern_Logistic = useBern_Logistic, 
                         useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs)
        
        # This version with dum2=dum2 as an input was created to test consistency with M. Arbeider code; 
        # ---- if wanting to run with geometric means calculated at subpopulation lvel, will beed to add AggEscp_gmBySP arguement to LRP_Mod function
        #LRP_Mod<-Run_LRP(EscDat=LBM_status_byCU,Mod = BMmodel, useBern_Logistic = useBern_Logistic, 
        #                 useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs, dum2=dum2)
        
        
        }
        
        
        
        
    # Integrated model (i.e., benchmark and LRP are estimated simultaneously in TMB) ===================================
    
    } else if(integratedModel == T){
      # Prep data frame to work with function
      
      # Dat$yr_num <- group_indices(Dat, BroodYear) -1 
      # LW - changed above line to below. Reason: ... argument of group_keys() is deprecated
      Dat$yr_num <- group_by(Dat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
      
      # Dat$CU_ID <- group_indices(Dat, CU_ID) - 1
      # LW - changed above line to below. Reason: ... argument of group_keys() is deprecated
      Dat$CU_ID <- group_by(Dat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
      
      EscDat <- EscpDat.yy %>%  right_join(unique(Dat[,c("CU_ID", "CU_Name")]))
      
      LRP_Mod <- Run_Ricker_LRP(SRDat = Dat, EscDat = EscDat, BMmodel = BMmodel, Bern_Logistic = useBern_Logistic, 
                     useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs)
      
      # Compile Ricker parameters from TMB outputs
      RickerEsts <- LRP_Mod$All_Ests %>% filter(Param %in% c("A", "B" , "Smsy", "Sgen", "sigma"))
      output <- RickerEsts %>% dplyr::select(-Std..Error) %>% 
                               pivot_wider( names_from = "Param", values_from = "Estimate", names_prefix = "est_")
      Sgen_SE <- RickerEsts %>% filter(Param == "Sgen") %>% pull(Std..Error)
      Sgens <- output$est_Sgen
      output$low_Sgen <- Sgens - 1.96*Sgen_SE
      output$up_Sgen <- Sgens + 1.96*Sgen_SE
      output$retroYr <- yearList[yy]
      output <- cbind(output, RunInfo)  
      
      # Put SR data into dataframes
      if (yy == 1) {
        outSR.df<-output
      } else {
        outSR.df<-full_join(outSR.df, output)
      }
      
    }
    
    newLRP.df <- data.frame(retroYear = yearList[yy], LRP = LRP_Mod$LRP$fit, 
                                        LRP_lwr = LRP_Mod$LRP$lwr, LRP_upr = LRP_Mod$LRP$upr)
    newLRP.df <- cbind(newLRP.df, RunInfo)
    
    if (yy == 1) {
      outLRP.df<-newLRP.df
    } else {
      outLRP.df <- full_join(outLRP.df, newLRP.df)
    }
    
    # If plotLRP=T, plot model fit and data  
    if (plotLRP == T) {
      
      figDir <- paste(outDir, "Figures/AnnualRetrospective", sep="/")
      if (file.exists(figDir) == FALSE){
        dir.create(figDir)
      } 
      
      figDir <- paste(figDir, RunName, sep="/")
      if (file.exists(figDir) == FALSE){
        dir.create(figDir)
      }
      
      plotLogistic(Data = LRP_Mod$Logistic_Data, Preds = LRP_Mod$Preds, 
                   LRP = LRP_Mod$LRP, useGenMean = useGenMean,
                   plotName = paste("LogisticMod", yearList[yy], sep ="_"), outDir = figDir,
                   p = p, useBern_Logistic = useBern_Logistic)
    }
    
    

  } # end of year loop


 
  # Save outputs and make plots if function is not being run as part of a bootstrap
  if (bootstrapMode == F) {
  
    #Create directories to store outputs in
    retroDir<-paste(outDir, "DataOut/AnnualRetrospective", sep="/")
        if (file.exists(retroDir) == FALSE){
      dir.create(retroDir)
    }
    
    outputDir <- paste(retroDir, RunName, sep="/")
    if (file.exists(outputDir) == FALSE){
      dir.create(outputDir)
    } 
    
    figDir <- paste(outDir, "Figures/AnnualRetrospective", sep="/")
    if (file.exists(figDir) == FALSE){
      dir.create(figDir)
    } 
    
    figDir <- paste(figDir, RunName, sep="/")
    if (file.exists(figDir) == FALSE){
      dir.create(figDir)
    }
   
    
    
    # Save annual retrospective outputs as csv files =========================
    if (integratedModel == T) write.csv(outSR.df, paste(outputDir,"/annualRetro__SRparsByCU.csv", sep=""))
    write.csv(outLRP.df, paste(outputDir,"/annualRetro_LRPs.csv", sep=""))
  
    # in final year also plot geometric mean of aggregate abundance + LRPs
      # (Note: set her to always plot reference points relative to geometric mean)
    plotAnnualRetro(Dat = EscpDat.yy, Name  = RunName ,outDir = outDir, useGenMean = T, genYrs = genYrs)
  
    
  }
    
  out <- list()
  if (integratedModel == T) out$SRparsByCU <- outSR.df
  out$LRPs <- outLRP.df

  out
  
} # end of function 






runNCUsRetro <- function(nCUList, EscpDat, SRDat, startYr, endYr, BroodYrLag, genYrs, p,
                                BMmodel, LRPmodel, integratedModel,useGenMean, TMB_Inputs, outDir, RunName) {
 
  CU_Names<-unique(EscpDat$CU_Name)
  
  for (nn in 1:length(nCUList)) {
    
    nCUs<-nCUList[nn]
    CU_combn<-combn(CU_Names,nCUs)
    
    for (ii in 1:ncol(CU_combn)) {
      
      # Select CUs for this replicate
      CUs.ii<-CU_combn[,ii]
      
      # Filter data to only include sampled CUs
      EscpDat.ii <- EscpDat %>% filter(CU_Name %in% CUs.ii)
      if (!is.null(SRDat)) {
        SRDat.ii <- SRDat %>% filter(CU_Name %in% CUs.ii)
      } else {SRDat.ii <- SRDat}
      
      # Calculate aggregate escapement (with geometric mean added) for sampled CUs
      AggEscp.ii <- EscpDat.ii %>% group_by(yr) %>% summarise(Agg_Escp = sum(Escp)) %>% 
        mutate(Gen_Mean = rollapply(Agg_Escp, 3, gm_mean, fill = NA, align="right"))
      
      TMB_Inputs.ii<-TMB_Inputs
      
      # Specify priors on Srep specific to nCUs, if model with SRep prior is used:
      if (is.null(TMB_Inputs$cap_mean) == FALSE) {
        TMB_Inputs.ii$cap_mean<-TMB_Inputs$cap_mean[which(CUs.ii %in% as.character(unique(EscpDat$CU_Name)))]
      }
      
      # Run annual retrospective analyses for sampled CUs:
        out<- runAnnualRetro(EscpDat=EscpDat.ii, SRDat=SRDat.ii, startYr=startYr, endYr=endYr, BroodYrLag=BroodYrLag, 
                            genYrs=genYrs, p = p, BMmodel = BMmodel, LRPmodel=LRPmodel, integratedModel=integratedModel,
                            useGenMean=useGenMean, TMB_Inputs=TMB_Inputs.ii, outDir=cohoDir, RunName = NA, 
                            bootstrapMode = T, plotLRP=F)
          

          # Add output to table
          # - if first retrospective run, create dataframe to store outputs
          if (ii == 1) {
            output.df<-left_join(out$LRPs,AggEscp.ii, by = c("retroYear" = "yr"))
            nCUs<-rep(nCUList[nn],nrow(output.df))
            iCombn<-rep(ii,nrow(output.df))
            status<-output.df$Gen_Mean/output.df$LRP
            se<-(output.df$LRP_upr-output.df$LRP) / 1.96 # back-calculate what SE on LRP estimate was
            cv<-se/output.df$LRP # calculate CV as SE / LRP
            output.df<-data.frame(iCombn, nCUs, output.df, status, se, cv) 
           } else {
            new.df<-left_join(out$LRPs,AggEscp.ii, by = c("retroYear" = "yr"))
            nCUs<-rep(nCUList[nn],nrow(new.df))
            iCombn<-rep(ii,nrow(new.df))
            status<-new.df$Gen_Mean/new.df$LRP
            se<-(new.df$LRP_upr-new.df$LRP) / 1.96 # back-calculate what SE on LRP estimate was
            cv<-se/new.df$LRP # calculate CV as SE / LRP
            new.df<-data.frame(iCombn, nCUs, new.df, status, se, cv)
            output.df<-rbind(output.df, new.df)
          }
          
          
          
#        } # End of LRP types loop
        
#      } # End of BM types loop
      
    } # End of nReps loop 
    
    outputDir<-paste(outDir,"/DataOut/nCUCombinations",sep="")
    
    if (file.exists(outputDir) == FALSE){
      dir.create(outputDir)
    } 
    
    write.csv(output.df, paste(outputDir,"/",RunName,"_",nCUList[nn],"CUs.csv", sep="")) 
    
    
  } # End of nCUs loop
  
  
}
