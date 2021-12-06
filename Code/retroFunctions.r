# ==================================================================

# Functions that loop over different subsets of data (either via annual retrospective analyses or bootstrap of years or CUs).
### For each subset of data tested, functions from LRPFunctions.r will be called and results will be saved.
#
# runAnnualRetro    -  Runs retrospective analysis using specified model type. Calls Run_LRP or Run_Ricker_LRP functions from LRPFunctions.r (latter if integratedModel=T) 
# runNCUsRetro      -  Runs retrospective analysis using specified model type. Calls runAnnualRetro function
#
# ==================================================================










#' runAnnualRetro
#'
#' Runs retrospective analysis using specified model type. Calls Run_LRP or Run_Ricker_LRP functions 
#' from LRPFunctions.r (latter if integratedModel=T) 
#' 
#' @param EscpDat data frame, ecapement data with the following columns: MU, CU_Name, CU, yr, Escp
#' @param SRDat data frame, SR dat with the following columns: CU_Name, CU_ID, BroodYear, Spawners, 
#' Recruits, Age_3_Recruits, Age_4_Recruits, STAS_Age_3, STAS_Age_4, ER_Age_3, ER_Age_4
#' @param startYr numerical, year to start the retrospective analysis (year that the shortest time 
#' series ends)
#' @param endYr numerical, year to end the retrospective analysis (most complete time series)
#' @param BroodYrLag the number of years to subtract from the current year to get to the most recent
#'  brood year
#' @param genYrs length of a generation
#' @param p default is 0.95 
#' @param BMmodel Type of model being run options are: if integratedModel = F, ThreshAbund_Subpop1000_ST, 
#' ThreshAbund_Subpop1000_LT, Percentile and if integratedModel = F other values passed to Run_Ricker_LRP 
#' @param LRPmodel BernLogistic or BinLogistic (deprecated)
#' @param LRPfile default NULL
#' @param integratedModel Logical, default FALSE, indicates if Ricker curve and logistic regression sre estimated in the same model 
#' @param useGenMean Logical, default TRUE, indicates if generational mean should e used 
#' @param TMB_Inputs List of parameter starting values and prior definitions to be passed to tmb model
#' @param outDir Case study results directory
#' @param RunName Name to identify run results and plote
#' @param bootstrapMode Logical, default False
#' @param plotLRP Logical default T
#' @param runLogisticDiag Logical default False
#' 
#' @export
#' 
runAnnualRetro<-function(EscpDat, SRDat, startYr, endYr, BroodYrLag, genYrs, p = 0.95,
                         BMmodel, LRPmodel=NULL, LRPfile=NULL, integratedModel=F, useGenMean=T, 
                         TMB_Inputs=NULL, outDir, RunName, bootstrapMode=F, plotLRP=T,
                         runLogisticDiag=F, codeDir) {
  
  yearList <- startYr:endYr
  
  # Put together Run info
  RunInfo <- data.frame(BMmodel, LRPmodel, useGenMean, integratedModel)
  
  
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
  
  # Loop over years ===============================================
  for (yy in 1:length(yearList)) {
  
      # Only use SR data for brood years that have recruited by current year yy 
       # (note: most recent brood year is calculated by subtracting BroodYearLag (e.g. 4 years) from current year)
      if (!is.null(SRDat)) {
        Dat <- SRDat %>%  filter(BroodYear <= (yearList[yy])-BroodYrLag)
      }
      EscpDat.yy <- EscpDat %>% filter(yr <= yearList[yy]) 
 
      # get trigger ready so can be used with both basic and TMB
      if (LRPmodel == "BernLogistic") useBern_Logistic <- TRUE
      if (LRPmodel == "BinLogistic") useBern_Logistic <- FALSE
      
      # Model Not Integrated (i.e., separate benchmark and LRP functions are being used) ===========================
      if(integratedModel == F){
        
        # 1) Call specified benchmark function:

        
        # # Case 1: LRP model is based on projections
        # if (LRPmodel == "Proj") {
        #   
        #   Dat$yr_num <- group_by(Dat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
        #   Dat$CU_ID <- group_by(Dat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
        #   EscDat <- EscpDat.yy %>%  right_join(unique(Dat[,c("CU_ID", "CU_Name")]))
        #   
        #   LRP_Mod<-Run_ProjRicker_LRP(SRDat = Dat, EscDat = EscDat, BMmodel = BMmodel, LRPmodel = LRPmodel, 
        #                               useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs, outDir=outDir)
        #   
        # }
        
      
        # Case 1: BM Model == 1000 fish threshold at sub-population level
        
        if (BMmodel == "ThreshAbund_Subpop1000_ST" | BMmodel == "ThreshAbund_Subpop1000_LT") { # Note that ST is for short-term threshold, LT is for long-term threshold

          # Calculate geometric means for subpopulation escapements (where, subpopulation is sum of tributaries) 
          EscpDat.cu <- EscpDat.yy %>% group_by(CU_Name, CU_ID, Subpop_Name, yr)  %>% summarise(SubpopEscp=sum(Escp)) %>%
            mutate(Gen_Mean = rollapply(SubpopEscp, genYrs, gm_mean, fill = NA, align="right"))
          # Add a column indicating whether geometric mean escapement is > 1000 fish threshold for each subpopulation 
          #Above1000<-ifelse(EscpDat.cu$Gen_Mean >= 1000, 1, 0)
          Above1000<-ifelse(EscpDat.cu$SubpopEscp >= 1000, 1, 0)
          EscpDat.cu<-EscpDat.cu %>% add_column(Above1000)

            # Calculate the number of subpopulations that are above 1000 fish threshold in each CU
          LBM_status_byCU <- EscpDat.cu %>% group_by(CU_Name, CU_ID, yr) %>% 
            summarise(Escp=sum(SubpopEscp), N = length(unique(Subpop_Name)),N_grThresh=sum(Above1000))
          
          if (BMmodel == "ThreshAbund_Subpop1000_ST") {
            # Add a column indicating whether >= 50% of subpopulations were above 1000 fish threshold in each CU (1 = yes, 0 = no) (short-term recovery goal)
            HalfGrThresh<-ifelse(LBM_status_byCU$N_grThresh >=  ceiling(LBM_status_byCU$N/2),1,0)
            LBM_status_byCU <- LBM_status_byCU %>% add_column(AboveBenchmark=HalfGrThresh)
          }
          if (BMmodel == "ThreshAbund_Subpop1000_LT") {
            # Add a column indicating whether all subpopulations were above 1000 fish threshold in each CU (long-term recovery goal) 
            AllGrThresh<-ifelse(LBM_status_byCU$N_grThresh == LBM_status_byCU$N,1,0)
            LBM_status_byCU <- LBM_status_byCU %>% add_column(AboveBenchmark=AllGrThresh)
          }
          
          # Model check: Added dum 2 as a test to confirm our results match M. Arbeider code when we take geometric mean at the sub-population level before aggregating escapements
          # --- June 5, 2020: Have confirmed that our results do match . Our results for long-term recovery target are the same as those 
          #         from the "IFC Conservation Goal Targets.Rmd" file Michael provided us.
          #dum<-EscpDat.yy %>% group_by(Subpop_Name,yr) %>% summarise(Escp.sp = sum(Escp)) %>% mutate(Gen_Mean = rollapply(Escp.sp, genYrs, gm_mean, fill = NA, align="right"))
          #dum2<- dum %>% group_by(yr) %>% summarise(AggEscp.gm=sum(Gen_Mean))
          

          # Call specified LRP function:
          LRP_Mod<-Run_LRP(Dat=LBM_status_byCU, Mod = LRPfile, useBern_Logistic = useBern_Logistic, 
                         useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs)
          # This version with dum2=dum2 as an input was created to test consistency with M. Arbeider code; 
          # ---- if wanting to run with geometric means calculated at subpopulation level, will need to add AggEscp_gmBySP arguement to LRP_Mod function
            #LRP_Mod<-Run_LRP(EscDat=LBM_status_byCU,Mod = BMmodel, useBern_Logistic = useBern_Logistic, 
            #                 useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs, dum2=dum2)
        }
       
        # Case 3: BM model is percentile-based benchmark model (e.g., for Inside South Coast Chum)
        #"if (BMmodel %in% c( "LRP_Logistic_Only", "LRP_Logistic_Only_LowAggPrior )) {
        
        if (BMmodel %in% c("Percentile")) {
          LBM_status_byCU <- EscpDat.yy %>% group_by(CU_Name) %>%  # make new column with 25% benchmark
            mutate(benchmark_perc_25= quantile(Escp, probs=0.25, na.rm=TRUE), benchmark_perc_50 = quantile(Escp, probs=0.5, na.rm=TRUE)) 
          # Need to end up with a data frame that has CU, year, and whether CU is above benchmark (1 means yes, 0 mean no)
          # new data frame that is the key to which percentiles to use for which CUs
          which_perc <- TMB_Inputs$perc_benchmark 
          # New variable with which percentile to use
          LBM_status_byCU$use_perc <- which_perc$percentile[match(x=LBM_status_byCU$CU_Name, table=which_perc$CU)]

          # make function to get status, will use correct benchmark (25 or 50)
          get_status <- function(escp, bm25, bm50, use_bm) {
            if(use_bm == 0.25)
              bm <- bm25
            if(use_bm == 0.5)
              bm <- bm50
            status <- as.integer(ifelse(escp >= bm, 1,0))
            status
          }
          # get status by mapping over four variables of the data frame
          LBM_status_byCU$AboveBenchmark <- pmap_int(list(LBM_status_byCU$Escp, LBM_status_byCU$benchmark_perc_25, 
                                                      LBM_status_byCU$benchmark_perc_50, LBM_status_byCU$use_perc),
                                                 get_status)

          # Call specified LRP function:
          LRP_Mod<-Run_LRP(Dat=LBM_status_byCU, Mod = LRPfile, useBern_Logistic = useBern_Logistic, 
                           useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs)
        }


       # Case 4: BM Model == Multi dimensional rapid tool assesment
        #multi dimensional status starts in 2000 but the data starts in 1998.
        
        if (BMmodel == "RapidAssessment_Ricker" | BMmodel == "RapidAssessment_Ricker_Cap") { 
        
         #Calculate or read in data from rapid assessment tool
          if(BMmodel == "RapidAssessment_Ricker" ){
            rapid_status<-read.csv(paste(outDir,"/DataOut/multiDimStatusEsts_Ricker.csv", sep=""))
          }
          if(BMmodel == "RapidAssessment_Ricker_Cap"){
            rapid_status<-read.csv(paste(outDir,"/DataOut/multiDimStatusEsts_Ricker_priorCap.csv", sep=""))
          }
          #filter years to allow for retrospective analysis
          rapid_status <- rapid_status %>% 
                  filter(yr <= yearList[yy]) %>% 
                  mutate(AboveBenchmark = ifelse(rapidStatus=="Red",0,1))
          
          RMA_status_byCU <- rapid_status %>% 
                              select(CU_Name, CU, yr, Escp, AboveBenchmark)%>%
                              rename(CU_ID=CU)
         

          #LBM_status_byCU <-  CU_Name       CU_ID    yr  Escp     N N_grThresh

         #run LRP logistic regression with otput from rapid assessment tool
          # Call specified LRP function:
          LRP_Mod<-Run_LRP(Dat=RMA_status_byCU, Mod = LRPfile, useBern_Logistic = useBern_Logistic, 
                         useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs)
          #

        }#end rapid assessment option

    # Integrated model (i.e., benchmark and LRP are estimated simultaneously in TMB) ===================================
    
    } else if(integratedModel == T){
      
      # Prep data frame to work with function
      
      Dat$yr_num <- group_by(Dat,BroodYear) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
      # FLAG: why isn't this in generic model run? Also need CU_ID column for percentile benchmarks, currently doing that in data processing code in runSouthCoastChum.R and runFraserCoho.r
      Dat$CU_ID <- group_by(Dat, CU_ID) %>% group_indices() - 1 # have to subtract 1 from integer so they start with 0 for TMB/c++ indexing
      
      EscDat <- EscpDat.yy %>%  right_join(unique(Dat[,c("CU_ID", "CU_Name")])) # This line is creating problems when I make CU_ID in data being read in
     
      LRP_Mod <- Run_Ricker_LRP(SRDat = Dat, EscDat = EscDat, BMmodel = BMmodel, Bern_Logistic = useBern_Logistic, 
                     useGenMean = useGenMean, genYrs = genYrs, p = p,  TMB_Inputs)
      
      # Compile Ricker parameters from TMB outputs
      RickerEsts <- LRP_Mod$All_Ests %>% filter(Param %in% c("A", "B" , "Smsy", "Sgen", "sigma"))
      # remove standard error, z and Pr z 2 values from values fron output, these are used for model diagnostics
      output <- RickerEsts %>% dplyr::select(-c(Std..Error, z.value, Pr...z.2..)) %>% 
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
    

    # Make output csv file that has 25% benchmarks and escapement for each year, for each retro year
    if (BMmodel %in% c( "Percentile" )) {
      new.perc.df <- LBM_status_byCU[, names(LBM_status_byCU) %in% c("CU_Name", "yr", "Escp", "CU", "benchmark_perc_25", "benchmark_perc_50", "use_perc", "AboveBenchmark") ]
      new.perc.df$retro_year <- yearList[yy] # add column with retrospective year
    
      if(yy==1){
        out.perc.df <- new.perc.df
      } else {
        out.perc.df <- rbind(out.perc.df, new.perc.df) # as looping through retro years, add rows with percentile benchmarks for each CU (with escapement too)
      }
    }
    
    
    if (plotLRP == T) {
      plotLogistic(Data = LRP_Mod$Logistic_Data, Preds = LRP_Mod$Preds, 
                   LRP = LRP_Mod$LRP, useGenMean = useGenMean,
                   plotName = paste("LogisticMod", yearList[yy], sep ="_"), outDir = figDir,
                   p = p, useBern_Logistic = useBern_Logistic)
    }

    
    
  
    # 
    # if (plotLRP == T & LRPmodel == "Proj") {
    # 
    #   
    #   plotProjected(Data = LRP_Mod$Proj, LRP = LRP_Mod$LRP,
    #                plotName = paste("ProjMod", yearList[yy], sep ="_"), 
    #                outDir = figDir, p = p)
    #   
    #   
    # }
    

  
    ## Run logistic model diagnostics =======
    if (runLogisticDiag==TRUE & is.na(LRP_Mod$LRP$fit) == FALSE) {
  
      SMUlogisticData <- LRP_Mod$Logistic_Data %>% rename(ppn=yy, SMU_Esc=xx, 
                                                          Years=yr)
      
<<<<<<< HEAD
      
      logisticDiagStats<-LRdiagnostics(SMUlogisticData, 
=======
      logisticDiagStats<-LRdiagnostics(SMUlogisticData=SMUlogisticData, 
>>>>>>> 8c530b5386e800ad936200b58b70fe3fd8412882
                 nCU=length(unique(EscpDat$CU_Name)), All_Ests = LRP_Mod$All_Ests,
                  p=p, Bern_logistic = useBern_Logistic, dir=paste(figDir,"/",sep=""),
                  plotname=paste("logisticFitDiag",yearList[yy], sep="_"),
                 caseStudy=unique(EscpDat$MU), codeDir)
    
      
      capture.output(logisticDiagStats, file = paste(outputDir,"/logisticFitDiag_",yearList[yy] ,".txt", sep=""))
      save(logisticDiagStats, file = paste(outputDir,"/logisticFitDiagStats_",yearList[yy] ,".rda", sep=""))
      save(LRP_Mod, file = paste(outputDir,"/logisticFit_",yearList[yy] ,".rda", sep=""))
    
    }
    
    
  } # end of year loop


 
  # Save outputs and make plots if function is not being run as part of a bootstrap
  if (bootstrapMode == F) {
  
    # Save annual retrospective outputs as csv files =========================
    write.csv(outLRP.df, paste(outputDir,"/annualRetro_LRPs.csv", sep=""))
    if (integratedModel == T) {
      write.csv(outSR.df, paste(outputDir,"/annualRetro_SRparsByCU.csv", sep=""))
    }
    if (BMmodel %in% c( "Percentile" )) { # save percentile benchmark data frame output for plotting
      write.csv(out.perc.df, paste(outputDir, "/annualRetro_perc_benchmarks.csv", sep=""))
    }
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
                                BMmodel, LRPmodel, LRPfile=NULL, integratedModel,useGenMean, TMB_Inputs, outDir, RunName,
                                  runLogisticDiag=F) {
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
      
     
      # For Subpopulation threshold models:
      if(BMmodel %in% c("ThreshAbund_Subpop1000_ST")) {
        
        # Calculate penalty for sub-population approach
        subPop_B_penalty_lwr<-1000 # set at abundance below which no one CU could be above subpop have at least half of subpops above 1000 fish
        subPop_B_penalty_upr<-ceiling(length(unique(EscpDat.ii$Subpop_Name))*0.5)*1000 # short-cut to get upper limit as half of subpopulations * 1000
        B_penalty_mu<-mean(c(subPop_B_penalty_lwr,subPop_B_penalty_upr))
        dum<-optim(par=200, fn = getSD, method="Brent",lower=1, upper=5000, low_lim=subPop_B_penalty_lwr, hi_lim=subPop_B_penalty_upr)
        B_penalty_sigma<-dum$par
        
        # Add penalty parts to TMB_Inputs for this combination of CUs
        TMB_Inputs.ii$B_penalty_mu<-B_penalty_mu
        TMB_Inputs.ii$B_penalty_sigma<-B_penalty_sigma
      
    }
      
# Run annual retrospective analyses for sampled CUs:
 
      out<- runAnnualRetro(EscpDat=EscpDat.ii, SRDat=SRDat.ii, startYr=startYr, endYr=endYr, BroodYrLag=BroodYrLag,
                           genYrs=genYrs, p = p, BMmodel = BMmodel, LRPmodel=LRPmodel, LRPfile=LRPfile, integratedModel=integratedModel,
                           useGenMean=useGenMean, TMB_Inputs=TMB_Inputs.ii, outDir=cohoDir, RunName = NA,
                           bootstrapMode = T, plotLRP=F, runLogisticDiag=runLogisticDiag)

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

      
    } # End of nReps loop 
    
    

    outputDir<-paste(outDir,"/DataOut/nCUCombinations",sep="")
    
    if (file.exists(outputDir) == FALSE){
      dir.create(outputDir)
    } 
    
    write.csv(output.df, paste(outputDir,"/",RunName,"_",nCUList[nn],"CUs.csv", sep="")) 
    
  
  } # End of nCUs loop
  
  
}
