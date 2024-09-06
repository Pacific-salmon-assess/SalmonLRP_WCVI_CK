---
output:
  pdf_document: default
  html_document: default
---
# Projection-based reference points for Pacific salmon

Primary contact: Carrie Holt and Kendra Holt, DFO (carrie.holt@dfo-mpo.gc.ca; Kendra.Holt@dfo-mpo.gc.ca). Code developed by Carrie Holt & Kendra Holt, DFO. 

### Overview

This readme file provides steps to run projection-based reference points for West coast Vancouver Island (WCVI) Chinook salmon. The files needed to run these analyses are primarily in this repository "SalmonLRP_WCVI_CK", but some inputs files are generated in a second repository, "Watershed-Area-Model" repository, which contains files and code required to estimate benchmarks for WCVI Chinook, as noted below. 

The "SalmonLRP_WCVI_CK" repository is forked from the "SalmonLRP_RetroEval" repository which provides code used to estimate projection-based reference points that are documented in Holt, K. et al. (2023), including WCVI Chinook as a case study.

The analyses and results from this repository are reported in Brown et al. (in revision). Citations are listed below.

Code and associated files are organized into the following sub-folders:  

1. **Code**
 * Contains files/functions used to estimate projection LRPs. All files in this folder are intended to be generic code that can be applied to any stock aggregate. This code was developed for Holt et al. (2023) to run for case studies on WCVI Chinook, Interior Fraser Coho, and South Coast Chum salmon, with if statements for individual SMUs were specificity is required.
 
2. **WCVIChinookStudy**
 * Contains files specific to West Coast Vancouver Island Chinook, which was a case study in Holt, K. et al. (2023), and the focus of Brown et al. (in revision). The Rproj files is run from this folder.
 This folder contains subfolders:
 
 **DataIn**: Folder of input data required for analyses described below in steps 3-7
 
 **DataOut**: Folder of output data derived from analyses described below in steps 3-7.
 
 **Figures**: Folder of figures generated from analyses described below in steps 3-7, and projection-based reference points in steps 9-12.
 
 **RmdReports**: Folder of R Markdown reports generated for Holt, K. et al (2023) case study analysis on WCVI Chinook
 
 **samSimInputs**: Folder of inputs required to projection-based reference points, executed using samSim R package (see Step 9 below).
 
 **samSimOutputs**: Folder of outputs generated from estimation of projection-based reference points (see Step 9 below)

3. **Documents**
 * Miscellaneous documents included for reference (currently not used).  



#### To run analyses, follow these steps:

**Step 1)** Open Rproj file in the 'SalmonLRP_WCVI_CK\WCVIChinookStudy' folder.
File: 'SalmonLRP_WCVI_CK\WCVIChinookStudy\WCVIChinookStudy.Rproj'

**Step 2)** Open 'runWCVIChinook_projLRP.r' file and run the 'Set-up:Libraries and Source Code' section
File: runWCVIChinook_projLRP.r (Section: Set-up:Libraries and Source Code)

The following steps (3-8) are required to generate inputs for projection-based reference points. These steps have been completed and saved to the repository, so analysts may skip to step 9 in this list below.


**Step 3)** Estimate CU-specific stock-recruitment parameters based on unpublished run reconstruction (D. Dobson pers. comm.). A standard Ricker model is applied.

*File*: 'runWCVICHinook_projLRP.R' (Section 13)

*Inputs*: 

'/DataIn/WCVI_SRbyCU.csv' (summary of data provided by D. Dobson, 2020)

*Outputs*: 

'/DataIn/riclogArr.csv' (log alpha values by CU based on run reconstruction, rr)

'pl' is the object that provides all Ricker parameter estimates for each CU.

*Required for*: File of inputs for projections 'samSimInputs/CUpars_AllExMH.csv', where AllExMH specifies the default indicators used in Brown et al. (in press) for projections, which is all indicators except those with major hatchery facilities. Ricker parameter values are manually input into 'CUpars_AllExMH.csv' (see step 7 below).

**Step 4)** Calculate mean proportions at age in recruitment and estimate variance in proportions at age (tau parameter of multivariate logistic distribution) from the run reconstruction for WCVI Chinook (D. Dobson pers. comm.) by CU. Note, all inlets within a given CU have the same tau parameter and mean proportions at age.

*File*: 'runWCVICHinook_projLRP.r' (Section 8)

*Inputs*: 

'samSimInputs/CUPars_AllExMH.csv' for Inlet names

'DataIn/CUages.csv' for proportions at age from run reconstruction

*Outputs*: 

'inletTau' object providing inlet-specific tau estimates. 

*Required for*: File of inputs for projections 'samSimInputs/CUpars_AllExMH.csv'. Tau values and the mean proportions of age by CU/inlet are manually input into samSimInputs/CUpars_AllExMH.csv. See step (8) below.

**Step 5)** Calculate aggregate spawner abundances summed across indicators (except major hatchery facilities) to inlets, used to estimate pairwise correlations in abundances between inlets required for projections (see Step 6 below).  The summed spawner abundances are calculated in the 'Watershed-Area-Model' repository. 

*File*: R/Inlet_Sum.R ('Watershed-Area-Model' repository)

*Inputs*: 

'DataIn/WCVIstocks.csv'(Watershed-Area-Model' repository) for list of indicator stocks

'DataIn/WCVIEsc.csv' (Watershed-Area-Model' repository) for spawner abundances by indicator population, provided by Diana McHugh (DFO South Coast Area, 2021, pers. comm.)
        
*Outputs*: 

'DataOut/Inlet_Sum_AllExMH.csv' (Watershed-Area-Model repository). This file is copied to the 'SalmonLRP_WCVI_CK' repository, location: 'DataIn/Inlet_Sum_AllExMH.csv'.

*Required for*: Estimating pairwise correlations in abundances between inlets required for projections (see Step 6 below)


**Step 6)** Generate a bubble plot of pairwise correlations in spawner abundances between inlets. 

*File*: runWCVIChinook_projLRP.R file (See top of Section 3, which refers to Fig. 6 in RmdReports/WCVI_Projections.Rmd for the creation of this plot).

*Inputs*:

'DataIn/Inlet_Sum_AllExMH.csv' for summed spawner abundances to the inlet level, generated in Step 5 above

'samSimInputs/CUPars_AllExMH.csv' for names of inlets

*Outputs*: 

'Figures/SpawnerCorrelation_AllExMH.png', a bubble plot of pairwise correlations in spawner abundances

*Required for*: Figure presented in Brown et al. (in revision). Note, the correlation matrix is calculated again when the projection-based reference points are generated (Step 9 below) based on the same input files.

**Step 7)** Manually complete input files of CU-specific parameters for projections.

*File*: samSimInput/CUPars_AllExMH.csv

*Inputs*: 

column manUnit = WCVIChinook (stock management unit, SMU)

column stk = (each inlet numbered sequentially)

column model = ricker (default stock-recruitment model assumption)

column minER = 0.05 (minimum exploitation rate assuming minimal unavoidable level of bycatch)

column usERscalar = 1 (scalar applied for additional US harvest, assumed default 1)

column canERscalar = 1 (scalar applied for additional CDN harvest, assumed default 1)

column tauCycAge = tau parameter of the multivariate logistic distribution of proportion of ages in recruitment, provided from Step 4 above, by inlet (all CUs within an inlet are identical)

column alpha = ln alpha parameter in Ricker model by inlet, provided from Step 3 above, by inlet (all CUs within an inlet are identical)

column beta0 = beta parameter of the Ricker model by inlet, calculated from lnalpha/SREP, where lnalpha is generated in Step 3 and SREP from accessible watershed area model for each inlet (file of accessible watershed-area based benchmarks: 'DataIn/WCVI_SMSY_AllExMH.csv').

column sigma = SD of Ricker residuals by inlet, provided from Step 3 above, by inlet (all CUs within an inlet are identical)

column aveGenTime = 4 (average generation time in years)

column ageFirstRec = 2 (age at first recruitment)

column ageMaxRec = 4 (generation time used to set length of the initialization period [=ageMaxRec*10], and is set to 4 here)

columns meanRec2, meanRec3, meanRec4, meanRec5, meanRec6 = mean proportion of ages in recruitment, provided from Step 4 above, by inlet (all CUs within an inlet are identical)

column Sinit = initial spawner abundances, set to spawners at equilibrium generated from accessible watershed area model, by inlet, predicted values from IWAM.r.

The following columns are not used:  domCycle, cvER, coef1, covarInit, mu_logCovar1, sig_logCovar1, min_logCovar, max_logCovar, larkAlpha, larkBeta0, larkBeta1, larkBeta2, larkBeta3, larkSigma, medianRec, lowQRec, highQRec, meanDBE, sdDBE, medMA, meanForecast, sdForecast 

*Output*: updated samSimInput/CUPars_AllExMH.csv

*Required for*: Generating random draws of Ricker model to be used in projections (Step 8 below) and running projections for projection-based reference points (Step 9 below)

**Step 8)** Generate a series of random draws for Ricker model to be used in the projections. These are generated and saved a priori to allow for comparison in outputs for scenarios when bias adjustment from back-transformation of log-normal deviates is not included (default is to include this bias adjustment) 

*File*: 'runWCVIChinook_projLRP.R' (Section 7) 

*Inputs*:

'samSimInputs/CUPars_AllExMH.csv' for inlet names, and productivity  (ln alpha) and sigma (SD of Ricker recruitment deviations) by inlet

'DataIn/WCVI_SMSY_AllExMH.csv' for SREP estimates based on accessible watershed area model

*Outputs*:

'SamSimInputs/Ricker_mcmc_AllExMH.csv' file of random draws for Ricker parameters

'Figures/AlphaDensity.png', figure of density in ln alpha values from random draws

'Figures/SREPDensity.png', figure of density of SREP values from random draws

*Required for*: Running projection-based reference points (Step 9 below)

**Step 9)** Run projections to derive projection-based reference points. This code requires R package, samSim (branch LRP), which is sourced in the code below. See samSim repository (https://github.com/Pacific-salmon-assess/samSim/tree/LRP) for code, and Holt, K. et al. (2023) Appendix B for equations and documentation.

*File*: runWCVIChinook_projLRP.R file (See Section 3). This section sources the generic function, run_scenarioProj() in the Code subfolder to generate projections.


*Inputs*:

'samSimInputs/CUPars_AllExMH.csv' (see Step 7 above)

'samSimInputs/Ricker_mcmc_AllExMH.csv' (see Step 8 above)

*Outputs*:

'SamSimOutputs/diagnostics/baseER_AllExMH/baseER_AllExMH_baseER_singleTrialFig.pdf' Example trajectory from a single Monte Carlo trial in sub-directory, labeled by scenarioName 

'SamSimOutputs/simData/projLRPDat_baseER_AllExMH.csv' Projection data to estimate projection based reference point, where file is labeled by the scenarioName

'SamSimOutputs/simData/projSpawnDat_baseER_AllExMH.csv' Projected spawner-recruit time-series by CU (or inlet), where file is labeled by the scenarioName

'SamSimOutputs/simData/baseER_AllExMH/' standard outputs from samSim, where the sub-directory is labelled by scenarioName. See repository, 'samSim' (LRP branch) for code, and Holt, K. et al. (2023) Appendix B for equations and documentation.

*Required for*: Generating projection-based reference points and figure of projections

**Step 10)** Compile output file documenting whether individual inlets were above their lower benchmarks for each MC trial (this step is supplemental and is only needed to plot individual inlets on aggregate reference points plot, see 'Required for' section below).

*File*: 'CUppnLB.R'

*Inputs*:

'SamSimOutputs/simData/baseER_AllExMH/baseER_AllExMH_baseER_CUaboveLB.RData', created in step 9 above.

'SamSimOutputs/simData/baseER_AllExMH/baseER_AllExMH_baseER_lrpdat.csv', created in step 9 above.

*Outputs*:

'SamSimOutputs/simData/projCUBenchDat_baseER_AllExMH.csv'

*Required for*: Plotting projection-based reference points with specified probability of all inlets being above lower benchmarks along a gradient in aggregate spawner abundances, with the probabilities of individual inlets being above their lower benchmarks overlaid (Step 11)

**Step 11)** Plot projection results: binned aggregate abundances against the proportion of Monte Carlo trials in that bin where all component inlets were above their lower benchmark. From this plot, projection-based reference points can be identified at various probability levels. Overlain on this plot is the probability of individual inlets being above their lower benchmark.

*File*: runWCVIChinook_projLRP.r (see Section 12).

*Inputs*: 

'SamSimOutputs/simData/projLRPDat_baseER_AllExMH.csv', output from samSim that shows in which Monte Carlo trial and year all inlets were above their lower benchmarks.

'SamSimOutputs/simData/projCUBenchDat_baseER_AllExMH.csv', output file derived in step 10.

'SamSimInputs/CUPars_AllExMH.csv' for list of inlet names

*Outputs*:

'Figures/ProjectedLRPs/baseER_AllExMH-ProjLRPCurve-ALLp.png' for figure showing projection-based reference point for various propbabilities of all inlets being above lower benchmarks, with individual inlet probabilities included.

'DataOut/ProjectedLRPs/ProjectedLRPsbaseER_AllExMH_ALLp.csv' for projection-based reference points at various probabilities

'DataOut/ProjectedLRPs/ProjectedLRP_databaseER_AllExMH_Allp.csv' for underlying binned data used to generate Figure of probabilities of all inlets being above lower benchmarks along gradient in aggregate abundances, above.

**Step 12)** Run sensitivity analyses for projection-based reference points, generated in Step 9. Note, these were run extensively in Holt, K. et al. (2023), but not for Brown et al. (in revision).

*File*: runWCVIChinook_projLRP.r (see Section 4)

*Inputs*: as in Step 9, for various alternative assumptions

*Outputs*: as in Step 9, for various alternative assumptions


#### Citations
Brown, N. et al. (in revision). West Coast of Vancouver Island Natural-Origin Chinook Salmon (Oncorhynchus tshawytscha) Stock Assessment. CSAS Working Paper20xx/nnn.

Holt, K.R., Holt, C.A., Warkentin, L., Wor, C., Davis, B., Arbeider, M., Bokvist, J., Crowley, S., Grant, S., Luedke, W., McHugh, D., Picco, C., and Van Will, P. 2023. Case Study Applications of LRP Estimation Methods to Pacific Salmon Stock Management Units. DFO Can. Sci. Advis. Sec. Res. Doc. 2023/010. iv+129p.


