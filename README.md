# Retrospective Analysis of LRP Options

Primary contact: Kendra Holt, DFO (Kendra.Holt@dfo-mpo.gc.ca). Code developed by Kendra Holt & Brooke Davis, DFO. 

### Overview

This repository contains the files necessary to run retrospective analyses of Limit Reference Point (LRP) options for Pacific Salmon stock management units (SMUs). Using this code, alternative LRP options can be evaluated based on the level of uncertainty associated with them and sensitivity to missing data. At present, the only case study included is Interior Fraser Coho, but this may be expanded in the future.

For the most part, LRPs options available for retrospective evaluation represent the aggregate SMU-level abundance level that has historically been associated with either (i) a specified probability that all CUs will be above their lower benchmarks or (ii) a specified proportion of CUs being above their lower benchmarks. LRPs are estimated using a logistic regression fit to historical data.   

Code and associated files are organized into the following sub-folders:  

1. **Code**
 * Contains files/functions used to (i) estimate LRPs, (ii) run retrospective analyses of LRP estimation over multiple years and numbers of CUs, and (iii) plot results  
 * All files in this folder are intended to be generic code that can be applied to any stock aggregate
 
2. **IFCohoStudy**
 * Contains files specific to Interior Fraser Coho case study, including input data, outputs, sensitivity analyses, and master file from which Interior Fraser Coho retrospective analyses are run ("runFraserCoho.r")

3. **Documents**
 * Miscellaneous word documents included for reference or created while formulating case studies 


### Interior Fraser Coho Case Study

*** Note: Some of the Interior Fraser Coho data are still provisional, and are not to be used for other purposes without permission from Michael Arbeider (Michael.Arbeider@dfo-mpo.gc.ca) ***

This case study compares multiple LRP options for the Interior Fraser Coho Stock Management Unit (SMU), which is made up of 5 Conservation Units (Middle Fraser, Fraser Canyon, Lower Thompson, North Thompson, South Thompson). Lower benchmarks for each CU are obtained by fitting stock-recruitment models to CU-level data in order to estimate Sgen, which is the spawner abundance from which the CU can recover to Smsy within one generation in the absence of fishing. 

Data for this case study covered brood years 1998-2014 (return years 2001-2018). Data for the 1998-2013 brood years were taken from the 2018 Interior Fraser Coho RPA report (Appendix 4 of Arebider et al. 2020, available at http://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2020/2020_025-eng.pdf), while one additional year of spawner-recruit data (2014 brood year) and sub-population level escapement series were provided by DFO’s Fraser River Stock Assessment Unit (M. Arbeider, pers. comm., 2020). 

Aggregate LRPs are estimated using an integrated model coded in TMB that simultaneoulsy fits (i) CU-level stock-recruit models to estimate Sgen and (ii) a SMU-level logsitic regression model that estimates the aggregate abundance that has historically been associated with a specified proportion of CUs being above Sgen (binomial model) or the aggregate abundance that has historically been associated with all CUs having a specified probability of being above Sgen (bernoulli model).  

Retrospective analyses for the Interior Fraser Coho case study focus on the following choices:

1. Type of Ricker stock recruitment model used to estimate Sgen (i.e., choice of whether to include an informative prior distribution on carrying capacity)
2. The use of hierarchical versus individual modelling approaches when fitting stock-recruit models
3. Type of logistic regression used to define aggregate escapement level for LRP (bernoulli vs. binomial)
4. The choice of threshold probability (or proportion) of CUs above Sgen that is used to identify the LRP

In addition, questions about the effect of missing data on LRP estimates are evaluated by looking at:

1. How LRP estimates vary with the number of years of available data
2. How LRP estimates vary when data is only available from a subset of CUs

The first two of these options (type of Ricker Model,  hierarchical vs. individual models) are implemented by calling one of four different TMB files to estimate LRPs (SR_HierRicker_Surv.cpp, SR_HierRicker_SurvCap.cpp, SR_IndivRicker_Surv.cpp, SR_IndivRicker_SurvCap.cpp; located in "LRP_RetroEval/Code/TMB_Files"). All other options are implemented by changing data and parameter inputs to these four TMB files.

#### Threshold benchmarks based on sub-population abundance

An alternative definition of CU-level benchmarks is also considered for Interior Fraser Coho to match work that has been undertaken by the **Interior Fraser Coho Recovery Team** and incorporated into the 2018 Recovery Potential Assessment (RPA; Arbeider et al. 2020). As part of this work, short- and long-term recovery targets are based on maintaining the diversity of 11 subpopulations nested within CUs. Using this approach, lower benchmarks are based on the proportion of sub-populations within a CU that are above a 1000-fish threshold. 

* A **short-term recovery target** was identified as the 3-year geometric average, natural-origin escapement for which *at least half* of the subpopulations within each of the five CUs exceeded 1,000 spawning Coho Salmon. 

* A **long-term recovery target** was identified as the 3-year geometric average, natural-origin escapement for which *all* the subpopulations within each of the five CUs exceeded 1,000 spawning Coho Salmon.     

The Interior Fraser Coho case study includes functions that retrospectively compare and plot LRPs based on sub-population diversity (using a 1000-fish threshold) with LRPs based on CU-level diversity (using Sgen).  TMB code to calculate LRPs based on sub-population diversity is located in "LRP_RetroEval/Code/TMB_Files/ThresholdAbund_Subpop1000.cpp".  


#### To run case study:

1) Set working directory to "LRP_RetroEval\IFCohoStudy"

2) Run code in "runFraserCoho.r". This is the master file for this analysis from which all other files are sourced, and in which scenarios and plot options can be specified.


### Analysis outputs

Current model outputs, including .csv files of estimated parameters and plots, are posted to a dropbox site that can be accessed at:
https://www.dropbox.com/sh/otiku88jc2eu8cx/AACagLQd85blX7jc8-yBYs4Na?dl=0

