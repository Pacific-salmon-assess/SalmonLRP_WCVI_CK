# Projection LRPs for Pacific salmon

Primary contact: Carrie Hold and Kendra Holt, DFO (carrie.holt@dfo-mpo.gc.ca; Kendra.Holt@dfo-mpo.gc.ca). Code developed by Carrie Holt, Kendra Holt & Brooke Davis, DFO. 

### Overview

This repository contains the files necessary to run projection-based Limit Reference Points (LRPs) for the West Coast Vancouver Island (WCVI) Chinook salmon stock management units (SMUs). Using this code, alternative LRP options can be derived based on the probability of all component inlets being above lower benchmarks.

Code and associated files are organized into the following sub-folders:  

1. **Code**
 * Contains files/functions used to (i) estimate LRPs, (ii) run retrospective analyses of LRP estimation over multiple years and numbers of CUs, and (iii) plot results  
 * All files in this folder are intended to be generic code that can be applied to any stock aggregate. Note, this code was developed to run for WCVI Chinook, Interior Fraser Coho, and South Coast Chum salmon, so is meant to be generic with if statements for those individual SMUs were specificity is required.
 
2. **WCVIChinookStudy**
 * Contains files specific to West Coast Vancouver Island Chinook, which was a case study in Holt, K. et al. (2023)

3. **Documents**
 * Miscellaneous word documents included for reference  



#### To run analyses, follow these steps:

1) Set working directory to "SalmonLRP_RetroEval\WCVIChinookStudy"

2) Run code in "runWCVIChinook_projLRP.r". See documentation within that file for further instructions.

#### Citations
Holt, K.R., Holt, C.A., Warkentin, L., Wor, C., Davis, B., Arbeider, M., Bokvist, J., Crowley, S., Grant, S., Luedke, W., McHugh, D., Picco, C., and Van Will, P. 2023. Case Study Applications of LRP Estimation Methods to Pacific Salmon Stock Management Units. DFO Can. Sci. Advis. Sec. Res. Doc. 2023/010. iv+129p.


