# Projection LRPs for Pacific salmon

Primary contact: Carrie Hold and Kendra Holt, DFO (carrie.holt@dfo-mpo.gc.ca; Kendra.Holt@dfo-mpo.gc.ca). Code developed by Carrie Holt, Kendra Holt & Brooke Davis, DFO. 

### Overview

This repository contains the files necessary to run projection-based Limit Reference Points (LRPs) for the West Coast Vancouver Island (WCVI) Chinook salmon stock management units (SMUs). Using this code, alternative LRP options can be derived based on probabilt of all component inlets being above lower benchmarks.

Code and associated files are organized into the following sub-folders:  

1. **Code**
 * Contains files/functions used to (i) estimate LRPs, (ii) run retrospective analyses of LRP estimation over multiple years and numbers of CUs, and (iii) plot results  
 * All files in this folder are intended to be generic code that can be applied to any stock aggregate. Note, this code was developed to run for WCVI Chinook, Interior Fraser Coho, and South Coast Chum salmon, so is meant to be generic with if statements for those individual SMUs were specificity is required.
 
2. **WCVIChinookStudy**
 * Contains files specific to West Coast Vancouver Island Chinook case study

3. **Documents**
 * Miscellaneous word documents included for reference or created while formulating case studies 



#### To run analyses, follow these steps:

1) Set working directory to "SalmonLRP_RetroEval\WCVIChinookStudy"

2) Run code in "runWCVIChinook_projLRP.r". See documentation within that file for further instructions.


