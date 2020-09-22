This folder contains all code for the LRP retrospective analysis.

Archive - This folder contains files that are no longer in use for the main analysis, but may be usefull when looking back on model development.

JAGS_Code - This folder contails the three JAGS/Winbugs models used in the Arbeider et al. Interior Fraser Coho RPA res. doc., that were provided with that document. It also includes a basic multi-stock hierarchical Ricker model written by B Davis that was compared against these models during intitial development.

TMB_Files - This folder contains 5 versions of the integrated assessment code that calculated Ricker parameters, Sgen benchmarks for each CU, and then fits a logistic model to estimate an aggregate LRP. Each version differs mainly in the form of stock-recruitment relationship used. TMB_Example.r is an R script that can be used to run each version of the TMB model.

The other files include variance functions, as well as main file from which the IFC analysis is run (runFraserCoho.r). 

StatusPlot.r should be integrated into plotFunctions.r and runFraserCoho.r, but since this figure was still in development, it is currently separate.



