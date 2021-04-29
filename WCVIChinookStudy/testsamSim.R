genericRecoverySim(simPars[i, ], cuPar=CUpars, srDat=cohoRecDatTrim,
        variableCU=FALSE, ricPars=mcmcOut, cuCustomCorrMat = corMatrix,
         nTrials=nProj, makeSubDirs=FALSE, random=FALSE, outDir=outDir)

simPar <- simPars[1,]
cuPar <- CUpars
srDat <- cohoRecDatTrim
variableCU <- FALSE
ricPars <- mcmcOut
cuCustomCorrMat <- corMatrix
nTrials <- nProj
makeSubDirs <- FALSE
random <- FALSE
outDir <- outDir
dirName <- outDir

catchDat <- NULL
larkPars <- NULL
erCorrMat <- NULL
uniqueProd <- TRUE
uniqueSurv <- FALSE

