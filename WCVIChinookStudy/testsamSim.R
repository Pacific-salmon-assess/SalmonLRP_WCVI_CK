genericRecoverySim(simPars[i, ], cuPar=CUpars, srDat=cohoRecDatTrim,
        variableCU=FALSE, ricPars=mcmcOut, cuCustomCorrMat = corMatrix,
         nTrials=nProj, makeSubDirs=FALSE, random=FALSE, outDir=outDir)

simPar <- simPars[4,]
cuPar <- CUpars
srDat <- recDatTrim%>%mutate(rec2=NA, rec3=NA, rec4=NA, rec5=NA, rec6=NA)
variableCU <- FALSE
ricPars <- NULL#mcmcOut
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


recDatTrim <- recDatTrim%>%mutate(rec2=NA, rec3=NA, rec4=NA, rec5=NA, rec6=NA)
