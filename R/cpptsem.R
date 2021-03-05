##### The following functions are used when translating a model object from
##### ctsemOMX to Rcpp

#### Core function ####

#' cpptsemFromCtsem
#'
#' transforms fitted ctsem model to cpptsem model
#'
#' @param ctsemModel fittet ctsem object
#'
#' @examples
#' library(regCtsem)
#'
#' addCINT <- FALSE
#' if(addCINT){
#'   CINT = matrix(c("cint1", "cint2"), nrow = 2, ncol = 1)
#' }else{
#'   CINT = matrix(c(0,0), nrow = 2, ncol = 1)
#' }
#' stationary <- c('T0TRAITEFFECT','T0TIPREDEFFECT')
#'
#' ## ctsem model without trait
#' AnomAuthmodel1 <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
#'                           Tpoints = 5, n.latent = 2, n.manifest = 2,
#'                           MANIFESTVAR=diag(0, 2),
#'                           TRAITVAR = NULL,
#'                           CINT = CINT)
#' AnomAuthfit1 <- ctFit(AnomAuth, AnomAuthmodel1, useOptimizer = F, stationary = stationary)
#' AnomAuthfit1$mxobj$fitfunction$result[[1]]
#' gradientModel1 <- OpenMx::mxRun(OpenMx::mxModel(AnomAuthfit1$mxobj,
#'                                                 OpenMx::mxComputeSequence(steps=list(
#'                                                   OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
#'                                                                                 hessian = FALSE))
#'                                                 )))
#' centralGrandients <- gradientModel1$compute$steps[[1]]$output[["gradient"]][,"central"]
#' names(centralGrandients) <- rownames(gradientModel1$compute$steps[[1]]$output[["gradient"]])
#' centralGrandients
#'
#' ## with cpptsem
#' cpptsemmodel1 <- cpptsemFromCtsem(ctsemModel = AnomAuthfit1)
#' cpptsemmodel1$computeRAM()
#' cpptsemmodel1$fitRAM()
#' cpptsemmodel1$m2LL
#' cpptsemmodel1$approxRAMGradients((1.1 * 10^(-16))^(1/3))[names(centralGrandients)]
#'
#' # change parameter values
#' AnomAuthfit1_1 <- ctFit(AnomAuth, AnomAuthmodel1, useOptimizer = T, stationary = stationary)
#' newParameters <- omxGetParameters(AnomAuthfit1_1$mxobj)
#' cpptsemmodel1$setParameterValues(newParameters, names(newParameters))
#' cpptsemmodel1$computeRAM()
#' cpptsemmodel1$fitRAM()
#' cpptsemmodel1$m2LL
#'
#' ## ctsem model with trait
#' AnomAuthmodel2 <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
#'                           Tpoints = 5, n.latent = 2, n.manifest = 2,
#'                           MANIFESTVAR=diag(0, 2),
#'                           TRAITVAR = "auto",
#'                           CINT = CINT)
#' AnomAuthfit2 <- ctFit(AnomAuth, AnomAuthmodel2, useOptimizer = F, stationary = stationary)
#' AnomAuthfit2$mxobj$fitfunction$result[[1]]
#' gradientModel2 <- OpenMx::mxRun(OpenMx::mxModel(AnomAuthfit2$mxobj,
#'                                                 OpenMx::mxComputeSequence(steps=list(
#'                                                   OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
#'                                                                                 hessian = FALSE))
#'                                                 )))
#' centralGrandients <- gradientModel2$compute$steps[[1]]$output[["gradient"]][,"central"]
#' names(centralGrandients) <- rownames(gradientModel2$compute$steps[[1]]$output[["gradient"]])
#' centralGrandients
#' ## with cpptsem
#' cpptsemmodel2 <- cpptsemFromCtsem(AnomAuthfit2)
#' cpptsemmodel2$computeRAM()
#' cpptsemmodel2$fitRAM()
#' cpptsemmodel2$m2LL
#' cpptsemmodel2$approxRAMGradients((1.1 * 10^(-16))^(1/3))[names(centralGrandients)]
#'
#'
#' ## Example 3: Kalman Filter
#' set.seed(175446)
#' ## define the population model:
#'
#' # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
#' ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)
#'
#' generatingModel<-ctsem::ctModel(Tpoints=20,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                                 MANIFESTVAR=diag(0,2),
#'                                 LAMBDA=diag(1,2),
#'                                 DRIFT=ct_drift,
#'                                 DIFFUSION=matrix(c(.5,0,0,.5),2),
#'                                 CINT=matrix(0,nrow = 2, ncol = 1),
#'                                 T0MEANS=matrix(0,ncol=1,nrow=2),
#'                                 T0VAR=diag(1,2), type = "omx")
#'
#' # simulate a training data and testing data set
#' traindata <- ctsem::ctGenerate(generatingModel,n.subjects = 20, wide = TRUE)
#'
#' ## Build the analysis model. Note that drift eta1_eta2 is freely estimated
#' # although it is 0 in the population.
#' myModel <- ctsem::ctModel(Tpoints=20,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                           LAMBDA=diag(1,2),
#'                           MANIFESTVAR=diag(0,2),
#'                           CINT=CINT,
#'                           DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
#'                           T0MEANS=matrix(0,ncol=1,nrow=2),
#'                           T0VAR="auto", type = "omx")
#' myModel <- ctFit(myModel, dat = traindata, objective = "Kalman", useOptimizer = F, stationary = stationary)
#' myModel$mxobj$fitfunction$result[[1]]
#'
#' gradientModel3 <- OpenMx::mxRun(OpenMx::mxModel(myModel$mxobj,
#'                                                 OpenMx::mxComputeSequence(steps=list(
#'                                                   OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
#'                                                                                 hessian = FALSE))
#'                                                 )))
#' centralGrandients <- gradientModel3$compute$steps[[1]]$output[["gradient"]][,"central"]
#' names(centralGrandients) <- rownames(gradientModel3$compute$steps[[1]]$output[["gradient"]])
#' centralGrandients
#' KalmanScores <- mxKalmanScores(myModel$mxobj)
#' KalmanScores$xUpdated[2:21,]
#'
#' ## with cpptsem
#' cpptsemmodel3 <- cpptsemFromCtsem(ctsemModel = myModel,wideData = traindata)
#' cpptsemmodel3$computeAndFitKalman()
#' cpptsemmodel3$m2LL
#' cpptsemmodel3$latentScores[1,]
#' cpptsemmodel3$approxKalmanGradients((1.1 * 10^(-16))^(1/3))[names(centralGrandients)]
#'
#' @export

cpptsemFromCtsem <- function(ctsemModel, wideData = NULL, removeD = TRUE){
  if(!"T0TRAITEFFECT" %in% ctsemModel$ctfitargs$stationary){
    stop("Non-stationary trait effect not yet implemented")
  }
  stationaryT0VAR <- "T0VAR" %in% ctsemModel$ctfitargs$stationary
  stationaryT0MEANS <- "T0MEANS" %in% ctsemModel$ctfitargs$stationary

  if(!is.null(ctsemModel$ctmodelobj$TDPREDEFFECT) | !is.null(ctsemModel$ctmodelobj$TDPREDMEANS) | !is.null(ctsemModel$ctmodelobj$TDPREDVAR)){
    stop("Time dependent predictors not yet implemented")
  }

  # step 1: extract nlatent, nmanifest and tPoints as well as mxObject
  nlatent <- ctsemModel$ctmodelobj$n.latent
  nmanifest <- ctsemModel$ctmodelobj$n.manifest
  Tpoints <- ctsemModel$ctmodelobj$Tpoints
  mxObject <- ctsemModel$mxobj

  # step 2: prepare dataset

  if(tolower(ctsemModel$ctfitargs$objective) == "mxram"){
    dataInformation <- constructDataset(wideData = mxObject$data$observed)
    dataForRAM <- prepareRAMData(dataset = dataInformation$dataset,
                                 individualMissingPatternID = dataInformation$individualMissingPatternID,
                                 uniqueMissingPatterns = dataInformation$uniqueMissingPatterns)
  }

  if(tolower(ctsemModel$ctfitargs$objective) == "kalman"){
    if(is.null(wideData)){
      stop("Kalman filter requires the wideData to be provided.")
    }
    dataInformation <- constructDatasetKalman(wideData = wideData)
    dataForKalman <- prepareKalmanData(dataset = dataInformation$dataset,
                                       nlatent = nlatent, nmanifest = nmanifest,
                                       dtIndicators = dataInformation$dtIndicators,
                                       Tpoints = Tpoints)
  }

  # step 3: extract continuous time matrices and parameters
  ctMatrices <- extractCtsemMatrices(mxObject = mxObject, nlatent = nlatent, nmanifest = nmanifest)

  # step 4: get parameter table from mxObject
  if(tolower(ctsemModel$ctfitargs$objective) == "mxram"){
    parameterTable <- extractParameterTableFromMx(mxObject = mxObject)
  }
  if(tolower(ctsemModel$ctfitargs$objective) == "kalman"){
    parameterTable <- extractParameterTableFromMx(mxObject = mxObject)
    if(removeD){
      # for some reason ctsem defines the D matrix with the same values as the ManifestMeans but does
      # use separate matrices for both. We remove the D matrix from the parameter table as the
      # implementation used by cpptsem simply adds the manifestmeans to the predicted observed variables.
      parameterTable <- subset(parameterTable, !parameterTable$matrix == "D")
    }
  }

  if(tolower(ctsemModel$ctfitargs$objective) == "mxram"){
    # step 5: prepare RAM matrices
    AMatrix <- regCtsem::prepareAMatrix(mxObject = mxObject, ctMatrices = ctMatrices, nlatent = nlatent, nmanifest = nmanifest,
                                        Tpoints = Tpoints, dT = dataInformation$dT)
    SMatrix <- regCtsem::prepareSMatrix(mxObject = mxObject, ctMatrices = ctMatrices, nlatent = nlatent, nmanifest = nmanifest,
                                        Tpoints = Tpoints, dT = dataInformation$dT, stationaryT0VAR = stationaryT0VAR)
    MMatrix <- regCtsem::prepareMMatrix(mxObject = mxObject, ctMatrices = ctMatrices, nlatent = nlatent, nmanifest = nmanifest,
                                        Tpoints = Tpoints, dT = dataInformation$dT, stationaryT0MEANS = stationaryT0MEANS)
    FMatrix <- mxObject$F$values

    # create cpptsem model
    cpptsem <- new(cpptsemRAMmodel, ctsemModel$mxobj$name, ctMatrices, parameterTable, stationaryT0VAR, stationaryT0MEANS)

    if(!is.null(AMatrix$discreteDRIFTUnique)){
      cpptsem$setDiscreteDRIFTUnique(AMatrix$discreteDRIFTUnique)
    }
    if(!is.null(AMatrix$discreteTRAITUnique)){
      cpptsem$setDiscreteTRAITUnique(AMatrix$discreteTRAITUnique)
    }
    if(!is.null(SMatrix$DRIFTHASHExponentialUnique)){
      cpptsem$setDRIFTHASHExponentialUnique(SMatrix$DRIFTHASHExponentialUnique)
    }
    if(!is.null(SMatrix$discreteDIFFUSIONUnique)){
      cpptsem$setDiscreteDIFFUSIONUnique(SMatrix$discreteDIFFUSIONUnique)
    }
    if(!is.null(MMatrix$discreteCINTUnique)){
      cpptsem$setDiscreteCINTUnique(MMatrix$discreteCINTUnique)
    }

    cpptsem$setRAMMatrices(AMatrix$AInitializer, SMatrix$SInitializer, MMatrix$MInitializer, FMatrix,
                           AMatrix$cppAParameterIndicators, SMatrix$cppSParameterIndicators, MMatrix$cppMParameterIndicators)
    cpptsem$setRAMData(dataForRAM)

    return(cpptsem)
  }

  if(tolower(ctsemModel$ctfitargs$objective) == "kalman"){
    # prepare discrete time parameters
    discreteTimeParameterConstructor <- prepareDiscreteElements(mxObject = mxObject, ctMatrices = ctMatrices,
                                                                nlatent = nlatent, nmanifest = nmanifest,
                                                                Tpoints = Tpoints, dT = as.vector(dataInformation$dT))

    discreteTimeParameterNames <- prepareDiscreteElementNames(ctMatrices, Tpoints, dT = dataInformation$dT)

    # prepare Kalman matrices to store the predicted values
    KalmanMatrices <- prepareKalmanMatrices(nlatent = nlatent, nmanifest = nmanifest,
                                            Tpoints = Tpoints, sampleSize = nrow(wideData))

    # create cpptsem model
    cpptsem <- new(cpptsemKalmanModel, ctsemModel$mxobj$name, ctMatrices, parameterTable, stationaryT0VAR, stationaryT0MEANS)

    if(!is.null(discreteTimeParameterConstructor$discreteDRIFTUnique)){
      cpptsem$setDiscreteDRIFTUnique(discreteTimeParameterConstructor$discreteDRIFTUnique)
    }
    if(!is.null(discreteTimeParameterConstructor$discreteTRAITUnique)){
      cpptsem$setDiscreteTRAITUnique(discreteTimeParameterConstructor$discreteTRAITUnique)
    }
    if(!is.null(discreteTimeParameterConstructor$DRIFTHASHExponentialUnique)){
      cpptsem$setDRIFTHASHExponentialUnique(discreteTimeParameterConstructor$DRIFTHASHExponentialUnique)
    }
    if(!is.null(discreteTimeParameterConstructor$discreteDIFFUSIONUnique)){
      cpptsem$setDiscreteDIFFUSIONUnique(discreteTimeParameterConstructor$discreteDIFFUSIONUnique)
    }
    if(!is.null(discreteTimeParameterConstructor$discreteCINTUnique)){
      cpptsem$setDiscreteCINTUnique(discreteTimeParameterConstructor$discreteCINTUnique)
    }

    cpptsem$setDiscreteTimeParameterNames(discreteTimeParameterNames)

    cpptsem$setKalmanData(dataForKalman)

    cpptsem$setKalmanMatrices(KalmanMatrices)

    return(cpptsem)
  }

  stop("Error: The objective function of the provided ctsem object was not found. Did you fit the model before providing it?")


}

#### Create Dataset ####

#' constructDataset
#'
#' Separates data and time intervals. Computes the number of unique time intervals and the unqiue missingness patterns
#' @param wideData dataset in wide format compatible with ctsem
#' @export
constructDataset <- function(wideData){
  # separate data from time intervals
  dT <- subset(wideData, select = grepl("dT", colnames(wideData)))

  if(any(apply(dT, 2, function(x) length(unique(x)))>1)){
    stop("Observed columns with more than one unqiue dT. This is not possible")
  }
  dT <- apply(dT, 2, unique)

  dataset <- subset(wideData, select = !grepl("dT", colnames(wideData)) & !grepl("intervalID", colnames(wideData)) )

  # extract unique dts
  dTUnique <- unique(dT)

  # extract unique missingness patterns

  isMissing <- is.na(dataset)
  uniqueMissingPatterns <- unique(isMissing)
  individualMissingPatternID <- c()
  for(mrow in 1:nrow(isMissing)){
    rowMissing <- apply(uniqueMissingPatterns, 1, function(x) all(isMissing[mrow,] == x))
    individualMissingPatternID <- c(individualMissingPatternID,
                                     (1:nrow(uniqueMissingPatterns))[rowMissing])
  }

  dataList <- list("dataset" = as.matrix(dataset),
                   "dT" = dT,
                   "dTUnique" = dTUnique,
                   "uniqueMissingPatterns" = uniqueMissingPatterns,
                   "individualMissingPatternID" = individualMissingPatternID
  )

  return(dataList)
}

#' prepareRAMData
#'
#' prepare data for RAM model
#'
#' @param dataset dataset
#' @param individualMissingPatternID IDs for individual missingness pattern
#' @param uniqueMissingPatterns unique missingness patterns
#' @export
prepareRAMData <- function(dataset, individualMissingPatternID, uniqueMissingPatterns){

  if(!is.matrix(uniqueMissingPatterns)){
    cols = length(uniqueMissingPatterns)
    uniqueMissingPatterns <- matrix(uniqueMissingPatterns, nrow = 1, ncol = cols)
  }

  uniqueMissingPatternIDs <- unique(individualMissingPatternID)
  patternNames <- paste0("missingnessPattern_", seq_len(length(uniqueMissingPatternIDs)))
  datasetRAM <- vector("list", length = length(uniqueMissingPatternIDs)+1)
  names(datasetRAM) <- c("names", patternNames)
  datasetRAM[["names"]] <- patternNames

  for(id in uniqueMissingPatternIDs){
    patternName <- patternNames[id]
    datasetRAM[[patternName]] <- list(
      "sampleSize" = NA,
      "nObservedVariables" = NA,
      "rawData" = NA,
      "observedMeans" = NA,
      "observedCov" = NA,
      "cppExpectedSelector" = NA,
      "expectedMeans" = NA,
      "expectedCovariance" = NA,
      "m2LL" = -999999.9)
    datasetRAM[[patternName]] [["rawData"]] <- matrix(dataset[individualMissingPatternID == id,!uniqueMissingPatterns[id,]], nrow = sum(individualMissingPatternID == id))
    datasetRAM[[patternName]] [["sampleSize"]] <- nrow(datasetRAM[[patternName]] [["rawData"]])
    datasetRAM[[patternName]] [["nObservedVariables"]] <- sum(!uniqueMissingPatterns[id,])

    if(datasetRAM[[patternName]] [["sampleSize"]] > 1){
      datasetRAM[[patternName]] [["observedMeans"]] <- matrix(apply(datasetRAM[[patternName]] [["rawData"]], 2, mean), ncol = 1)
      datasetRAM[[patternName]] [["observedCov"]] <- ((datasetRAM[[patternName]] [["sampleSize"]]-1)/(datasetRAM[[patternName]] [["sampleSize"]]))*cov(datasetRAM[[patternName]] [["rawData"]])
    }else{
      datasetRAM[[patternName]] [["rawData"]] <- t(datasetRAM[[patternName]] [["rawData"]])
      datasetRAM[[patternName]] [["observedMeans"]] <- matrix(-999999.9, nrow = datasetRAM[[patternName]] [["nObservedVariables"]], ncol = 1)
      datasetRAM[[patternName]] [["observedCov"]] <- matrix(-999999.9, nrow = datasetRAM[[patternName]] [["nObservedVariables"]], ncol = datasetRAM[[patternName]] [["nObservedVariables"]])
    }

    datasetRAM[[patternName]] [["cppExpectedSelector"]] <- regCtsem::getRowsAndColsOfNonMissing(uniqueMissingPatterns[id,])

    datasetRAM[[patternName]] [["expectedMeans"]] <- matrix(-999999.9, nrow = datasetRAM[[patternName]] [["nObservedVariables"]], ncol = 1)
    datasetRAM[[patternName]] [["expectedCovariance"]] <- matrix(-999999.9, nrow = datasetRAM[[patternName]] [["nObservedVariables"]], ncol = datasetRAM[[patternName]] [["nObservedVariables"]])
  }

  return(datasetRAM)
}

#' constructDatasetKalman
#'
#' Separates data and time intervals. Computes the number of unique time intervals and the unqiue missingness patterns
#' @param wideData dataset in wide format compatible with ctsem
#' @import mgcv
#' @export
constructDatasetKalman <- function(wideData){
  wideData <- as.matrix(wideData)
  # separate data from time intervals
  dT <- subset(wideData, select = grepl("dT", colnames(wideData)))

  dataset <- subset(wideData, select = !grepl("dT", colnames(wideData)) & !grepl("intervalID", colnames(wideData)) )

  # extract unique dts
  dTUnique <- unique(as.vector(dT))

  # make dt indicators
  dtIndicators <- matrix(NA, nrow = nrow(dT), ncol = ncol(dT))
  for(i in seq_len(length(dTUnique))){
    dtIndicators[dT == dTUnique[i]] <- i-1
  }

  # extract unique missingness patterns

  isMissing <- is.na(dataset)
  uniqueMissingPatterns <- mgcv::uniquecombs(isMissing, ordered=FALSE)
  individualMissingPatternID <- attr(uniqueMissingPatterns,"index")

  dataList <- list("dataset" = as.matrix(dataset),
                   "dT" = dT,
                   "dTUnique" = dTUnique,
                   "dtIndicators" = dtIndicators,
                   "uniqueMissingPatterns" = uniqueMissingPatterns,
                   "individualMissingPatternID" = individualMissingPatternID
  )

  return(dataList)
}

#' prepareKalmanData
#'
#' prepare data for Kalman model
#'
#' @param dataset data set
#' @param nlatent number of latent variables
#' @param nmanifest number of manifest variables
#' @param dtIndicators indicators for discrete time parameters
#' @param Tpoints number of time points
#' @import OpenMx
#' @export
prepareKalmanData <- function(dataset, nlatent, nmanifest, dtIndicators, Tpoints){

  sampleSize <- nrow(dataset)

  dataList <- vector("list", Tpoints)

  return(list("dataset" = dataset,
              "sampleSize" = sampleSize,
              "nlatent" = nlatent,
              "nmanifest" = nmanifest,
              "dtIndicators" = dtIndicators,
              "Tpoints" = Tpoints))
}


#### Extract CT Parameter ####

#' extractCtsemMatrices
#'
#' Extract DRIFT, DIFFUSION, etc from an mxObject
#' @param mxObject mxObject from ctsemOMX
#' @param nlatent number of latent variables
#' @param nmanifest number of manifest variables
#' @import OpenMx
#' @export
extractCtsemMatrices <- function(mxObject, nlatent, nmanifest){
  ctMatrices <- list()
  ## latent
  # T0MEANS
  cat("Translating model to C++. Found the following elements: ")
  if(!is.null(mxObject$T0MEANS)){
    cat("T0MEANS, ")
    T0MEANS <- list("values" = NULL, "names" = NULL)
    T0MEANS$values <- deepCopyNumericMatrix(mxObject$T0MEANS$values)
    T0MEANS$names <- deepCopyStringMatrix(mxObject$T0MEANS$labels)
    ctMatrices[["T0MEANS"]] <- T0MEANS
  }

  # T0VARbase
  if(!is.null(mxObject$T0VARbase)){
    cat("T0VARbase, ")
    T0VARbase <- list("values" = NULL, "names" = NULL)
    T0VARbase$values <- deepCopyNumericMatrix(mxObject$T0VARbase$values)
    T0VARbase$names <- deepCopyStringMatrix(mxObject$T0VARbase$labels)
    ctMatrices[["T0VARbase"]] <- T0VARbase
  }


  # DRIFT
  if(!is.null(mxObject$DRIFT)){
    cat("DRIFT, ")
    DRIFT <- list("values" = NULL, "names" = NULL)
    DRIFT$values <- deepCopyNumericMatrix(mxObject$DRIFT$values)
    DRIFT$names <- deepCopyStringMatrix(mxObject$DRIFT$labels)
    ctMatrices[["DRIFT"]] <- DRIFT
  }


  # DIFFUSIONbase
  if(!is.null(mxObject$DIFFUSIONbase)){
    cat("DIFFUSIONbase, ")
    DIFFUSIONbase <- list("values" = NULL, "names" = NULL)
    DIFFUSIONbase$values <- deepCopyNumericMatrix(mxObject$DIFFUSIONbase$values)
    DIFFUSIONbase$names <- deepCopyStringMatrix(mxObject$DIFFUSIONbase$labels)
    ctMatrices[["DIFFUSIONbase"]] <- DIFFUSIONbase
  }

  # TRAITVARbase
  if(!is.null(mxObject$TRAITVARbase$values)){
    cat("TRAITVARbase, ")
    TRAITVARbase <- list("values" = NULL, "names" = NULL)
    TRAITVARbase$values <- deepCopyNumericMatrix(mxObject$TRAITVARbase$values)
    TRAITVARbase$names <- deepCopyStringMatrix(mxObject$TRAITVARbase$labels)
    ctMatrices[["TRAITVARbase"]] <- TRAITVARbase
  }

  # T0TRAITEFFECT
  if(!is.null(mxObject$T0TRAITEFFECT)){
    cat("T0TRAITEFFECT, ")
    T0TRAITEFFECT <- list("values" = NULL, "names" = NULL)
    T0TRAITEFFECT$values <- deepCopyNumericMatrix(mxObject$T0TRAITEFFECT$values)
    T0TRAITEFFECT$names <- deepCopyStringMatrix(mxObject$T0TRAITEFFECT$labels)
    ctMatrices[["T0TRAITEFFECT"]] <- T0TRAITEFFECT
  }

  # CINT
  if(!is.null(mxObject$CINT)){
    cat("CINT, ")
    CINT <- list("values" = NULL, "names" = NULL)
    CINT$values <- deepCopyNumericMatrix(mxObject$CINT$values)
    CINT$names <- deepCopyStringMatrix(mxObject$CINT$labels)
    ctMatrices[["CINT"]] <- CINT
  }

  ## manifest

  # MANIFESTMEANS
  if(!is.null(mxObject$MANIFESTMEANS)){
    cat("MANIFESTMEANS, ")
    MANIFESTMEANS <- list("values" = NULL, "names" = NULL)
    MANIFESTMEANS$values <- deepCopyNumericMatrix(mxObject$MANIFESTMEANS$values)
    MANIFESTMEANS$names <- deepCopyStringMatrix(mxObject$MANIFESTMEANS$labels)
    ctMatrices[["MANIFESTMEANS"]] <- MANIFESTMEANS
  }

  # LAMBDA
  if(!is.null(mxObject$LAMBDA)){
    cat("LAMBDA, ")
    LAMBDA <- list("values" = NULL, "names" = NULL)
    LAMBDA$values <- deepCopyNumericMatrix(mxObject$LAMBDA$values)
    LAMBDA$names <- deepCopyStringMatrix(mxObject$LAMBDA$labels)
    ctMatrices[["LAMBDA"]] <- LAMBDA
  }

  # MANIFESTVAR
  if(!is.null(mxObject$MANIFESTVARbase$values)){
    cat("MANIFESTVARbase.")
    MANIFESTVARbase <- list("values" = NULL, "names" = NULL)
    MANIFESTVARbase$values <- deepCopyNumericMatrix(mxObject$MANIFESTVARbase$values)
    MANIFESTVARbase$names <- deepCopyStringMatrix(mxObject$MANIFESTVARbase$labels)
    ctMatrices[["MANIFESTVARbase"]] <- MANIFESTVARbase
  }

  cat("\n")
  return(ctMatrices)

}

#' extractParameterTableFromMx
#'
#' wrapper for the omxLocateParameters function from OpenMx
#' returns a data frame with parameter labels, location and values
#'
#' @param mxObject Object of type mxObject
#' @import OpenMx
#' @export
extractParameterTableFromMx <- function(mxObject){
  parameterTable <- OpenMx::omxLocateParameters(mxObject)
  return(parameterTable)
}

#' getRowsAndColsOfNonMissing
#'
#' Build indicators for missingness patterns
#' @param uniqueMissingPattern unique missingness patterns
#' @export
getRowsAndColsOfNonMissing <- function(uniqueMissingPattern){
  nonMissing <- !uniqueMissingPattern
  indicators <- c()
  for(i in seq_len(length(nonMissing))){
    if(nonMissing[i]){
      indicators <- c(indicators, i)
    }
  }
  # for cpp:
  indicators = indicators -1
  return(indicators)
}

#### Prepare Discrete Time Parameters ####

#' prepareDiscreteElementNames
#'
#' prepare the names of the discrete parameters (e.g. discrete time drift)
#'
#' @param ctMatrices continuous time matrices
#' @param Tpoints number of time points
#' @param dT time differences
#' @import OpenMx
#' @export
prepareDiscreteElementNames <- function(ctMatrices, Tpoints, dT){

  retList <- list()

  # identify unique dTs. This will later on be used to compute the corresponding discrete time values just once
  dTUnique <- unique(as.vector(dT))
  uniqueNumbers <- seq_len(length(dTUnique))
  discreteNumbering <- matrix(NA, nrow = nrow(dT), ncol = ncol(dT))
  for(i in seq_len(length(dTUnique))){
    discreteNumbering[dT == dTUnique[i]] <- uniqueNumbers[i]
  }

  ctParameterNames <- names(ctMatrices)

  if("DRIFT" %in% ctParameterNames){
    discreteDRIFTNames <- matrix(paste0("discreteDRIFT_", discreteNumbering),
                                 nrow = nrow(discreteNumbering),
                                 ncol = ncol(discreteNumbering))
    retList[["discreteDRIFTNames"]] <- discreteDRIFTNames
  }

  if("T0TRAITEFFECT" %in% ctParameterNames){
    discreteTRAITNames <- matrix(paste0("discreteTRAIT_", discreteNumbering),
                                 nrow = nrow(discreteNumbering),
                                 ncol = ncol(discreteNumbering))
    retList[["discreteTRAITNames"]] <- discreteTRAITNames
  }

  if("DIFFUSIONbase" %in% ctParameterNames){
    discreteDIFFUSIONNames <- matrix(paste0("discreteDIFFUSION_", discreteNumbering),
                                     nrow = nrow(discreteNumbering),
                                     ncol = ncol(discreteNumbering))
    DRIFTHASHExponentialNames <- matrix(paste0("DRIFTHASHExponential_", discreteNumbering),
                                        nrow = nrow(discreteNumbering),
                                        ncol = ncol(discreteNumbering))

    retList[["discreteDIFFUSIONNames"]] = discreteDIFFUSIONNames
    retList[["DRIFTHASHExponentialNames"]] = DRIFTHASHExponentialNames
  }

  if("CINT" %in% ctParameterNames){
    discreteCINTNames <- matrix(paste0("discreteCINT_", discreteNumbering),
                                nrow = nrow(discreteNumbering),
                                ncol = ncol(discreteNumbering))

    retList[["discreteCINTNames"]] <- discreteCINTNames
  }


  return(retList)
}

#' prepareDiscreteElements
#'
#' prepare matrices for discrete time elements
#'
#' @param mxObject Object of type mxObject
#' @param ctMatrices continuous time matrices
#' @param nlatent number of latent variables
#' @param nmanifest number of manifest variables
#' @param Tpoints number of time points
#' @param dT time differences
#' @import OpenMx
#' @export
prepareDiscreteElements <- function(mxObject, ctMatrices,
                                    nlatent, nmanifest, Tpoints, dT){

  retList <- list()

  # identify unique dTs. This will later on be used to compute the corresponding discrete time values just once
  dTUnique <- unique(dT)
  uniqueNumbers <- seq_len(length(dTUnique))
  discreteNumbering <- rep(NA, length(dT))
  for(i in seq_len(length(dTUnique))){
    discreteNumbering[dT == dTUnique[i]] <- uniqueNumbers[i]
  }

  ctParameterNames <- names(ctMatrices)

  if("DRIFT" %in% ctParameterNames){
    discreteDRIFTNames <- paste0("discreteDRIFT_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteDRIFTUnique <- vector("list", length = length(unique(discreteDRIFTNames))+2)
    names(discreteDRIFTUnique) <- c("labels", "dT", unique(discreteDRIFTNames))
    discreteDRIFTUnique$labels <- unique(discreteDRIFTNames)
    discreteDRIFTUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteDRIFTNames)))){
      discreteDRIFTUnique[[unique(discreteDRIFTNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }
    retList[["discreteDRIFTUnique"]] <- discreteDRIFTUnique
  }

  if("T0TRAITEFFECT" %in% ctParameterNames){
    T0TRAITEFFECTName <- "T0TRAITEFFECT"

    discreteTRAITNames <- paste0("discreteTRAIT_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteTRAITUnique <- vector("list", length = length(unique(discreteTRAITNames))+2)
    names(discreteTRAITUnique) <- c("labels", "dT", unique(discreteTRAITNames))
    discreteTRAITUnique$labels <- unique(discreteTRAITNames)
    discreteTRAITUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteTRAITNames)))){
      discreteTRAITUnique[[unique(discreteTRAITNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }
    retList[["discreteTRAITUnique"]] <- discreteTRAITUnique
  }

  if("DIFFUSIONbase" %in% ctParameterNames){
    discreteDIFFUSIONNames <- paste0("discreteDIFFUSION_", discreteNumbering)
    DRIFTHASHExponentialNames <- paste0("DRIFTHASHExponential_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteDIFFUSIONUnique <- vector("list", length = length(unique(discreteDIFFUSIONNames))+2)
    names(discreteDIFFUSIONUnique) <- c("labels", "dT", unique(discreteDIFFUSIONNames))
    discreteDIFFUSIONUnique$labels <- unique(discreteDIFFUSIONNames)
    discreteDIFFUSIONUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteDIFFUSIONNames)))){
      discreteDIFFUSIONUnique[[unique(discreteDIFFUSIONNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }

    DRIFTHASHExponentialUnique <- vector("list", length = length(unique(DRIFTHASHExponentialNames))+2)
    names(DRIFTHASHExponentialUnique) <- c("labels", "dT", unique(DRIFTHASHExponentialNames))
    DRIFTHASHExponentialUnique$labels <- unique(DRIFTHASHExponentialNames)
    DRIFTHASHExponentialUnique$dT <- dTUnique
    for(i in seq_len(length(unique(DRIFTHASHExponentialNames)))){
      DRIFTHASHExponentialUnique[[unique(DRIFTHASHExponentialNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }

    retList[["discreteDIFFUSIONUnique"]] = discreteDIFFUSIONUnique
    retList[["DRIFTHASHExponentialUnique"]] = DRIFTHASHExponentialUnique
  }

  if("CINT" %in% ctParameterNames){
    discreteCINTNames <- paste0("discreteCINT_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteCINTUnique <- vector("list", length = length(unique(discreteCINTNames))+2)
    names(discreteCINTUnique) <- c("labels", "dT", unique(discreteCINTNames))
    discreteCINTUnique$labels <- unique(discreteCINTNames)
    discreteCINTUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteCINTNames)))){
      discreteCINTUnique[[unique(discreteCINTNames)[i]]] <- matrix(0, nrow = nlatent, ncol = 1)
    }

    retList[["discreteCINTUnique"]] <- discreteCINTUnique
  }


  return(retList)
}


#### Prepare RAM Model ####

#' prepareAMatrix
#'
#' prepare matrix with directed effects for RAM model
#'
#' @param mxObject Object of type mxObject
#' @param ctMatrices continuous time matrices
#' @param nlatent number of latent variables
#' @param nmanifest number of manifest variables
#' @param Tpoints number of time points
#' @param dT time differences
#' @import OpenMx
#' @export
prepareAMatrix <- function(mxObject, ctMatrices, nlatent, nmanifest, Tpoints, dT){

  retList <- list()
  # AParameterIndicators will be used to identify the segments in the
  # A matrix which belong to the discreteDRIFTs discreteTRAITs and LAMBDAS

  AParameterIndicators <- data.frame()
  AParameterIndicatorsLabels <- c()

  # identify unique dTs. This will later on be used to compute the corresponding discrete time values just once
  dTUnique <- unique(dT)
  uniqueNumbers <- seq_len(length(dTUnique))
  discreteNumbering <- rep(NA, length(dT))
  for(i in seq_len(length(dTUnique))){
    discreteNumbering[dT == dTUnique[i]] <- uniqueNumbers[i]
  }

  ctParameterNames <- names(ctMatrices)

  if("DRIFT" %in% ctParameterNames){
    discreteDRIFTNames <- paste0("discreteDRIFT_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteDRIFTUnique <- vector("list", length = length(unique(discreteDRIFTNames))+2)
    names(discreteDRIFTUnique) <- c("labels", "dT", unique(discreteDRIFTNames))
    discreteDRIFTUnique$labels <- unique(discreteDRIFTNames)
    discreteDRIFTUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteDRIFTNames)))){
      discreteDRIFTUnique[[unique(discreteDRIFTNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }

    AParameterIndicatorsLabels <- c(AParameterIndicatorsLabels, discreteDRIFTNames)
  }


  if("T0TRAITEFFECT" %in% ctParameterNames){
    T0TRAITEFFECTName <- "T0TRAITEFFECT"

    discreteTRAITNames <- paste0("discreteTRAIT_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteTRAITUnique <- vector("list", length = length(unique(discreteTRAITNames))+2)
    names(discreteTRAITUnique) <- c("labels", "dT", unique(discreteTRAITNames))
    discreteTRAITUnique$labels <- unique(discreteTRAITNames)
    discreteTRAITUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteTRAITNames)))){
      discreteTRAITUnique[[unique(discreteTRAITNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }

    AParameterIndicatorsLabels <- c(AParameterIndicatorsLabels, discreteTRAITNames)
  }

  if("LAMBDA" %in% ctParameterNames){
    LAMBDANames <- rep("LAMBDA", Tpoints)

    AParameterIndicatorsLabels <- c(AParameterIndicatorsLabels, LAMBDANames)
  }



  AParameterIndicators <- data.frame("label" = AParameterIndicatorsLabels,
                                     "row_start" = NA,
                                     "row_end" = NA,
                                     "col_start" = NA,
                                     "col_end" = NA
  )

  # A matrix
  #AInitializer <- matrix(0, nrow = (nlatent)*Tpoints +nlatent# DRIFTPart and TRAITPart
  #                       + nmanifest*Tpoints # MANIFESTPart
  #                       ,
  #                       ncol = (nlatent)*Tpoints +nlatent # eta <- eta and eta <- trait
  #                       + nmanifest*Tpoints # y <- eta
  #)

  AInitializer <- mxObject$A$values

  ## segment A

  counterOffset <- 0

  # segment DRIFTs
  if("DRIFT" %in% ctParameterNames){
    discreteDRIFTRowStart <- seq(from = nlatent+1, # skip initial time point
                                 to = (Tpoints-1)*nlatent+1,
                                 by = nlatent)
    discreteDRIFTColStart <- seq(from = 1,
                                 to = (Tpoints-2)*nlatent+1,
                                 by = nlatent)

    for(discreteDRIFT in seq_len(length(discreteDRIFTNames)) + counterOffset){
      if(!(AParameterIndicators[discreteDRIFT,"label"] == discreteDRIFTNames[discreteDRIFT-counterOffset])){
        stop("Wrong element selected in segment DRIFTs")
      }
      AParameterIndicators[discreteDRIFT,"row_start"] <- discreteDRIFTRowStart[discreteDRIFT-counterOffset]
      AParameterIndicators[discreteDRIFT,"row_end"] <- discreteDRIFTRowStart[discreteDRIFT-counterOffset] + nlatent -1
      AParameterIndicators[discreteDRIFT,"col_start"] <- discreteDRIFTColStart[discreteDRIFT-counterOffset]
      AParameterIndicators[discreteDRIFT,"col_end"] <- discreteDRIFTColStart[discreteDRIFT-counterOffset] + nlatent -1
    }
    counterOffset <- length(discreteDRIFTNames)

    retList[["discreteDRIFTUnique"]] <- discreteDRIFTUnique
  }

  # segment TRAITs
  if("T0TRAITEFFECT" %in% ctParameterNames){
    addTrait <- nlatent
    if(any(mxObject$T0TRAITEFFECT$free)){
      stop("free T0TRAITEFFECT detected. This is not implemented.")
    }
    AInitializer[1:nlatent, (nlatent*Tpoints+1) : (nlatent*Tpoints+nlatent)] <- ctMatrices$T0TRAITEFFECT$values

    discreteTRAITRowStart <- seq(from = nlatent+1, # skip initial time point
                                 to = (Tpoints-1)*nlatent+1,
                                 by = nlatent)
    discreteTRAITColStart <- rep(nlatent*Tpoints+1, Tpoints-1)

    for(discreteTRAIT in seq_len(length(discreteTRAITNames)) + counterOffset){
      if(!(AParameterIndicators[discreteTRAIT,"label"] == discreteTRAITNames[discreteTRAIT-counterOffset])){
        stop("Wrong element selected in segment TRAITs")
      }
      AParameterIndicators[discreteTRAIT,"row_start"] <- discreteTRAITRowStart[discreteTRAIT-counterOffset]
      AParameterIndicators[discreteTRAIT,"row_end"] <- discreteTRAITRowStart[discreteTRAIT-counterOffset] + nlatent -1
      AParameterIndicators[discreteTRAIT,"col_start"] <- discreteTRAITColStart[discreteTRAIT-counterOffset]
      AParameterIndicators[discreteTRAIT,"col_end"] <- discreteTRAITColStart[discreteTRAIT-counterOffset] + nlatent -1
    }
    counterOffset <- counterOffset + length(discreteTRAITNames)

    retList[["discreteTRAITUnique"]] <- discreteTRAITUnique
  }else{
    addTrait <- 0
  }

  if("LAMBDA" %in% ctParameterNames){
    LAMBDARowStart <- seq(from = Tpoints*(nlatent) + addTrait + 1,
                          to = Tpoints*(nlatent) + addTrait + (Tpoints-1)*nmanifest + 1,
                          by = nmanifest)
    LAMBDAColStart <- seq(from = 1,
                          to = (Tpoints-1)*(nlatent) + 1,
                          by = nlatent)

    for(LAMBDA in seq_len(length(LAMBDANames)) + counterOffset){
      if(!(AParameterIndicators[LAMBDA,"label"] == LAMBDANames[LAMBDA-counterOffset])){
        stop("Wrong element selected in segment TRAITs")
      }
      AParameterIndicators[LAMBDA,"row_start"] <- LAMBDARowStart[LAMBDA-counterOffset]
      AParameterIndicators[LAMBDA,"row_end"] <- LAMBDARowStart[LAMBDA-counterOffset] + nmanifest -1
      AParameterIndicators[LAMBDA,"col_start"] <- LAMBDAColStart[LAMBDA-counterOffset]
      AParameterIndicators[LAMBDA,"col_end"] <- LAMBDAColStart[LAMBDA-counterOffset] + nlatent -1
    }
  }

  # cpp indices start at 0
  cppAParameterIndicators <- AParameterIndicators
  cppAParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] <- cppAParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] - 1

  retList[["AParameterIndicators"]] = AParameterIndicators
  retList[["cppAParameterIndicators"]] = cppAParameterIndicators
  retList[["AInitializer"]] = AInitializer

  return(retList)
}

#' prepareSMatrix
#'
#' prepare matrix with covariances for RAM model
#'
#' @param mxObject Object of type mxObject
#' @param ctMatrices continuous time matrices
#' @param nlatent number of latent variables
#' @param nmanifest number of manifest variables
#' @param dtIndicators indicators for discrete time parameters
#' @param Tpoints number of time points
#' @param dT time intervals
#' @param stationaryT0VAR boolean: are variances stationary?
#' @export
prepareSMatrix <- function(mxObject, ctMatrices, nlatent, nmanifest, Tpoints, dT, stationaryT0VAR){
  retList <- list()
  # SParameterIndicators will be used to identify the segments in the
  # S matrix which belong to the discreteDIFFUSION MANIFESVAR and T0VAR

  SParameterIndicators <- data.frame()
  SParameterIndicatorsLabels <- c()

  # identify unique dTs. This will later on be used to compute the corresponding discrete time values just once
  dTUnique <- unique(dT)
  uniqueNumbers <- seq_len(length(dTUnique))
  discreteNumbering <- rep(NA, length(dT))
  for(i in seq_len(length(dTUnique))){
    discreteNumbering[dT == dTUnique[i]] <- uniqueNumbers[i]
  }

  ctParameterNames <- names(ctMatrices)

  if("DIFFUSIONbase" %in% ctParameterNames){
    discreteDIFFUSIONNames <- paste0("discreteDIFFUSION_", discreteNumbering)
    DRIFTHASHExponentialNames <- paste0("DRIFTHASHExponential_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteDIFFUSIONUnique <- vector("list", length = length(unique(discreteDIFFUSIONNames))+2)
    names(discreteDIFFUSIONUnique) <- c("labels", "dT", unique(discreteDIFFUSIONNames))
    discreteDIFFUSIONUnique$labels <- unique(discreteDIFFUSIONNames)
    discreteDIFFUSIONUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteDIFFUSIONNames)))){
      discreteDIFFUSIONUnique[[unique(discreteDIFFUSIONNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }

    DRIFTHASHExponentialUnique <- vector("list", length = length(unique(DRIFTHASHExponentialNames))+2)
    names(DRIFTHASHExponentialUnique) <- c("labels", "dT", unique(DRIFTHASHExponentialNames))
    DRIFTHASHExponentialUnique$labels <- unique(DRIFTHASHExponentialNames)
    DRIFTHASHExponentialUnique$dT <- dTUnique
    for(i in seq_len(length(unique(DRIFTHASHExponentialNames)))){
      DRIFTHASHExponentialUnique[[unique(DRIFTHASHExponentialNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }

    SParameterIndicatorsLabels <- c(SParameterIndicatorsLabels, discreteDIFFUSIONNames)
  }

  if("MANIFESTVARbase"%in% ctParameterNames){
    MANIFESTVARNames <- rep("MANIFESTVAR", Tpoints)
    SParameterIndicatorsLabels <- c(SParameterIndicatorsLabels, MANIFESTVARNames)
  }
  if("T0VARbase" %in% ctParameterNames | stationaryT0VAR){
    T0VARNames <- "T0VAR"
    SParameterIndicatorsLabels <- c(SParameterIndicatorsLabels, T0VARNames)
  }

  if("TRAITVARbase" %in% ctParameterNames){
    TRAITVARNames <- "TRAITVAR"
    SParameterIndicatorsLabels <- c(SParameterIndicatorsLabels, TRAITVARNames)
  }

  SParameterIndicators <- data.frame("label" = SParameterIndicatorsLabels,
                                     "row_start" = NA,
                                     "row_end" = NA,
                                     "col_start" = NA,
                                     "col_end" = NA
  )

  ## Initialize S

  #SInitializer <- matrix(0, nrow = (nlatent)*Tpoints +nlatent # DRIFTPart and TRAITPart
  #                       + nmanifest*Tpoints # MANIFESTPart
  #                       ,
  #                       ncol = (nlatent)*Tpoints + nlatent # eta <-> eta and eta <- trait
  #                       + nmanifest*Tpoints # y <-> eta
  #)
  SInitializer <- mxObject$S$values

  ## segment S

  counterOffset <- 0

  if("DIFFUSIONbase" %in% ctParameterNames){
    # segment discreteDIFFUSION
    discreteDIFFUSIONsRowStart <- seq(from = nlatent+1, # skip initial time point
                                      to = (Tpoints-1)*nlatent+1,
                                      by = nlatent)
    discreteDIFFUSIONsColStart <- seq(from = nlatent+1, # skip initial time point
                                      to = (Tpoints-1)*nlatent+1,
                                      by = nlatent)

    for(discreteDIFFUSIONs in seq_len(length(discreteDIFFUSIONNames)) + counterOffset){
      if(!(SParameterIndicators[discreteDIFFUSIONs,"label"] == discreteDIFFUSIONNames[discreteDIFFUSIONs-counterOffset])){
        stop("Wrong element selected in segment discreteDIFFUSIONs")
      }
      SParameterIndicators[discreteDIFFUSIONs,"row_start"] <- discreteDIFFUSIONsRowStart[discreteDIFFUSIONs-counterOffset]
      SParameterIndicators[discreteDIFFUSIONs,"row_end"] <- discreteDIFFUSIONsRowStart[discreteDIFFUSIONs-counterOffset] + nlatent -1
      SParameterIndicators[discreteDIFFUSIONs,"col_start"] <- discreteDIFFUSIONsColStart[discreteDIFFUSIONs-counterOffset]
      SParameterIndicators[discreteDIFFUSIONs,"col_end"] <- discreteDIFFUSIONsColStart[discreteDIFFUSIONs-counterOffset] + nlatent -1
    }
    counterOffset <- length(discreteDIFFUSIONNames)

    retList[["discreteDIFFUSIONUnique"]] = discreteDIFFUSIONUnique
    retList[["DRIFTHASHExponentialUnique"]] = DRIFTHASHExponentialUnique
  }

  if("MANIFESTVARbase" %in% ctParameterNames){
    if("TRAITVARbase" %in% ctParameterNames){
      traitVarOffset <- nlatent
    }else{
      traitVarOffset <- 0
    }

    # segment MANIFESTVARs
    MANIFESTVARRowStart <- seq(from = Tpoints*(nlatent) + traitVarOffset + 1,
                               to = Tpoints*(nlatent) + traitVarOffset + (Tpoints-1)*nmanifest + 1,
                               by = nmanifest)
    MANIFESTVARColStart <- seq(from = Tpoints*(nlatent) + traitVarOffset + 1,
                               to = Tpoints*(nlatent) + traitVarOffset + (Tpoints-1)*nmanifest + 1,
                               by = nmanifest)

    for(MANIFESTVAR in seq_len(length(MANIFESTVARNames)) + counterOffset){
      if(!(SParameterIndicators[MANIFESTVAR,"label"] == MANIFESTVARNames[MANIFESTVAR-counterOffset])){
        stop("Wrong element selected in segment TRAITs")
      }
      SParameterIndicators[MANIFESTVAR,"row_start"] <- MANIFESTVARRowStart[MANIFESTVAR-counterOffset]
      SParameterIndicators[MANIFESTVAR,"row_end"] <- MANIFESTVARRowStart[MANIFESTVAR-counterOffset] + nmanifest -1
      SParameterIndicators[MANIFESTVAR,"col_start"] <- MANIFESTVARColStart[MANIFESTVAR-counterOffset]
      SParameterIndicators[MANIFESTVAR,"col_end"] <- MANIFESTVARColStart[MANIFESTVAR-counterOffset] + nmanifest -1
    }

    counterOffset <- counterOffset + length(MANIFESTVARNames)
  }

  if("T0VARbase" %in% ctParameterNames | stationaryT0VAR){
    # segment T0VAR
    SParameterIndicators[counterOffset+1,"row_start"] <- 1
    SParameterIndicators[counterOffset+1,"row_end"] <- nlatent
    SParameterIndicators[counterOffset+1,"col_start"] <- 1
    SParameterIndicators[counterOffset+1,"col_end"] <- nlatent
  }

  if("TRAITVARbase" %in% ctParameterNames){
    # segment TRAITVAR
    SParameterIndicators[counterOffset+2,"row_start"] <- nlatent*Tpoints + 1
    SParameterIndicators[counterOffset+2,"row_end"] <- nlatent*Tpoints + nlatent
    SParameterIndicators[counterOffset+2,"col_start"] <- nlatent*Tpoints + 1
    SParameterIndicators[counterOffset+2,"col_end"] <- nlatent*Tpoints + nlatent
  }

  # cpp indices start at 0
  cppSParameterIndicators <- SParameterIndicators
  cppSParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] <- cppSParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] - 1

  retList[["SParameterIndicators"]] = SParameterIndicators
  retList[["cppSParameterIndicators"]] = cppSParameterIndicators
  retList[["SInitializer"]] = SInitializer

  return(retList)

}

#' prepareFMatrix
#'
#' prepare filter matrix for RAM model
#'
#' @param nlatent number of latent variables
#' @param nmanifest number of manifest variables
#' @param Tpoints number of time points
#' @import OpenMx
#' @export
prepareFMatrix <- function(nlatent, nmanifest, Tpoints){

  Fmatrix <- matrix(0,
                    nrow = nmanifest*Tpoints,
                    ncol = Tpoints*nlatent + nlatent + Tpoints*nmanifest)

  diag(Fmatrix[,(Tpoints*nlatent + nlatent + 1):(Tpoints*nlatent + nlatent + Tpoints*nmanifest)]) <- 1

  return(Fmatrix)
}


#' prepareMMatrix
#'
#' prepare matrix with means for RAM model
#'
#' @param mxObject Object of type mxObject
#' @param ctMatrices continuous time matrices
#' @param nlatent number of latent variables
#' @param nmanifest number of manifest variables
#' @param dtIndicators indicators for discrete time parameters
#' @param Tpoints number of time points
#' @param dT time intervals
#' @param stationaryT0MEANS boolean: are Means stationary?
#' @export
prepareMMatrix <- function(mxObject, ctMatrices, nlatent, nmanifest, Tpoints, dT, stationaryT0MEANS){

  retList <- list()
  # MParameterIndicators will be used to identify the segments in the
  # M matrix which belong to the CINT and MANIFESTMEANS
  MParameterIndicators <- data.frame()
  MParameterIndicatorsLabels <- c()

  # identify unique dTs. This will later on be used to compute the corresponding discrete time values just once
  dTUnique <- unique(dT)
  uniqueNumbers <- seq_len(length(dTUnique))
  discreteNumbering <- rep(NA, length(dT))
  for(i in seq_len(length(dTUnique))){
    discreteNumbering[dT == dTUnique[i]] <- uniqueNumbers[i]
  }

  ctParameterNames <- names(ctMatrices)

  if("T0MEANS" %in% ctParameterNames | stationaryT0MEANS){
    MParameterIndicatorsLabels <- c(MParameterIndicatorsLabels, "T0MEANS")
  }

  if("CINT" %in% ctParameterNames){
    discreteCINTNames <- paste0("discreteCINT_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteCINTUnique <- vector("list", length = length(unique(discreteCINTNames))+2)
    names(discreteCINTUnique) <- c("labels", "dT", unique(discreteCINTNames))
    discreteCINTUnique$labels <- unique(discreteCINTNames)
    discreteCINTUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteCINTNames)))){
      discreteCINTUnique[[unique(discreteCINTNames)[i]]] <- matrix(0, nrow = nlatent, ncol = 1)
    }
    MParameterIndicatorsLabels <- c(MParameterIndicatorsLabels, discreteCINTNames)
  }

  if("MANIFESTMEANS" %in% ctParameterNames){
    MANIFESTMEANSNames <- rep("MANIFESTMEANS", Tpoints)
    MParameterIndicatorsLabels <- c(MParameterIndicatorsLabels, MANIFESTMEANSNames)
  }


  ## initialize M

  #MInitializer <- matrix(0, nrow = (nlatent)*Tpoints + nlatent + nmanifest*Tpoints, ncol = 1)
  MInitializer  <- matrix(mxObject$M$values, ncol = 1)
  rownames(MInitializer) <- names(mxObject$M$values)
  # create indicators for T0MEANS and CINT
  MParameterIndicators <- data.frame("label" = MParameterIndicatorsLabels,
                                     "row_start" = NA,
                                     "row_end" = NA,
                                     "col_start" = NA,
                                     "col_end" = NA
  )

  if("T0MEANS" %in% ctParameterNames | stationaryT0MEANS){

    # T0MEANS
    MParameterIndicators[1,"row_start"] <- 1
    MParameterIndicators[1,"row_end"] <- nlatent
    MParameterIndicators[1,"col_start"] <- 1
    MParameterIndicators[1,"col_end"] <- 1

    offset <- 1
  }else{
    offset <- 0
  }

  if("CINT" %in% ctParameterNames){

    # discreteCINT
    MParameterIndicators[1:length(discreteCINTNames)+offset,"row_start"] <- seq(nlatent+1, Tpoints*nlatent, nlatent)
    MParameterIndicators[1:length(discreteCINTNames)+offset,"row_end"] <- seq(nlatent+1, Tpoints*nlatent, nlatent)+nlatent-1
    MParameterIndicators[1:length(discreteCINTNames)+offset,"col_start"] <- 1
    MParameterIndicators[1:length(discreteCINTNames)+offset,"col_end"] <- 1

    retList[["discreteCINTUnique"]] <- discreteCINTUnique
    offset <- offset + length(discreteCINTNames)
  }

  if("MANIFESTMEANS" %in% ctParameterNames){
    if(!is.null(mxObject$T0TRAITEFFECT)){
      traitOffset <- nlatent
    }else{
      traitOffset <- 0
    }
    # MANIFESTMEANS
    MParameterIndicators[(offset+1) : (offset + length(MANIFESTMEANSNames) ),"row_start"] <- seq(from = Tpoints * nlatent + traitOffset + 1,
                                                                                                 to = Tpoints * nlatent + traitOffset + nmanifest*(Tpoints-1) + 1,
                                                                                                 by = nmanifest)
    MParameterIndicators[(offset+1) : (offset + length(MANIFESTMEANSNames) ),"row_end"] <- seq(from = Tpoints * nlatent + traitOffset + 1,
                                                                                               to = Tpoints * nlatent + traitOffset + nmanifest*(Tpoints-1) + 1,
                                                                                               by = nmanifest) + nmanifest - 1
    MParameterIndicators[(offset+1) : (offset + length(MANIFESTMEANSNames) ),"col_start"] <- 1
    MParameterIndicators[(offset+1) : (offset + length(MANIFESTMEANSNames) ),"col_end"] <- 1

  }

  # cpp indices start at 0
  cppMParameterIndicators <- MParameterIndicators
  cppMParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] <- cppMParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] - 1

  retList[["MParameterIndicators"]] = MParameterIndicators
  retList[["cppMParameterIndicators"]] = cppMParameterIndicators
  retList[["MInitializer"]] = MInitializer

  return(retList)

}




#### Prepare Kalman Model ####

#' prepareKalmanMatrices
#'
#' prepare matrices for Kalman model
#'
#' @param nlatent number of latent variables
#' @param nmanifest number of manifest variables
#' @param dtIndicators indicators for discrete time parameters
#' @param Tpoints number of time points
#' @param sampleSize number of persons
#' @export
prepareKalmanMatrices <- function(nlatent, nmanifest, Tpoints, sampleSize){
  latentScores <- matrix(-9999999, nrow = sampleSize, ncol = nlatent * Tpoints)
  predictedManifestValues <- matrix(-9999999, nrow = sampleSize, ncol = nmanifest * Tpoints)

  return(list(
    "latentScores" = latentScores,
    "predictedManifestValues" = predictedManifestValues
  ))
}


#### Optimization ####

#' approx_KalmanM2LLCpptsem
#'
#' computes gradients for an approximate optimization of regularized ctsem based on cpptsem with Full Information Maximum Likelihood objective
#' @param parameters paramneter values
#' @param cpptsemmodel model from cpptsem
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @author Jannik Orzek
#' @export
approx_RAMM2LLCpptsem <- function(parameters, cpptsemmodel, failureReturns){
  cpptsemmodel$setParameterValues(parameters, names(parameters))
  # catching all errors from cpptsemmodel
  # when parameter values are impossible
  invisible(capture.output(RAM <- try(cpptsemmodel$computeRAM(),
                                      silent = TRUE),
                           type = "message"))
  invisible(capture.output(FIT <- try(cpptsemmodel$fitRAM(),
                                      silent = TRUE),
                           type = "message"))
  if(class(RAM) == "try-error" | class(FIT) == "try-error"){
    return(failureReturns)
  }
  m2LL <- cpptsemmodel$m2LL
  if(is.na(m2LL) | is.infinite(m2LL)){
    return(failureReturns)
  }
  return(m2LL)
}


#' approx_RAMRegM2LLCpptsem
#'
#' approximates the regularized likelihood function using cpptsem and Full Information Maximum Likelihood
#' @param parameters parameter values
#' @param cpptsemmodel model returned from cpptsem
#' @param adaptiveLassoWeights vector with weights of the adaptive lasso
#' @param N sample size
#' @param lambda tuning parameter lambda
#' @param regIndicators string vector with names of regularized parameters
#' @param epsilon tuning parameter for epsL1 approximation
#' @param objective ML or Kalman
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @author Jannik Orzek
#' @export
approx_RAMRegM2LLCpptsem <- function(parameters, cpptsemmodel, adaptiveLassoWeights, N, lambda, regIndicators, epsilon, objective, failureReturns){
  cpptsemmodel$setParameterValues(parameters, names(parameters))
  # catching all errors from cpptsemmodel
  # when parameter values are impossible
  invisible(capture.output(RAM <- try(cpptsemmodel$computeRAM(),
                                      silent = TRUE),
                           type = "message"))
  invisible(capture.output(FIT <- try(cpptsemmodel$fitRAM(),
                                      silent = TRUE),
                           type = "message"))
  if(class(RAM) == "try-error" | class(FIT) == "try-error"){
    return(failureReturns)
  }
  m2LL <- cpptsemmodel$m2LL
  regM2LL <- m2LL + sum(N*lambda*adaptiveLassoWeights[regIndicators] * abs(parameters[regIndicators]))
  if(is.na(regM2LL) | is.infinite(regM2LL)){
    return(failureReturns)
  }
  return(regM2LL)
}


#' approx_KalmanM2LLCpptsem
#'
#' computes gradients for an approximate optimization of regularized ctsem based on cpptsem with Kalman objective
#' @param parameters paramneter values
#' @param cpptsemmodel model from cpptsem
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @author Jannik Orzek
#' @export
approx_KalmanM2LLCpptsem <- function(parameters, cpptsemmodel, failureReturns){
  cpptsemmodel$setParameterValues(parameters, names(parameters))
  # catching all errors from cpptsemmodel
  # when parameter values are impossible
  invisible(capture.output(FIT <- try(cpptsemmodel$computeAndFitKalman(),
                                      silent = TRUE),
                           type = "message"))
  if(class(FIT) == "try-error"){
    return(failureReturns)
  }
  m2LL <- cpptsemmodel$m2LL
  if(is.na(m2LL) | is.infinite(m2LL)){
    return(failureReturns)
  }
  return(m2LL)
}


#' approx_KalmanRegM2LLCpptsem
#'
#' approximates the regularized likelihood function using cpptsem and Kalman
#' @param parameters parameter values
#' @param cpptsemmodel model returned from cpptsem
#' @param adaptiveLassoWeights vector with weights of the adaptive lasso
#' @param N sample size
#' @param lambda tuning parameter lambda
#' @param regIndicators string vector with names of regularized parameters
#' @param epsilon tuning parameter for epsL1 approximation
#' @param objective ML or Kalman
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @author Jannik Orzek
#' @export
approx_KalmanRegM2LLCpptsem <- function(parameters, cpptsemmodel, adaptiveLassoWeights, N, lambda, regIndicators, epsilon, objective, failureReturns){
  cpptsemmodel$setParameterValues(parameters, names(parameters))
  # catching all errors from cpptsemmodel
  # when parameter values are impossible
  invisible(capture.output(FIT <- try(cpptsemmodel$computeAndFitKalman(),
                                      silent = TRUE),
                           type = "message"))
  if(class(FIT) == "try-error"){
    return(failureReturns)
  }
  m2LL <- cpptsemmodel$m2LL
  regM2LL <- m2LL + sum(N*lambda*adaptiveLassoWeights[regIndicators] * abs(parameters[regIndicators]))
  if(is.na(regM2LL) | is.infinite(regM2LL)){
    return(failureReturns)
  }
  return(regM2LL)
}

#' exact_getCppGradients
#'
#' exact_getCppGradients will try to compute gradients with decreasing precision starting from the default in OpenMx. Sometimes the gradients will result in NA for a very specific setting of the precisions. Then it can help to slightly alter the precision
#'
#' @param cppmodel model of type cpptsem
#' @param objective ml or Kalman
#' @export
exact_getCppGradients <- function(cppmodel, objective){
  # will try different precisions for the gradients
  defaultPrecision <- OpenMx::imxAutoOptionValue("Gradient step size")
  defaultPrecision2 <- (1.1 * 10^(-16))^(1/3)
  precisions <- rev(seq(defaultPrecision, defaultPrecision2, length.out = 5))
  for(precision in precisions){
    if(tolower(objective) == "ml"){
      invisible(capture.output(gradients <- try(cppmodel$approxRAMGradients(precision), silent = T), type = "message"))
    }else{
      invisible(capture.output(gradients <- try(cppmodel$approxKalmanGradients(precision), silent = T), type = "message"))
    }
    if(!(any(class(gradients) == "try-error")) &
       !anyNA(gradients)){
      break
    }
  }
  return(gradients)
}


#' approx_gradCpptsem
#'
#' computes gradients for an approximate optimization of regularized ctsem based on cpptsem
#' @param parameters parameter values
#' @param adaptiveLassoWeights vector with weights of the adaptive lasso
#' @param N sample size
#' @param lambda tuning parameter lambda
#' @param regIndicators string vector with names of regularized parameters
#' @param epsilon tuning parameter for epsL1 approximation
#' @param maxit maximal number of iterations
#' @param objective ML or Kalman
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @param testGradients should be tested if the final parameters result in NA gradients?
#' @author Jannik Orzek
#' @export
approx_gradCpptsem <- function(parameters, cpptsemmodel, adaptiveLassoWeights, N, lambda, regIndicators, epsilon, objective, failureReturns){
  invisible(capture.output(grad <- try(exact_getCppGradients(cppmodel = cpptsemmodel, objective = objective),
                                       silent = TRUE),
                           type = "message"))
  if(class(grad) == "try-error" || anyNA(grad)){
    ret <- rep(failureReturns, length = length(parameters))
    names(ret) <- names(parameters)
    return(ret)
  }
  grad[regIndicators] <- grad[regIndicators] +
    N*lambda*adaptiveLassoWeights[regIndicators] * parameters[regIndicators]/sqrt(parameters[regIndicators]^2+epsilon)
  return(grad[names(parameters)])
}


#' approx_cpptsemOptim
#'
#' creates an approximate solution to regularized ctsem using optim with BFGS
#' @param cpptsemmodel model returned from cpptsem
#' @param regM2LLCpptsem regularized fitting function
#' @param gradCpptsem function for computing the gradients of regM2LLCpptsem
#' @param startingValues starting values for optimization
#' @param adaptiveLassoWeights vector with weights of the adaptive lasso
#' @param N sample size
#' @param lambda tuning parameter lambda
#' @param regIndicators string vector with names of regularized parameters
#' @param epsilon tuning parameter for epsL1 approximation
#' @param maxit maximal number of iterations
#' @param objective ML or Kalman
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @param testGradients should be tested if the final parameters result in NA gradients?
#' @author Jannik Orzek
#' @export
approx_cpptsemOptim <- function(cpptsemmodel,
                                regM2LLCpptsem,
                                gradCpptsem,
                                startingValues,
                                adaptiveLassoWeights,
                                N, lambda,
                                regIndicators,
                                epsilon,
                                maxit,
                                objective,
                                failureReturns = NA,
                                testGradients){
  CpptsemFit <- optim(par = startingValues,
                      fn = regM2LLCpptsem,
                      gr = gradCpptsem,
                      cpptsemmodel, adaptiveLassoWeights, N, lambda, regIndicators, epsilon, objective, failureReturns,
                      method = "BFGS",
                      control = list(maxit = maxit))

  if(testGradients){
    grad <- try(gradCpptsem(parameters = CpptsemFit$par,
                            cpptsemmodel = cpptsemmodel,
                            adaptiveLassoWeights = adaptiveLassoWeights,
                            N =  N,lambda =  lambda, regIndicators = regIndicators,
                            epsilon = epsilon, objective =  objective,
                            failureReturns =  failureReturns))
    if(any(class(grad) == "try-error") || anyNA(grad)){
      stop("NA in gradients")
    }
  }

  return(list("parameters" = CpptsemFit$par,
              "regM2LL" = CpptsemFit$value))
}


#' approx_cpptsemSolnp
#'
#' creates an approximate solution to regularized ctsem using solnp
#' @param cpptsemmodel model returned from cpptsem
#' @param regM2LLCpptsem regularized fitting function
#' @param gradCpptsem function for computing the gradients of regM2LLCpptsem
#' @param startingValues starting values for optimization
#' @param adaptiveLassoWeights vector with weights of the adaptive lasso
#' @param N sample size
#' @param lambda tuning parameter lambda
#' @param regIndicators string vector with names of regularized parameters
#' @param epsilon tuning parameter for epsL1 approximation
#' @param maxit maximal number of iterations
#' @param objective ML or Kalman
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @param testGradients should be tested if the final parameters result in NA gradients?
#' @author Jannik Orzek
#' @import Rsolnp
#' @export
approx_cpptsemSolnp <- function(cpptsemmodel,
                                regM2LLCpptsem,
                                gradCpptsem,
                                startingValues,
                                adaptiveLassoWeights,
                                N, lambda,
                                regIndicators,
                                epsilon,
                                maxit,
                                objective,
                                failureReturns = 1e24,
                                testGradients){
  CpptsemFit <- Rsolnp::solnp(pars = startingValues,
                              fun = regM2LLCpptsem,
                              eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL,
                              ineqUB = NULL, LB = NULL, UB = NULL, control = list(trace = 0, outer.iter = maxit),
                              cpptsemmodel, adaptiveLassoWeights, N, lambda, regIndicators, epsilon, objective, failureReturns)

  if(testGradients){
    grad <- try(gradCpptsem(parameters = CpptsemFit$pars,
                            cpptsemmodel = cpptsemmodel,
                            adaptiveLassoWeights = adaptiveLassoWeights,
                            N =  N,lambda =  lambda, regIndicators = regIndicators,
                            epsilon = epsilon, objective =  objective,
                            failureReturns =  failureReturns))
    if(any(class(grad) == "try-error") || anyNA(grad)){
      stop("NA in gradients")
    }
  }

  return(list("parameters" = CpptsemFit$pars,
              "regM2LL" = CpptsemFit$values[length(CpptsemFit$values)]))
}

#' testall_cpptsem
#'
#' test cpptsem
#' @return
testall_cpptsem <- function(){
  library(cpptsem)
  library(ctsemOMX)
  data(AnomAuth)
  AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
                           Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = "auto")
  AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel, useOptimizer = F)
  mxObject <- AnomAuthfit$mxobj

  # construct dataset and time intervals
  dat <- regCtsem::constructDataset(wideData = AnomAuth)
  datRAM <- regCtsem::prepareRAMData(dataset = dat$dataset, individualMissingPatternID = dat$individualMissingPatternID, uniqueMissingPatterns = dat$uniqueMissingPatterns)
  # extract ct matrices
  ctMatrices <- regCtsem::extractCtsemMatrices(AnomAuthfit$mxobj, nlatent = 2, nmanifest = 2)
  # extract parameter table
  parTab <- regCtsem::extractParameterTableFromMx(AnomAuthfit$mxobj)
  # prepare RAM matrices
  AMatrix <- regCtsem::prepareAMatrix(mxObject = mxObject, ctMatrices = ctMatrices, nlatent = 2, nmanifest = 2, Tpoints = 5, dT = dat$dT)
  SMatrix <- regCtsem::prepareSMatrix(mxObject = mxObject, ctMatrices = ctMatrices, nlatent = 2, nmanifest = 2, Tpoints = 5, dT = dat$dT)
  MMatrix <- regCtsem::prepareMMatrix(mxObject = mxObject, ctMatrices = ctMatrices, nlatent = 2, nmanifest = 2, Tpoints = 5, dT = dat$dT)
  FMatrix <- mxObject$F$values

  # show model class:
  show(cpptsemmodel)

  # create new model
  model1 <- new(cpptsemmodel, "mymodel", ctMatrices, parTab)

  # add data
  model1$setData(as.matrix(dat$dataset), dat$dT)

  # check if parameters can be changed:
  currentParameters <- omxGetParameters(AnomAuthfit$mxobj)
  parameterLabels <- names(currentParameters)
  parameterValues <- rnorm(length(currentParameters))

  model1$setParameterValues(parameterValues, parameterLabels)

  model1$ctMatrixList$DRIFT$values[1,2] == parameterValues[parameterLabels == "drift_eta1_eta2"]

  # reset parameters
  model1$setParameterValues(omxGetParameters(AnomAuthfit$mxobj), parameterLabels)

  # check DRIFT computation
  discreteDRIFT <- regCtsem::computeDiscreteDRIFTs(DRIFTValues = AnomAuthfit$mxobj$DRIFT$values, discreteDRIFTUnique = AMatrix$discreteDRIFTUnique)
  round(discreteDRIFT$discreteDRIFT_1 - expm(AnomAuthfit$mxobj$DRIFT$values),5) == 0
  round(discreteDRIFT$discreteDRIFT_2 - expm(AnomAuthfit$mxobj$DRIFT$values*2),5) == 0

  # check TRAIT computation
  discreteTRAIT <- regCtsem::computeDiscreteTRAITs(discreteDRIFTUnique = AMatrix$discreteDRIFTUnique, discreteTRAITUnique = AMatrix$discreteTRAITUnique)
  round(discreteTRAIT$discreteTRAIT_1 - mxObject$discreteTRAIT_T1$result, 5) == 0
  round(discreteTRAIT$discreteTRAIT_2 - mxObject$discreteTRAIT_T3$result, 5) == 0

  # check A population
  A = regCtsem::fillA(A = AMatrix$AInitializer,
                      hasDiscreteDRIFTUnique = TRUE, discreteDRIFTUnique = discreteDRIFT,
                      hasDiscreteTRAITUnique = T, discreteTRAITUnique = discreteTRAIT,
                      LAMBDA = ctMatrices$LAMBDA$values,
                      AParameterIndicators = AMatrix$cppAParameterIndicators)

  # check compute function
  model1$setDiscreteDRIFTUnique(AMatrix$discreteDRIFTUnique)
  model1$setDiscreteTRAITUnique(AMatrix$discreteTRAITUnique)
  model1$setDRIFTHASHExponentialUnique(SMatrix$DRIFTHASHExponentialUnique)
  model1$setDiscreteDIFFUSIONUnique(SMatrix$discreteDIFFUSIONUnique)
  model1$setDiscreteCINTUnique(MMatrix$discreteCINTUnique)
  model1$setRAMMatrices(AMatrix$AInitializer, SMatrix$SInitializer, MMatrix$MInitializer, FMatrix,
                        AMatrix$cppAParameterIndicators, SMatrix$cppSParameterIndicators, MMatrix$cppMParameterIndicators)
  model1$setRAMData(datRAM)
  model1$computeRAM()
  # check DRIFT
  model1$DRIFTValues
  mxObject$DRIFT$values
  # check DIFFUSION
  model1$DIFFUSIONValues
  mxObject$DIFFUSION$result
  # check T0VAR
  model1$T0VARValues
  mxObject$T0VAR$result
  # show discrete drift
  model1$discreteDRIFTUnique
  # show discrete trait
  model1$discreteTRAITUnique
  # show A
  model1$A

  # show S
  model1$S

  # show M
  model1$M

  # show expected covariance
  model1$expectedCovariance

  # show expected means
  model1$expectedMeans

  # compute m2LL
  model1$fitRAM()
  model1$m2LL
  mxObject$fitfunction$result[[1]]

  # compute gradients
  parBefore = model1$getParameterValues()
  gradients = model1$approxRAMGradients(.000001)
  parAfter = model1$getParameterValues()
  if(any((parAfter - parBefore)>0)){
    stop("Parametervalues changed when computing the gradient!")
  }

  gradientModel <- mxRun(OpenMx::mxModel(mxObject,
                                         OpenMx::mxComputeSequence(steps=list(OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
                                                                                                            hessian = FALSE))
                                         )))

  mxGradients <- gradientModel$compute$steps[[1]]$output$gradient[,"central"]
  names(mxGradients) <- rownames( gradientModel$compute$steps[[1]]$output$gradient)
  if(any(!(round(gradients - mxGradients[names(gradients)], 2) == 0))){
    stop("Different gradients from OpenMx and cpptsem!")
  }
}


