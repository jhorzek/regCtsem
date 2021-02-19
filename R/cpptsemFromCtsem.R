#' cpptsemFromCtsem
#'
#' transforms fitted ctsem model to cpptsem model
#'
#' @param ctsemModel fittet ctsem object
#'
#' @examples
#' library(ctsemOMX)
#' library(cpptsem)
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

