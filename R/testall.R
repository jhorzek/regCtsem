testall <- function(){
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

