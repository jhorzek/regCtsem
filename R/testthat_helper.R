#' checkAutoKFold
#'
#' used in testthat to check the automatic k-fold cross-validation feature of regCtsem
#' @param ctInit init object from ctsem
#' @param regCtsemObject object from regCtsem
checkAutoKFold <- function(ctInit, regCtsemObject){
  if(regCtsemObject$setup$autoCV != "kFold") stop("Only valid for kFold cross-validation")

  nModels <- regCtsemObject$setup$k

  for(model in 1:nModels){
    print(paste("Checking model", model, "of", nModels))
    testSamples <- regCtsemObject$cvFoldsAndModels$cvFolds[[model]]
    trainSamples <- (1:nrow(regCtsemObject$setup$dataset))[!(1:nrow(regCtsemObject$setup$dataset) %in% testSamples)]
    testSet <- regCtsemObject$cvFoldsAndModels$testSets[[model]]

    # check samples
    fullData <- regCtsemObject$setup$dataset

    if(any(!regCtsemObject$subModels[[model]]$setup$dataset == fullData[trainSamples,])){
      stop("Wrong sample size in train set")
    }

    trainModel <- ctFit(dat = fullData[trainSamples,], ctmodelobj = ctInit, useOptimizer = FALSE, objective = ifelse(regCtsemObject$setup$objective == "ML", "mxRAM", "Kalman"))

    cvModel <- ctFit(dat = fullData[testSamples,], ctmodelobj = ctInit, useOptimizer = FALSE, objective = ifelse(regCtsemObject$setup$objective == "ML", "mxRAM", "Kalman"))

    works <- checkFI(mxObject = trainModel$mxobj, regCtsemObject = regCtsemObject$subModels[[model]], cvModel = cvModel$mxobj)
  }

  print("Automatic Cross-Validation Works")
  return(works)
}

#' checkFI
#'
#' used in testthat to check the computation of fit indices of regCtsem
#' @param mxObject object from openmx
#' @param regCtsemObject object from regCtsem
#' @param cvModel optional cross-validation model
checkFI <- function(mxObject, regCtsemObject, cvModel = NULL){
  works <- FALSE
  ## get parameter labels
  parameterLabels <- names(OpenMx::omxGetParameters(mxObject))
  # separate fit and parameters
  parameters <- regCtsemObject$parameterEstimatesRaw[parameterLabels,]
  fit <- regCtsemObject$fit

  # iterate over lambdas
  for(i in 1:ncol(fit)){
    currentModel <- OpenMx::omxSetParameters(mxObject, labels = parameterLabels, values = parameters[,i])
    currentModel <- OpenMx::mxRun(currentModel, useOptimizer = F, silent = T)
    m2LL <- currentModel$fitfunction$result[[1]]
    if(abs(fit["m2LL",i] - m2LL)>.01){
      warning(paste0("Wrong m2LL for lambda = ", colnames(fit)[i]))
      return(works)
    }
    if(!is.null(cvModel)){

      cvModel <- OpenMx::omxSetParameters(cvModel, labels = parameterLabels, values = parameters[,i])

      cvModel <- OpenMx::mxRun(cvModel, useOptimizer = F, silent = T)
      cvM2LL <- cvModel$fitfunction$result[[1]]
      if(abs(fit["cvM2LL",i] - cvM2LL)>.01){
        warning("Wrong cvM2LL")
        return(works)
      }
    }
  }
  works <- TRUE

  print("FI Computation Works!")
  return(works)
}

