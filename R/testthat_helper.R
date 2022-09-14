#' checkAutoKFold
#'
#' used in testthat to check the automatic k-fold cross-validation feature of regCtsem
#' @param ctInit init object from ctsem
#' @param regCtsemObject object from regCtsem
#' @param threshold how close to zero should differences be to be treated as zero?
#' @param testIC should information criteria be tested?
#' @keywords internal
checkAutoKFold <- function(ctInit,
                           regCtsemObject,
                           threshold,
                           testIC){
  if(regCtsemObject$setup$autoCV != "kFold") stop("Only valid for kFold cross-validation")

  nModels <- regCtsemObject$setup$k

  works <- TRUE

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

    cvModel <- ctFit(dat = fullData[testSamples,],
                     ctmodelobj = ctInit,
                     useOptimizer = FALSE,
                     objective = ifelse(regCtsemObject$setup$objective == "ML",
                                        "mxRAM",
                                        "Kalman"))

    works <- works & checkFI(mxObject = trainModel$mxobj,
                     regCtsemObject = regCtsemObject$subModels[[model]],
                     cvModel = cvModel$mxobj,
                     threshold,
                     testIC)
  }

  print("Automatic Cross-Validation Works")
  return(works)
}

#' checkAutoBlocked
#'
#' used in testthat to check the automatic k-fold blocked cross-validation feature of regCtsem
#' @param ctInit init object from ctsem
#' @param regCtsemObject object from regCtsem
#' @param threshold how close to zero should differences be to be treated as zero?
#' @param testIC should information criteria be tested?
#' @keywords internal
checkAutoBlocked <- function(ctInit,
                             regCtsemObject,
                             threshold,
                             testIC){
  if(regCtsemObject$setup$autoCV != "Blocked") stop("Only valid for Blocked cross-validation")

  nModels <- regCtsemObject$setup$k

  TPoints <- regCtsemObject$setup$ctsemObject$ctmodelobj$Tpoints

  works <- TRUE

  for(model in 1:nModels){
    print(paste("Checking model", model, "of", nModels))
    testSamples <- regCtsemObject$cvFoldsAndModels$cvFolds[[model]]
    trainSamples <- (1:TPoints)[!(1:TPoints %in% testSamples)]
    testSet <- regCtsemObject$cvFoldsAndModels$testSets[[model]]

    # check samples
    fullData <- regCtsemObject$setup$dataset
    trainSet <- fullData
    trainSet[,grepl(pattern = paste0(paste0("[A-Z0-9]+_T", testSamples, "$"), collapse = "|"), colnames(trainSet))] <- NA

    if(any(regCtsemObject$cvFoldsAndModels$trainSets[[model]] != trainSet, na.rm = TRUE)){
      stop("Wrong sample size in train set")
    }

    trainModel <- ctFit(dat = trainSet,
                        ctmodelobj = ctInit,
                        useOptimizer = FALSE,
                        objective = ifelse(regCtsemObject$setup$objective == "ML", "mxRAM", "Kalman"))

    cvModel <- ctFit(dat = testSet,
                     ctmodelobj = ctInit,
                     useOptimizer = FALSE,
                     objective = ifelse(regCtsemObject$setup$objective == "ML",
                                        "mxRAM",
                                        "Kalman"))

    works <- works & checkFI(mxObject = trainModel$mxobj,
                     regCtsemObject = regCtsemObject$subModels[[model]],
                     cvModel = cvModel$mxobj,
                     threshold,
                     testIC)
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
#' @param threshold how close to zero should differences be to be treated as zero?
#' @param testIC should information criteria be tested?
#' @keywords internal
checkFI <- function(mxObject,
                    regCtsemObject,
                    cvModel = NULL,
                    threshold,
                    testIC){
  works <- FALSE
  ## get parameter labels
  parameterLabels <- names(OpenMx::omxGetParameters(mxObject))
  # separate fit and parameters
  parameters <- regCtsemObject$parameterEstimatesRaw[parameterLabels,]
  fit <- regCtsemObject$fit

  # iterate over lambdas
  for(i in 1:ncol(fit)){
    free <- rep(TRUE, length(parameterLabels))
    names(free) <- parameterLabels
    free[parameterLabels[parameters[parameterLabels,i] == 0]] <- FALSE
    currentModel <- OpenMx::omxSetParameters(mxObject,
                                             labels = parameterLabels,
                                             values = parameters[,i],
                                             free = free
    )
    currentModel <- OpenMx::mxRun(currentModel,
                                  useOptimizer = F,
                                  silent = T)
    m2LL <- currentModel$fitfunction$result[[1]]
    AICmx <- stats::AIC(currentModel)

    if(!is(object = currentModel$expectation, class2 = "MxExpectationStateSpace")){
      BICmx <- stats::BIC(currentModel)
    }else{
      BICmx <- currentModel$fitfunction$result[[1]] + log(length(unique(currentModel$data$observed[,"id"]))) * length(omxGetParameters(currentModel))
    }

    if(abs(fit["m2LL",i] - m2LL)>threshold){
      warning(paste0("Wrong m2LL for lambda = ", colnames(fit)[i]))
      return(works)
    }
    if(testIC && abs(fit["AIC",i] - AICmx)>threshold){
      warning(paste0("Wrong AIC for lambda = ", colnames(fit)[i]))
      return(works)
    }
    if(testIC && abs(fit["BIC",i] - BICmx)>threshold){
      warning(paste0("Wrong BIC for lambda = ", colnames(fit)[i]))
      return(works)
    }

    if(!is.null(cvModel)){

      cvModel <- OpenMx::omxSetParameters(cvModel,
                                          labels = parameterLabels,
                                          values = parameters[,i]
      )

      cvModel <- OpenMx::mxRun(cvModel, useOptimizer = F, silent = T)
      cvM2LL <- cvModel$fitfunction$result[[1]]
      if(abs(fit["cvM2LL",i] - cvM2LL)>threshold){
        warning("Wrong cvM2LL")
        return(works)
      }
    }
  }
  works <- TRUE

  print("FI Computation Works!")
  return(works)
}

