#' iterateOverCVFolds
#'
#' computes results for automativ cross-validation
#'
#' @param argsIn list of parameters passed to regCtsem
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param optimization type of optimization. Either exact or approx
iterateOverCVFolds <- function(argsIn, objective = "ML", optimization){
  # create folds
  cvFoldsAndModels <- createCVFoldsAndModels(mxObject = argsIn$mxObject, dataset = argsIn$dataset, k = argsIn$k)
  cvFolds <- cvFoldsAndModels$cvFolds
  fullData <- cvFoldsAndModels$fullData
  testSets <- cvFoldsAndModels$testSets
  trainSets <- cvFoldsAndModels$trainSets
  cvModels <-  cvFoldsAndModels$cvModels

  # get regValues if regValues == "auto"
  if(any(argsIn$regValues == "auto")){
    maxRegValue <- getMaxRegValueCV(ctsemObject = argsIn$ctsemObject,
                                    mxObject = argsIn$mxObject,
                                    fullData = fullData,
                                    trainSets = trainSets,
                                    KalmanStartValues = argsIn$KalmanStartValues,
                                    regOn = argsIn$regOn,
                                    regIndicators = argsIn$regIndicators,
                                    penalty = argsIn$penalty,
                                    adaptiveLassoWeights = argsIn$adaptiveLassoWeights,
                                    standardizeDrift = argsIn$standardizeDrift,
                                    k = argsIn$k,
                                    cvFolds = cvFolds,
                                    scaleLambdaWithN = argsIn$scaleLambdaWithN,
                                    objective = objective,
                                    optimization = optimization,
                                    differenceApprox = ifelse(is.null(argsIn$differenceApprox), "central", argsIn$differenceApprox)
    )
    if(optimization == "exact"){
      sparseParameterMatrix <- maxRegValue$sparseParameterMatrix
    }
    maxRegValue <- maxRegValue$maxRegValue
    argsIn$regValues <- seq(0, maxRegValue, length.out = argsIn$regValuesAutoLength)
  }

  # create fit table
  fit <- createCVFitTable(k = argsIn$k, regValues = argsIn$regValues)

  for(foldNumber in 1:argsIn$k){
    if(is.vector(fullData[cvFolds[[foldNumber]],])){
      # if a single row is selected, R extracts this row as vector, not as matrix
      testSets[[foldNumber]] <-  t(as.matrix(fullData[cvFolds[[foldNumber]],]))
    }else{
      testSets[[foldNumber]] <-  fullData[cvFolds[[foldNumber]],]
    }
    if(is.vector(fullData[-cvFolds[[foldNumber]],])){
      # if a single row is selected, R extracts this row as vector, not as matrix
      trainSets[[foldNumber]] <-  t(as.matrix(fullData[-cvFolds[[foldNumber]],]))
    }else{
      trainSets[[foldNumber]] <- fullData[-cvFolds[[foldNumber]],]
    }

    if(tolower(objective) == "ml"){
      currentModel <- argsIn$mxObject
      currentModel$data <- OpenMx::mxData(trainSets[[foldNumber]], type = "raw")
      currentModel <- mxRun(currentModel, useOptimizer = F, silent = T)

      # set input arguments
      currentModelArgsIn <- argsIn
      currentModelArgsIn$ctsemObject$mxobj <- currentModel
      currentModelArgsIn$mxObject <- currentModel
      currentModelArgsIn$autoCV <- FALSE
      currentModelArgsIn$cvSample <- OpenMx::mxData(testSets[[foldNumber]], type = "raw")
      currentModelArgsIn$returnFitIndices <- FALSE

      if(optimization == "exact"){
        if(exists("sparseParameterMatrix")){
          sparseParameterLabels <- rownames(sparseParameterMatrix)
          currentModelArgsIn$sparseParameters <- sparseParameterMatrix[,foldNumber]
          names(currentModelArgsIn$sparseParameters) <- sparseParameterLabels
        }
        currentModelFit <- try(do.call(regCtsem::exact_regCtsem, currentModelArgsIn))
      }else if(optimization == "approx"){
        currentModelFit <- try(do.call(regCtsem::approx_regCtsem, currentModelArgsIn))
      }
    }else if(tolower(objective) == "kalman"){
      # set input arguments
      currentModelArgsIn <- argsIn
      currentModelArgsIn$autoCV <- FALSE
      currentModelArgsIn$dataset <- trainSets[[foldNumber]]
      currentModelArgsIn$cvSample <- testSets[[foldNumber]]
      currentModelArgsIn$returnFitIndices <- FALSE
      currentModelArgsIn$mxObject <- createKalmanMultiSubjectModel(ctsemObject = argsIn$ctsemObject,
                                                                   dataset = trainSets[[foldNumber]],
                                                                   useOptimizer = FALSE,
                                                                   KalmanStartValues = argsIn$KalmanStartValues)

      if(optimization == "exact"){
        if(exists("sparseParameterMatrix")){
          sparseParameterLabels <- rownames(sparseParameterMatrix)
          currentModelArgsIn$sparseParameters <- sparseParameterMatrix[,foldNumber]
          names(currentModelArgsIn$sparseParameters) <- sparseParameterLabels
        }
        currentModelFit <- try(do.call(exact_regCtsem, currentModelArgsIn))
      }else if(optimization == "approx"){
        currentModelFit <- try(do.call(approx_regCtsem, currentModelArgsIn))
      }
    }

    cvModels[[foldNumber]] <- currentModelFit

    if(!any(class(currentModelFit) == "try-error")){
      if(any(!(colnames(currentModelFit$fitAndParameters) %in% colnames(fit) ))){
        warning("Error when binding cv results. Returing results as global variable currentModelFitError")
        currentModelFitError <<- currentModelFit
        next
      }
      if(! "cvM2LL" %in% rownames(currentModelFit$fitAndParameters)){
        warning("Error when binding cv results. Returing results as global variable currentModelFitError")
        currentModelFitError <<- currentModelFit
        next
      }
      fit[foldNumber,colnames(currentModelFit$fitAndParameters)] <- currentModelFit$fitAndParameters["cvM2LL",colnames(currentModelFit$fitAndParameters)]
      cat("\n")
      print(paste("Finished CV", foldNumber, "of", argsIn$k))
    }
  }
  return(list("fit" = fit, "folds" = list("models" = cvModels, "members" = cvFolds)))
}



