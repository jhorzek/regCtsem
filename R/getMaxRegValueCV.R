#' getMaxRegValueCV
#'
#'
#' computes an approximation of the lowest regValue which will set all regularized parameters to zero. This function is adapted from Murphy (2012) Machine learning: a probabilistic perspective. See p. 434 for more details.
#'
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param fullData Dataset for all samples combined
#' @param KalmanStartValues Optional starting values for the parameters when using Kalman filter
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso and adaptiveLasso for optimization = exact
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the unregularized parameter estimates.
#' @param standardizeDrift Boolean: Should Drift parameters be standardized automatically using T0VAR?
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @param cvFolds list of numeric vectors indicating which data row in fullData belongs to which cv sample
#' @param trainSets empty list of the same length as k. The data sets for training the models will be saved here
#' @param regIndicators Labels of the regularized parameters (e.g. drift_eta1_eta2)
#' @param differenceApprox Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the unregularized parameter estimates.
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param optimization which optimization procedure should be used. Possible are  "exact" or "approx".
#' @param differenceApprox Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#'
#' @author Jannik Orzek
#' @import OpenMx
#' @export

getMaxRegValueCV <- function(ctsemObject, mxObject, fullData, KalmanStartValues, regOn, regIndicators, penalty, adaptiveLassoWeights, standardizeDrift, k, cvFolds, trainSets, scaleLambdaWithN, objective, optimization, differenceApprox = "central"){

  maxRegValue <- 0

  for(foldNumber in 1:k){
    cat(paste0("Sample ", foldNumber, ": "))
    if(is.vector(fullData[-cvFolds[[foldNumber]],])){
      # if a single row is selected, R extracts this row as vector, not as matrix
      trainSets[[foldNumber]] <-  t(as.matrix(fullData[-cvFolds[[foldNumber]],]))
    }else{
      trainSets[[foldNumber]] <- fullData[-cvFolds[[foldNumber]],]
    }
    sampleSize <- length(unlist((cvFolds))) - length(cvFolds[[foldNumber]])

    if(tolower(objective) == "ml"){
      currentModel <- mxObject
      currentModel$data <- OpenMx::mxData(trainSets[[foldNumber]], type = "raw")
      currentModel <- mxRun(currentModel, silent = TRUE)

      currentAdaptiveLassoWeights <- getAdaptiveLassoWeights(mxObject = currentModel, penalty = penalty, adaptiveLassoWeights = adaptiveLassoWeights, standardizeDrift = standardizeDrift)

    }else if(tolower(objective) == "kalman"){
      # set input arguments
      currentModel <- createKalmanMultiSubjectModel(ctsemObject = ctsemObject,
                                                    dataset = trainSets[[foldNumber]],
                                                    useOptimizer = TRUE,
                                                    KalmanStartValues = KalmanStartValues)

      currentAdaptiveLassoWeights <- getAdaptiveLassoWeights(mxObject = currentModel,
                                                             penalty = penalty,
                                                             adaptiveLassoWeights = adaptiveLassoWeights,
                                                             standardizeDrift = standardizeDrift)
    }
    if(tolower(optimization) == "approx"){
      regIndicatorsString <- mxObject[[regOn]]$labels[regIndicators == 1]
      currentMaxRegValue <- getMaxRegValue(mxObject = currentModel,
                                           regIndicators = regIndicatorsString,
                                           differenceApprox = differenceApprox,
                                           adaptiveLassoWeights = currentAdaptiveLassoWeights)$maxRegValue
      sparseParameterMatrix <- NULL
    }else{
      parameterLabels <- names(OpenMx::omxGetParameters(currentModel))
      currentMaxRegValue <- getMaxRegValue(mxObject = currentModel,
                                           regIndicators = regIndicators,
                                           differenceApprox = differenceApprox,
                                           adaptiveLassoWeights = currentAdaptiveLassoWeights)
      currentSparseParameters <- currentMaxRegValue$sparseParameters
      currentMaxRegValue <- currentMaxRegValue$maxRegValue
      if(foldNumber == 1){
        sparseParameterMatrix <- matrix(NA, nrow = length(parameterLabels), ncol = k)
        rownames(sparseParameterMatrix) <- parameterLabels
        sparseParameterMatrix[parameterLabels,foldNumber] <- currentSparseParameters[parameterLabels]
      }else{
        sparseParameterMatrix[parameterLabels,foldNumber] <- currentSparseParameters[parameterLabels]
      }

    }

    if(scaleLambdaWithN){
      currentMaxRegValue <- currentMaxRegValue/sampleSize
    }

    maxRegValue <- max(maxRegValue,
                       currentMaxRegValue)
  }

  return(list("maxRegValue" = maxRegValue, "sparseParameterMatrix" = sparseParameterMatrix))
}



