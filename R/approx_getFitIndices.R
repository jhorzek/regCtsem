#' approx_getFitIndices
#'
#' computes fit indices
#'
#' @param fittedRegModel fitted regCtsem object
#' @param ctsemObject Only for Kalman: Model to set up the Kalman object
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param cvSample cross-validation sample. Has to be of type mxData
#' @param zeroThresh threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
#' @param sampleSize sample size
#' @param objective Kalman or ML
#' @author Jannik Orzek
#' @export
approx_getFitIndices <- function(fittedRegModel,
                                 ctsemObject = NULL,
                                 sampleSize,
                                 # penalty settings
                                 regOn = "DRIFT",
                                 regIndicators,
                                 # fit settings
                                 returnFitIndices = TRUE,
                                 cvSample = NULL,
                                 # optimization settings
                                 zeroThresh = .001,
                                 objective = "ML"){
  fitLabels <- c("regM2LL", "m2LL")
  if(returnFitIndices){
    fitLabels <- c(fitLabels, "AIC", "BIC", "estimatedParameters")
  }
  if(!is.null(cvSample)){
    fitLabels <- c(fitLabels, "cvM2LL")
  }

  fitTable <- matrix(NA, nrow = length(fitLabels), ncol = 1)
  colnames(fitTable) <- fittedRegModel$regValue$values
  rownames(fitTable) <- fitLabels

  fitTable["regM2LL",] <- fittedRegModel$fitfunction$result[[1]]
  fitTable["m2LL",] <- fittedRegModel$submodel$fitfunction$result[[1]]

  if(returnFitIndices){
    currentParameters <- OpenMx::omxGetParameters(fittedRegModel)

    driftLabels <- fittedRegModel$submodel[[regOn]]$labels

    regularizedDrifts <- driftLabels[regIndicators == 1]
    setZero <- abs(currentParameters[regularizedDrifts]) <= zeroThresh
    estimatedParameters <- length(currentParameters)-sum(setZero)

    fitTable["AIC", as.character(fittedRegModel$regValue$values)] <- fitTable["m2LL",] + 2*estimatedParameters
    fitTable["BIC", as.character(fittedRegModel$regValue$values)] <- fitTable["m2LL",] + log(sampleSize)*estimatedParameters
    fitTable["estimatedParameters", as.character(fittedRegModel$regValue$values)] <- estimatedParameters
  }

  return(fitTable)

}


