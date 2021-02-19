#' exact_getFitIndices
#'
#' computes fit indices for optimization = "exact"
#'
#' @param mxObject Fitted object of class MxObject
#' @param parameterLabels labels of optimized parameters
#' @param fitAndParameters table with fit and parameter values
#' @param regValues vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param sampleSize sample size
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
exact_getFitIndices <- function(mxObject, parameterLabels, fitAndParameters, regValues, sampleSize){
  fits <- c("AIC", "BIC", "estimatedParameters")
  fitTable <- matrix(NA, nrow = length(fits), ncol = length(regValues), dimnames = list(fits, regValues))

  for(regValue in 1:length(regValues)){

    currentParameterValues <- fitAndParameters[parameterLabels,regValue]
    currentM2LL <- fitAndParameters["m2LL",regValue]

    if(any(is.na(currentParameterValues))){next}

    # set free = TRUE to free = FALSE for zeroed parameters
    currentParameterFree <- !(currentParameterValues == 0)
    fitTable["AIC", regValue] <- currentM2LL + 2*sum(currentParameterFree)
    fitTable["BIC", regValue] <- currentM2LL + log(sampleSize)*sum(currentParameterFree)
    fitTable["estimatedParameters", regValue] <- sum(currentParameterFree)


    #fitModel <- mxObject
    #fitModel <- OpenMx::omxSetParameters(model = fitModel, labels = parameterLabels, free = currentParameterFree, values = currentParameterValues)

    # run Model
    #fitFitModel <- OpenMx::mxRun(fitModel, useOptimizer = F, silent = T)

    #fitTable["AIC", regValue] <- AIC(fitFitModel)
    #fitTable["BIC", regValue] <- BIC(fitFitModel)
    #fitTable["estimatedParameters", regValue] <- length(OpenMx::omxGetParameters(fitFitModel))

  }

  return(fitTable)

}




