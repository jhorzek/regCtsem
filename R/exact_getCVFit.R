#' exact_getCVFit
#'
#' computes cross-validation fit for exact optimization
#'
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param ctsemObject if objective = "ML": Fitted object of class ctsem. If you want to use objective = "Kalman", pass an object of type ctsemInit from ctModel
#' @param mxObject Fitted object of class MxObject
#' @param parameterLabels labels of optimized parameters
#' @param parameterValuesTable table with parameter values+
#' @param regValues vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param cvSample cross-validation sample. Has to be of type mxData
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
exact_getCVFit <- function(objective, ctsemObject, mxObject, parameterLabels,
                           parameterValuesTable, regValues,
                           cvSample){
  fits <- c("cvM2LL")
  fitTable <- matrix(NA, nrow = length(fits), ncol = length(regValues), dimnames = list(fits, regValues))

  if(tolower(objective) == "kalman"){
    # create Model
    cvModel <- createKalmanMultiSubjectModel(ctsemObject = ctsemObject, dataset = cvSample, useOptimizer = FALSE, silent = TRUE)
  }else{
    # create Model
    cvModel <- mxObject
    # replace sample
    cvModel$data <- cvSample
  }


  for(regValue in 1:length(regValues)){
    if(any(is.na(parameterValuesTable[,regValue]))){
      next
    }

    # set parameters to iteration parameters
    cvModel <- try(OpenMx::omxSetParameters(model = cvModel, labels = parameterLabels, values = parameterValuesTable[,regValue]))

    if(!any(class(cvModel) == "try-error")){
      # compute -2LL
      fitCvModel <- OpenMx::mxRun(cvModel, useOptimizer = FALSE, silent = TRUE)
      fitTable["cvM2LL",regValue] <- fitCvModel$fitfunction$result[[1]]
    }

  }
  return(fitTable)
}



