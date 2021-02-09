#' approx_getCVFit
#'
#' computes cv fit for approximate optimization
#' @param ctsemObject Fitted object of class ctsemOMX
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param cvSample only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param fitAndParameters table with fit results and parameter values
#' @param parameterLabels vector with parameter labels
#' @param regValues vector with regValues
#' @author Jannik Orzek
#' @export
approx_getCVFit <- function(mxObject, ctsemObject, cvSample, objective, fitAndParameters, parameterLabels, regValues){
  if(tolower(objective) == "kalman"){
    # create Model
    cvModel <- createKalmanMultiSubjectModel(ctsemObject = ctsemObject, dataset = cvSample, useOptimizer = FALSE, silent = TRUE)

  }else{
    # create Model
    cvModel <- mxObject
    # replace sample
    cvModel$data <- cvSample
  }
  for(iteration in 1:length(regValues)){

    currentParameters <- fitAndParameters[parameterLabels,iteration]
    # set parameters to iteration parameters
    cvModel <- try(OpenMx::omxSetParameters(model = cvModel, labels = parameterLabels, values = currentParameters))

    # compute -2LL
    if(any(class(cvModel) == "try-error")){next}
    fitCvModel <- OpenMx::mxRun(cvModel, useOptimizer = F, silent = T)
    fitAndParameters["cvM2LL", iteration] <- fitCvModel$fitfunction$result[[1]]

  }

  return(fitAndParameters)
}



