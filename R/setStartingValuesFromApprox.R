#' setStartingValuesFromApprox
#'
#' set the parameter values of mxObject to the values in approx_regModel for the specified regValue
#'
#' @param approx_regModel fitted regCtsem with optimization = "approx" and without automatic cross-validation
#' @param mxObject Fitted object of class MxObject
#' @param regValue single regValue
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
setStartingValuesFromApprox <- function(approx_regModel, mxObject, regValue){
  parameterLabels <- names(OpenMx::omxGetParameters(mxObject))
  parameterValues <- approx_regModel$fitAndParameters[parameterLabels, as.character(regValue)]

  newModel <- OpenMx::omxSetParameters(model = mxObject, labels = parameterLabels, values = parameterValues)
  return(newModel)
}


