#' separateFitAndParameters
#'
#' separates fit results from parameter estimates
#'
#' @param regCtsemObject Object from exact_regCtsem or approx_regCtsem
#' @import OpenMx
#' @author Jannik Orzek
#' @export
separateFitAndParameters <- function(regCtsemObject){
  # get parameter labels:

  parameterLabels <- names(OpenMx::omxGetParameters(regCtsemObject$setup$mxObject))
  fitAndParametersLabels <- rownames(regCtsemObject$fitAndParameters)
  fitLabels <- fitAndParametersLabels[!(fitAndParametersLabels %in% parameterLabels)]

  fit <- regCtsemObject$fitAndParameters[fitLabels,]
  parameterEstimates <- regCtsemObject$fitAndParameters[parameterLabels,]

  if(!is.matrix(parameterEstimates)){
    parameterEstimates <- matrix(parameterEstimates, nrow = length(parameterLabels))
    rownames(parameterEstimates) <- parameterLabels
  }
  if(!is.matrix(fit)){
    fit <- matrix(fit, nrow = length(fitAndParametersLabels[fitLabels]))
    rownames(fit) <- fitLabels
  }
  return(list("fit" = fit, "parameterEstimates" = parameterEstimates))

}


