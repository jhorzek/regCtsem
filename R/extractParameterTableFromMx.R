#' extractParameterTableFromMx
#'
#' wrapper for the omxLocateParameters function from OpenMx
#' returns a data frame with parameter labels, location and values
#'
#' @param mxObject Object of type mxObject
#' @import OpenMx
extractParameterTableFromMx <- function(mxObject){
  parameterTable <- OpenMx::omxLocateParameters(mxObject)
  return(parameterTable)
}

