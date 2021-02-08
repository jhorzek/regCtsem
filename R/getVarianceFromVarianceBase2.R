#' getVarianceFromVarianceBase2
#'
#' returns variance from varianceBase
#'
#' @param varianceBaseValues values of varianceBase
#' @import OpenMx
#' @author Jannik Orzek
#' @export
getVarianceFromVarianceBase2 <- function(varianceBaseValues){
  varianceCholValues <- OpenMx::vec2diag(exp(OpenMx::diag2vec(varianceBaseValues))) +
    varianceBaseValues - OpenMx::vec2diag(OpenMx::diag2vec(varianceBaseValues))
  varianceValues <- varianceCholValues %*% t(varianceCholValues)
  return(varianceValues)
}


