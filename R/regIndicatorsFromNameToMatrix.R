#' regIndicatorsFromNameToMatrix
#'
#' transforms string indicated regIndicators to matrix indicated regIndicators
#'
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
regIndicatorsFromNameToMatrix <- function(mxObject, regOn, regIndicators){
  regIndicatorsMatrix <- matrix(as.numeric(c(mxObject[[regOn]]$labels %in% regIndicators)),
                                ncol = ncol(mxObject[[regOn]]$labels),
                                nrow = nrow(mxObject[[regOn]]$labels))
  return(regIndicatorsMatrix)
}




