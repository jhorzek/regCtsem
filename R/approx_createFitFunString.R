#' approx_createFitFunString
#'
#' Creates the new approx_createFitFunString for specified regularization
#'
#' @param penalty type of penalty. Currently supported are lasso and ridge
#' @param elastic_alpha PLACEHOLDER. Used in elastic net. Not yet implemented
#' @param elastic_gamma PLACEHOLDER. Used in elastic net. Not yet implemented
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentible one
#' @param silent Boolean indicating if print statements will be suppressed
#' @author Jannik Orzek
#' @keywords internal
#' @export
approx_createFitFunString = function(penalty, elastic_alpha = NULL, elastic_gamma = NULL, epsilon = .00001, silent = FALSE){
  ## unregularized fitfunction:
  fitFunString <- "submodel.fitfunction"

  ## add penalty
  if(tolower(penalty) == "lasso" | tolower(penalty) == "adaptivelasso"){
    addPenalty <- paste("numObs*(regValue*(omxSelectCols(rowVecDriftWeights,t(cvectorize(selectedDRIFTValues)))%*%omxSelectRows(cvectorize((submodel.DRIFT^2 + ",epsilon,")^(1/2)), cvectorize(selectedDRIFTValues))))", sep = "")
  }else if(penalty == "ridge"){
    addPenalty <- "numObs*(regValue*(omxSelectCols(rowVecDriftWeights,t(cvectorize(selectedDRIFTValues)))%*%omxSelectRows(cvectorize((submodel.DRIFT)^2), cvectorize(selectedDRIFTValues))))"

  }else if(penalty == "elastic"){
    stop("not implemented")
  }else{
    stop(paste("could not find pealty function ", penalty, ". Currently supported are lasso, ridge, and adaptiveLasso.", sep = ""))
  }
  fitFunString <- paste(fitFunString, addPenalty, sep = "+")
  return(fitFunString)
}



