#' getFlatStdizer
#'
#' returns the standardizer for the standardized drift parameters if standardizeDrift = TRUE. Computes T0SD_predictor/T0SD_dependent
#' @param T0VAR matrix with T0VAR values
#' @param driftLabels vector with drift names
#' @export
getFlatStdizer <- function(T0VAR, driftLabels){
  stdizer <- matrix(1, nrow = nrow(T0VAR), ncol = ncol(T0VAR))
  stdizer <- stdizer%*%diag(sqrt(diag(T0VAR))) # times predictor sd
  stdizer <- t(t(stdizer)%*%diag(sqrt(diag(T0VAR))^(-1))) # divided by dependent sd
  flatStdizer <- OpenMx::cvectorize(stdizer) # flatten
  rownames(flatStdizer) <- as.vector(driftLabels)
  return(flatStdizer)
}


