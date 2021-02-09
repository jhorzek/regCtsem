#' createCVFitTable
#'
#' sets up a table for CV fit values
#'
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @param regValues vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @author Jannik Orzek
#' @export
createCVFitTable <- function(k, regValues){
  CVFitTable <- matrix(NA, nrow = k+1, ncol = length(regValues))
  rownames(CVFitTable) <- c(paste0("fold", 1:k), "mean CV fit")
  colnames(CVFitTable) <- regValues
  return(CVFitTable)
}



