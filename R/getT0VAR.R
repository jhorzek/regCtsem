#' getT0VAR
#'
#' extracts the T0VAR from an mxObject
#' @param mxObject mxObject
#' @export
getT0VAR <- function(mxObject){

  # case 1: T0VAR not stationary
  if(!is.null(mxObject$T0VAR$result)){
    return(mxObject$T0VAR$result)
  }
  # case 2: T0VAR stationary
  stop("stop from getT0VAR: Automatic standardization not yet implemented for stationarity T0VAR")

}



