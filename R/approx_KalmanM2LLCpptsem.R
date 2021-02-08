#' approx_KalmanM2LLCpptsem
#'
#' computes gradients for an approximate optimization of regularized ctsem based on cpptsem with Kalman objective
#' @param parameters paramneter values
#' @param cpptsemmodel model from cpptsem
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @author Jannik Orzek
#' @import cpptsem
#' @export
approx_KalmanM2LLCpptsem <- function(parameters, cpptsemmodel, failureReturns){
  cpptsemmodel$setParameterValues(parameters, names(parameters))
  # catching all errors from cpptsemmodel
  # when parameter values are impossible
  invisible(capture.output(FIT <- try(cpptsemmodel$computeAndFitKalman(),
                                      silent = TRUE),
                           type = "message"))
  if(class(FIT) == "try-error"){
    return(failureReturns)
  }
  m2LL <- cpptsemmodel$m2LL
  if(is.na(m2LL) | is.infinite(m2LL)){
    return(failureReturns)
  }
  return(m2LL)
}


