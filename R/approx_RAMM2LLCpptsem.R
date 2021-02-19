#' approx_KalmanM2LLCpptsem
#'
#' computes gradients for an approximate optimization of regularized ctsem based on cpptsem with Full Information Maximum Likelihood objective
#' @param parameters paramneter values
#' @param cpptsemmodel model from cpptsem
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @author Jannik Orzek
#' @export
approx_RAMM2LLCpptsem <- function(parameters, cpptsemmodel, failureReturns){
  cpptsemmodel$setParameterValues(parameters, names(parameters))
  # catching all errors from cpptsemmodel
  # when parameter values are impossible
  invisible(capture.output(RAM <- try(cpptsemmodel$computeRAM(),
                                      silent = TRUE),
                           type = "message"))
  invisible(capture.output(FIT <- try(cpptsemmodel$fitRAM(),
                                      silent = TRUE),
                           type = "message"))
  if(class(RAM) == "try-error" | class(FIT) == "try-error"){
    return(failureReturns)
  }
  m2LL <- cpptsemmodel$m2LL
  if(is.na(m2LL) | is.infinite(m2LL)){
    return(failureReturns)
  }
  return(m2LL)
}




