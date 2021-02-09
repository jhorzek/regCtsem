#' approx_RAMRegM2LLCpptsem
#'
#' approximates the regularized likelihood function using cpptsem and Full Information Maximum Likelihood
#' @param parameters parameter values
#' @param cpptsemmodel model returned from cpptsem
#' @param adaptiveLassoWeights vector with weights of the adaptive lasso
#' @param N sample size
#' @param lambda tuning parameter lambda
#' @param regIndicators string vector with names of regularized parameters
#' @param epsilon tuning parameter for epsL1 approximation
#' @param objective ML or Kalman
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @author Jannik Orzek
#' @import cpptsem
#' @export
approx_RAMRegM2LLCpptsem <- function(parameters, cpptsemmodel, adaptiveLassoWeights, N, lambda, regIndicators, epsilon, objective, failureReturns){
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
  regM2LL <- m2LL + sum(N*lambda*adaptiveLassoWeights[regIndicators] * abs(parameters[regIndicators]))
  if(is.na(regM2LL) | is.infinite(regM2LL)){
    return(failureReturns)
  }
  return(regM2LL)
}



