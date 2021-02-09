#' approx_gradCpptsem
#'
#' computes gradients for an approximate optimization of regularized ctsem based on cpptsem
#' @param parameters parameter values
#' @param adaptiveLassoWeights vector with weights of the adaptive lasso
#' @param N sample size
#' @param lambda tuning parameter lambda
#' @param regIndicators string vector with names of regularized parameters
#' @param epsilon tuning parameter for epsL1 approximation
#' @param maxit maximal number of iterations
#' @param objective ML or Kalman
#' @param failureReturns value which is returned if regM2LLCpptsem or gradCpptsem fails
#' @param testGradients should be tested if the final parameters result in NA gradients?
#' @author Jannik Orzek
#' @import cpptsem
#' @export
approx_gradCpptsem <- function(parameters, cpptsemmodel, adaptiveLassoWeights, N, lambda, regIndicators, epsilon, objective, failureReturns){
  invisible(capture.output(grad <- try(exact_getCppGradients(cppmodel = cpptsemmodel, objective = objective),
                                       silent = TRUE),
                           type = "message"))
  if(class(grad) == "try-error" || anyNA(grad)){
    ret <- rep(failureReturns, length = length(parameters))
    names(ret) <- names(parameters)
    return(ret)
  }
  grad[regIndicators] <- grad[regIndicators] +
    N*lambda*adaptiveLassoWeights[regIndicators] * parameters[regIndicators]/sqrt(parameters[regIndicators]^2+epsilon)
  return(grad[names(parameters)])
}



