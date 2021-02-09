#' approx_cpptsemOptim
#'
#' creates an approximate solution to regularized ctsem using optim with BFGS
#' @param cpptsemmodel model returned from cpptsem
#' @param regM2LLCpptsem regularized fitting function
#' @param gradCpptsem function for computing the gradients of regM2LLCpptsem
#' @param startingValues starting values for optimization
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
approx_cpptsemOptim <- function(cpptsemmodel,
                                regM2LLCpptsem,
                                gradCpptsem,
                                startingValues,
                                adaptiveLassoWeights,
                                N, lambda,
                                regIndicators,
                                epsilon,
                                maxit,
                                objective,
                                failureReturns = NA,
                                testGradients){
  CpptsemFit <- optim(par = startingValues,
                      fn = regM2LLCpptsem,
                      gr = gradCpptsem,
                      cpptsemmodel, adaptiveLassoWeights, N, lambda, regIndicators, epsilon, objective, failureReturns,
                      method = "BFGS",
                      control = list(maxit = maxit))

  if(testGradients){
    grad <- try(gradCpptsem(parameters = CpptsemFit$par,
                            cpptsemmodel = cpptsemmodel,
                            adaptiveLassoWeights = adaptiveLassoWeights,
                            N =  N,lambda =  lambda, regIndicators = regIndicators,
                            epsilon = epsilon, objective =  objective,
                            failureReturns =  failureReturns))
    if(any(class(grad) == "try-error") || anyNA(grad)){
      stop("NA in gradients")
    }
  }

  return(list("parameters" = CpptsemFit$par,
              "regM2LL" = CpptsemFit$value))
}



