#' approx_cpptsemSolnp
#'
#' creates an approximate solution to regularized ctsem using solnp
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
#' @import cpptsem Rsolnp
#' @export
approx_cpptsemSolnp <- function(cpptsemmodel,
                                regM2LLCpptsem,
                                gradCpptsem,
                                startingValues,
                                adaptiveLassoWeights,
                                N, lambda,
                                regIndicators,
                                epsilon,
                                maxit,
                                objective,
                                failureReturns = 1e24,
                                testGradients){
  CpptsemFit <- Rsolnp::solnp(pars = startingValues,
          fun = regM2LLCpptsem,
          eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL,
          ineqUB = NULL, LB = NULL, UB = NULL, control = list(trace = 0, outer.iter = maxit),
          cpptsemmodel, adaptiveLassoWeights, N, lambda, regIndicators, epsilon, objective, failureReturns)

  if(testGradients){
    grad <- try(gradCpptsem(parameters = CpptsemFit$pars,
                            cpptsemmodel = cpptsemmodel,
                            adaptiveLassoWeights = adaptiveLassoWeights,
                            N =  N,lambda =  lambda, regIndicators = regIndicators,
                            epsilon = epsilon, objective =  objective,
                            failureReturns =  failureReturns))
    if(any(class(grad) == "try-error") || anyNA(grad)){
      stop("NA in gradients")
    }
  }

  return(list("parameters" = CpptsemFit$pars,
              "regM2LL" = CpptsemFit$values[length(CpptsemFit$values)]))
}



