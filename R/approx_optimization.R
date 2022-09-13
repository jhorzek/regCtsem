#### The following functions implement approximate optimization in OpenMx


#### Main functions ####

#' approx_iterateOverLambdas
#'
#' loops over lambdas if optimization = "approx"
#'
#' NOTE: Function located in file approx_optimization.R
#'
#' @param cpptsemObject Fitted object of type cpptsem
#' @param dataset wide data set
#' @param sampleSize sample size
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso for optimization = exact
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param targetVector named vector with values towards which the parameters are regularized (Standard is regularization towards zero)
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param BICWithNAndT Boolean: TRUE = Use N and T in the formula for the BIC (-2log L + log(N+T)*k, where k is the number of parameters in the model). FALSE = Use both N in the formula for the BIC (-2log L + log(N))
#' @param Tpoints Number of time points (used for BICWithNAndT)
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
#' @param zeroThresh threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
#' @param controlApproxOptimizer settings passed to optimx or Rsolnp
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#'
#' @author Jannik Orzek
#' @import ctsemOMX
#' @keywords internal
approx_iterateOverLambdas <- function(  # model
  cpptsemObject,# = NULL,
  dataset,# = NULL,
  sampleSize,# = NULL,
  # penalty settings
  regIndicators,
  lambdas,
  penalty,# = "lasso",
  adaptiveLassoWeights,# = NULL,
  targetVector,
  # fit settings
  returnFitIndices,# = TRUE,
  BICWithNAndT,# = TRUE,
  Tpoints,# = NULL,
  # optimization settings
  objective,# = "ML",
  epsilon,# = .001,
  zeroThresh,# = .001,
  controlApproxOptimizer,
  # additional settings
  scaleLambdaWithN,
  verbose# = 0
  ){

  parameterLabels <- names(cpptsemObject$getParameterValues())
  initialParameterValues <- cpptsemObject$getParameterValues()
  startingValues <- initialParameterValues

  fitLabels <- c("regM2LL", "m2LL")
  if(returnFitIndices){
    fitLabels <- c(fitLabels, "AIC", "BIC", "estimatedParameters")
  }

  fitAndParameters <- matrix(NA, nrow = length(parameterLabels)+length(fitLabels), ncol = length(lambdas))
  rownames(fitAndParameters) <- c(fitLabels, parameterLabels)
  colnames(fitAndParameters) <- lambdas

  pbar <- utils::txtProgressBar(min = 0, max = length(lambdas), initial = 0, style = 3)

  for(iteration in 1:length(lambdas)){

    ## if (adaptive) lasso
    # create model
    if(penalty != "ridge"){
      if(controlApproxOptimizer$package == "optimx"){
      optimized <- try(approx_cpptsemOptimx(cpptsemmodel = cpptsemObject,
                                            regM2LLCpptsem = ifelse(tolower(objective) == "ml",
                                                                    approx_RAMRegM2LLCpptsem,
                                                                    approx_KalmanRegM2LLCpptsem),
                                            gradCpptsem = approx_gradCpptsem,
                                            startingValues = startingValues,
                                            adaptiveLassoWeights = adaptiveLassoWeights,
                                            N = ifelse(scaleLambdaWithN, sampleSize, 1), lambda = lambdas[iteration],
                                            regIndicators = regIndicators,
                                            targetVector = targetVector,
                                            epsilon = epsilon,
                                            objective = objective,
                                            testGradients = TRUE,
                                            controlOptimx =  controlApproxOptimizer,
                                            silent = FALSE), silent = TRUE)
      }else if(tolower(controlApproxOptimizer$package) == "rsolnp"){
        optimized <- try(approx_cpptsemRsolnp(cpptsemmodel = cpptsemObject,
                                              regM2LLCpptsem = ifelse(tolower(objective) == "ml",
                                                                      approx_RAMRegM2LLCpptsem,
                                                                      approx_KalmanRegM2LLCpptsem),
                                              gradCpptsem = approx_gradCpptsem,
                                              startingValues = startingValues,
                                              adaptiveLassoWeights = adaptiveLassoWeights,
                                              N = ifelse(scaleLambdaWithN, sampleSize, 1), lambda = lambdas[iteration],
                                              regIndicators = regIndicators,
                                              targetVector = targetVector,
                                              epsilon = epsilon,
                                              objective = objective,
                                              testGradients = TRUE,
                                              controlRsolnp =  controlApproxOptimizer,
                                              silent = FALSE), silent = TRUE)
      }else{
        stop("No package in controlApproxOptimizer")
      }
    }else{
      if(controlApproxOptimizer$package == "optimx"){
      optimized <- try(approx_cpptsemOptimx(cpptsemmodel = cpptsemObject,
                                            regM2LLCpptsem = ifelse(tolower(objective) == "ml",
                                                                    ridgeRAMRegM2LLCpptsem,
                                                                    ridgeKalmanRegM2LLCpptsem),
                                            gradCpptsem = approx_gradCpptsem,
                                            startingValues = startingValues,
                                            adaptiveLassoWeights = adaptiveLassoWeights,
                                            N = ifelse(scaleLambdaWithN, sampleSize, 1), lambda = lambdas[iteration],
                                            regIndicators = regIndicators,
                                            targetVector = targetVector,
                                            epsilon = epsilon,
                                            objective = objective,
                                            testGradients = TRUE,
                                            controlOptimx =  controlApproxOptimizer,
                                            silent = FALSE), silent = TRUE)
      }else if(tolower(controlApproxOptimizer$package) == "rsolnp"){
        optimized <- try(approx_cpptsemRsolnp(cpptsemmodel = cpptsemObject,
                                              regM2LLCpptsem = ifelse(tolower(objective) == "ml",
                                                                      ridgeRAMRegM2LLCpptsem,
                                                                      ridgeKalmanRegM2LLCpptsem),
                                              gradCpptsem = approx_gradCpptsem,
                                              startingValues = startingValues,
                                              adaptiveLassoWeights = adaptiveLassoWeights,
                                              N = ifelse(scaleLambdaWithN, sampleSize, 1), lambda = lambdas[iteration],
                                              regIndicators = regIndicators,
                                              targetVector = targetVector,
                                              epsilon = epsilon,
                                              objective = objective,
                                              testGradients = TRUE,
                                              controlRsolnp =  controlApproxOptimizer,
                                              silent = FALSE), silent = TRUE)

      }else{
        stop("No package in controlApproxOptimizer")
      }
    }
    # check if model returned errors:
    if(any(class(optimized) == "try-error")){
      # reset starting values
      startingValues <- initialParameterValues
      next
    } # go to next iteration if errors were encountered

    ## if no errors were produced
    # save current values as new starting values
    startingValues <- optimized$parameters

    # extract parameter estimates
    fitAndParameters[parameterLabels, as.character(lambdas[iteration])] <- optimized$parameters[parameterLabels]
    fits <- approx_getFitIndices(m2LL = optimized$m2LL,
                                           regM2LL = optimized$regM2LL,
                                           lambda = lambdas[iteration],
                                           parameterValues = optimized$parameters,
                                           targetVector = targetVector,
                                           sampleSize = ifelse(BICWithNAndT,Tpoints*sampleSize, sampleSize),
                                           # penalty settings
                                           regIndicators = regIndicators,
                                           # fit settings
                                           returnFitIndices = returnFitIndices,
                                           # optimization settings
                                           zeroThresh = zeroThresh)
    # extract fit
    fitAndParameters[fitLabels,as.character(lambdas[iteration])] <- fits[fitLabels,]

    # set progress bar for
    utils::setTxtProgressBar(pbar,iteration)
  }
  return(fitAndParameters)

}


#### Fit ####


#' approx_getFitIndices
#'
#' computes fit indices
#'
#' NOTE: Function located in file approx_optimization.R
#'
#' @param m2LL -2 log Likelihood
#' @param regM2LL regularized -2 log Likelihood
#' @param lambda tuning parameter lambda
#' @param parameterValues current parameter values
#' @param targetVector named vector with values towards which the parameters are regularized (Standard is regularization towards zero)
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param zeroThresh threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
#' @param sampleSize sample size
#' @author Jannik Orzek
#' @keywords internal
approx_getFitIndices <- function(m2LL,
                                 regM2LL,
                                 lambda,
                                 parameterValues,
                                 targetVector,
                                 sampleSize,
                                 # penalty settings
                                 regIndicators,
                                 # fit settings
                                 returnFitIndices,# = TRUE,
                                 # optimization settings
                                 zeroThresh # = .001
                                 ){
  fitLabels <- c("regM2LL", "m2LL")
  if(returnFitIndices){
    fitLabels <- c(fitLabels, "AIC", "BIC", "estimatedParameters")
  }

  fitTable <- matrix(NA, nrow = length(fitLabels), ncol = 1)
  colnames(fitTable) <- lambda
  rownames(fitTable) <- fitLabels

  fitTable["regM2LL",] <- regM2LL
  fitTable["m2LL",] <- m2LL

  if(returnFitIndices){
    if(lambda > 0){
      setZero <- abs(parameterValues[regIndicators] - targetVector[regIndicators]) <= zeroThresh
    }else{
      setZero <- 0
    }
    estimatedParameters <- length(parameterValues)-sum(setZero)

    fitTable["AIC", 1] <- m2LL + 2*estimatedParameters
    fitTable["BIC", 1] <- m2LL + log(sampleSize)*estimatedParameters
    fitTable["estimatedParameters", 1] <- estimatedParameters
  }

  return(fitTable)

}



