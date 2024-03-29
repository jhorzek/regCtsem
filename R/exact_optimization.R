### The following functions are used in the exact optimization. Note that the actual optimizers
### are defined in separate files (GIST and GLMNET)

#' exact_getCVFit
#'
#' computes cross-validation fit for exact optimization.
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param objective Kalman or mx
#' @param cvSampleCpptsemObject Fitted object of class cpptsem
#' @param parameterLabels labels of optimized parameters
#' @param parameterValuesTable table with parameter values+
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param cvSample cross-validation sample. Has to be of type mxData
#' @author Jannik Orzek
#' @import ctsemOMX
#' @keywords internal
exact_getCVFit <- function(objective, cvSampleCpptsemObject, parameterLabels,
                           parameterValuesTable, lambdas,
                           cvSample){
  fits <- c("cvM2LL")
  fitTable <- matrix(NA, nrow = length(fits), ncol = length(lambdas), dimnames = list(fits, lambdas))

  for(lambda in 1:length(lambdas)){
    if(any(is.na(parameterValuesTable[,lambda]))){
      next
    }

    # set parameters to iteration parameters
    cvModelFit <- try(fitCpptsem(parameterValues = parameterValuesTable[parameterLabels,lambda], cpptsemObject = cvSampleCpptsemObject, objective = objective, failureReturns = NA))

    if(!any(class(cvModelFit) == "try-error")){
      fitTable["cvM2LL",lambda] <- cvModelFit
    }

  }
  return(fitTable)
}


#' exact_getFitIndices
#'
#' computes fit indices for optimization = "exact"
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param parameterLabels labels of optimized parameters
#' @param fitAndParameters table with fit and parameter values
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param sampleSize sample size
#' @author Jannik Orzek
#' @import ctsemOMX
#' @keywords internal
exact_getFitIndices <- function(parameterLabels, fitAndParameters, lambdas, sampleSize){
  fits <- c("AIC", "BIC", "estimatedParameters")
  fitTable <- matrix(NA, nrow = length(fits), ncol = length(lambdas), dimnames = list(fits, lambdas))

  for(lambda in 1:length(lambdas)){

    currentParameterValues <- fitAndParameters[parameterLabels,lambda]
    currentM2LL <- fitAndParameters["m2LL",lambda]

    if(any(is.na(currentParameterValues))){next}

    # set free = TRUE to free = FALSE for zeroed parameters
    currentParameterFree <- !(currentParameterValues == 0)
    fitTable["AIC", lambda] <- currentM2LL + 2*sum(currentParameterFree)
    fitTable["BIC", lambda] <- currentM2LL + log(sampleSize)*sum(currentParameterFree)
    fitTable["estimatedParameters", lambda] <- sum(currentParameterFree)
  }

  return(fitTable)

}


#' exact_getFitIndicesWithTarget
#'
#' computes fit indices for optimization = "exact"
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param parameterLabels labels of optimized parameters
#' @param regIndicators Vector with names of regularized parameters
#' @param fitAndParameters table with fit and parameter values
#' @param targetVector named vector with values towards which the parameters are regularized
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param sampleSize sample size
#' @author Jannik Orzek
#' @import ctsemOMX
#' @keywords internal
exact_getFitIndicesWithTarget <- function(parameterLabels, regIndicators, fitAndParameters, targetVector, lambdas, sampleSize){
  fits <- c("AIC", "BIC", "estimatedParameters")
  fitTable <- matrix(NA, nrow = length(fits), ncol = length(lambdas), dimnames = list(fits, lambdas))

  for(lambda in 1:length(lambdas)){

    currentParameterValues <- fitAndParameters[parameterLabels,lambda]
    currentM2LL <- fitAndParameters["m2LL",lambda]

    if(any(is.na(currentParameterValues))){next}

    # set free = TRUE to free = FALSE for zeroed parameters
    currentParameterFree <- length(currentParameterValues) - sum(currentParameterValues[regIndicators] == targetVector[regIndicators])
    fitTable["AIC", lambda] <- currentM2LL + 2*sum(currentParameterFree)
    fitTable["BIC", lambda] <- currentM2LL + log(sampleSize)*sum(currentParameterFree)
    fitTable["estimatedParameters", lambda] <- sum(currentParameterFree)

  }

  return(fitTable)

}



#' exact_getPenaltyValue
#'
#' computes sum(lambda*abs(regularized Values))
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param lambda Penalty value
#' @param theta parameter values
#' @param regIndicators Names of regularized parameters
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @keywords internal
exact_getPenaltyValue <- function(lambda, theta, regIndicators, adaptiveLassoWeights){

  regVal <- 0
  if(!is.vector(theta)){
    thetaNames <- rownames(theta)
    # iterate over theta
    for(th in thetaNames){
      # if theta is regularized
      if(th %in% regIndicators){
        regVal <- regVal + lambda*abs(theta[th,])*ifelse(is.null(adaptiveLassoWeights),1,adaptiveLassoWeights[th])
      }
    }
    names(regVal) <- NULL
    return(regVal)
  }

  thetaNames <- names(theta)
  # iterate over theta
  for(th in thetaNames){
    # if theta is regularized
    if(th %in% regIndicators){
      regVal <- regVal + lambda*abs(theta[th])*ifelse(is.null(adaptiveLassoWeights),1,adaptiveLassoWeights[th])
    }
  }
  names(regVal) <- NULL
  return(regVal)
}

#' exact_getPenaltyValueWithTarget
#'
#' computes sum(lambda*abs(regularized Values))
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param lambda Penalty value
#' @param theta parameter values
#' @param regIndicators Names of regularized parameters
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param targetVector named vector with values towards which the parameters are regularized
#' @keywords internal
exact_getPenaltyValueWithTarget <- function(lambda, theta, regIndicators, targetVector, adaptiveLassoWeights){

  regVal <- 0
  if(!is.vector(theta)){
    thetaNames <- rownames(theta)
    # iterate over theta
    for(th in thetaNames){
      # if theta is regularized
      if(th %in% regIndicators){
        regVal <- regVal + lambda*abs(theta[th,] - targetVector[th])*ifelse(is.null(adaptiveLassoWeights),1,adaptiveLassoWeights[th])
      }
    }
    names(regVal) <- NULL
    return(regVal)
  }

  thetaNames <- names(theta)
  # iterate over theta
  for(th in thetaNames){
    # if theta is regularized
    if(th %in% regIndicators){
      regVal <- regVal + lambda*abs(theta[th] - targetVector[th])*ifelse(is.null(adaptiveLassoWeights),1,adaptiveLassoWeights[th])
    }
  }
  names(regVal) <- NULL
  return(regVal)
}


#' exact_getSubgradients
#'
#' Computes the subgradients of a lasso penalized likelihood f(theta) = L(theta)+lambda*p(theta)
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param theta vector with named parameters
#' @param jacobian derivative of L(theta)
#' @param regIndicators names of regularized parameters
#' @param lambda lambda value
#' @param adaptiveLassoWeightsMatrix matrix with adaptive lasso weights
#' @keywords internal
exact_getSubgradients <- function(theta, jacobian, regIndicators, lambda, adaptiveLassoWeightsMatrix){

  # first part: derivative of Likelihood
  subgradient <- jacobian

  # second part: derivative of penalty term
  ## Note: if adaptiveLassoWeightsMatrix*parameter != 0, the derivative is: lambda*adaptiveLassoWeight*sign(parameter)
  ## if adaptiveLassoWeightsMatrix*parameter == 0, the derivative is in [-lambda*adaptiveLassoWeight, +lambda*adaptiveLassoWeight]
  for(regularizedParameterLabel in regIndicators){
    absoluteValueOf <- theta[regularizedParameterLabel,]
    if(absoluteValueOf != 0){
      penaltyGradient <- adaptiveLassoWeightsMatrix[regularizedParameterLabel,regularizedParameterLabel]*lambda*sign(absoluteValueOf)
      subgradient[regularizedParameterLabel,] <- subgradient[regularizedParameterLabel,] + penaltyGradient
    }else{
      penaltyGradient <- c(-adaptiveLassoWeightsMatrix[regularizedParameterLabel,regularizedParameterLabel]*lambda, adaptiveLassoWeightsMatrix[regularizedParameterLabel,regularizedParameterLabel]*lambda)
      # check if likelihood gradient is within interval:
      setZero <- (subgradient[regularizedParameterLabel,] > penaltyGradient[1]) && (subgradient[regularizedParameterLabel,] < penaltyGradient[2])
      subgradient[regularizedParameterLabel,] <- ifelse(setZero,
                                                        0,
                                                        sign(subgradient[regularizedParameterLabel,])*(abs(subgradient[regularizedParameterLabel,]) -
                                                                                                         adaptiveLassoWeightsMatrix[regularizedParameterLabel,
                                                                                                                                    regularizedParameterLabel]*lambda))
    }
  }
  return(subgradient)
}

#' exact_tryStartingValues
#'
#' tries different starting values to find a good starting point for exact optimization
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param startingValues start values for optimization
#' @param returnAs "vector" or "matrix"
#' @param approxFirst boolean: should the solution be approximated
#' @param numStart number of starting values
#' @param controlApproxOptimizer settings passed to optimx or Rsolnp
#' @param lambda tuning parameter value
#' @param cpptsemObject object of type cpptsem
#' @param regIndicators Labels for the regularized parameters (e.g. drift_eta1_eta2)
#' @param targetVector named vector with values towards which the parameters are regularized
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param objective which objective is used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param sparseParameters labeled vector with parameter estimates of the most sparse model.
#' @author Jannik Orzek
#' @import ctsemOMX
#' @keywords internal
exact_tryStartingValues <- function(startingValues,
                                    returnAs,
                                    approxFirst,
                                    numStart,
                                    controlApproxOptimizer,
                                    lambda,
                                    cpptsemObject,
                                    regIndicators,
                                    targetVector,
                                    adaptiveLassoWeights,
                                    objective,
                                    sparseParameters){
  if(is.null(sparseParameters)){
    if(numStart > 0){
      warning("Can not try different starting values because sparseParameters is empty.")
    }
    invisible(utils::capture.output(newStartingValues <- tryApproxFirst(startingValues = startingValues,
                                                                 returnAs = returnAs,
                                                                 approxFirst = approxFirst,
                                                                 numStart = numStart,
                                                                 controlApproxOptimizer = controlApproxOptimizer,
                                                                 lambda = lambda,
                                                                 cpptsemObject = cpptsemObject,
                                                                 regIndicators = regIndicators,
                                                                 targetVector = targetVector,
                                                                 adaptiveLassoWeights = adaptiveLassoWeights,
                                                                 objective = objective,
                                                                 sparseParameters = sparseParameters), type = "message"))
    if(any(class(newStartingValues) == "try-error")){
      return(startingValues)
    }
    return(newStartingValues)
  }

  # generate multiple starting values
  weights <- seq(-5.2,5.6,length.out = numStart)
  weights <- c(1,1/(1+exp(weights)), 0)
  parameterNames <- names(startingValues)
  startValuesTable <- matrix(NA, nrow = length(startingValues), ncol = length(weights))
  regM2LLs <- rep(NA, length(weights))
  rownames(startValuesTable) <- parameterNames

  for(i in seq_len(length(weights))){
    # middle ground between current parameters and sparse parameters
    #weight <- (a^weights[i] - 1)/(a-1)
    weight <- weights[i]
    parameterValues_i <- weight*startingValues[parameterNames]+(1-weight)*sparseParameters[parameterNames]
    startValuesTable[parameterNames,i] <- parameterValues_i[parameterNames]
    # set parameters
    invisible(utils::capture.output(optimized <- try(tryApproxFirst(startingValues = parameterValues_i,
                                                             returnAs = "full",
                                                             approxFirst = approxFirst,
                                                             numStart = numStart,
                                                             controlApproxOptimizer = controlApproxOptimizer,
                                                             lambda = lambda,
                                                             cpptsemObject = cpptsemObject,
                                                             regIndicators = regIndicators,
                                                             targetVector = targetVector,
                                                             adaptiveLassoWeights = adaptiveLassoWeights,
                                                             objective = objective,
                                                             sparseParameters = sparseParameters),
                                              silent = TRUE), type = "message"))
    if(any(class(optimized) == "try-error")){
      next
    }

    regM2LLs[i] <- optimized$regM2LL
    startValuesTable[parameterNames,i] <- optimized$parameters[parameterNames]

  }

  if(all(is.na(regM2LLs))){stop()}
  startValues <- startValuesTable[parameterNames,which(regM2LLs == min(regM2LLs, na.rm = T))[1]]
  names(startValues) <- parameterNames
  return(startValues)
}

#' tryApproxFirst
#'
#' Approximates the solution to the lasso regularized fitting function using an established optimizer
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param startingValues start values for optimization
#' @param returnAs "vector" or "matrix"
#' @param approxFirst boolean: should the solution be approximated
#' @param numStart number of starting values
#' @param controlApproxOptimizer settings passed to optimx or Rsolnp
#' @param lambda tuning parameter value
#' @param cpptsemObject object of type cpptsem
#' @param regIndicators Labels for the regularized parameters (e.g. drift_eta1_eta2)
#' @param targetVector named vector with values towards which the parameters are regularized
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param objective which objective is used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param sparseParameters labeled vector with parameter estimates of the most sparse model.
#' @author Jannik Orzek
#' @import ctsemOMX
#' @keywords internal
tryApproxFirst <- function(startingValues, returnAs,
                           approxFirst,
                           numStart,
                           controlApproxOptimizer,
                           lambda,
                           cpptsemObject,
                           regIndicators,
                           targetVector,
                           adaptiveLassoWeights,
                           objective,
                           sparseParameters
){

  # define a target vector if none is provided
  if(is.null(targetVector)){
    targetVector <- rep(0, length(regIndicators))
    names(targetVector) <- regIndicators
  }

  if(approxFirst && !is.null(cpptsemObject)){
    if(!is.null(controlApproxOptimizer$nudgeVariancesLambda) && controlApproxOptimizer$nudgeVariancesLambda > 0){
      varianceLabels <- cpptsemObject$parameterTable$label[grepl("base", cpptsemObject$parameterTable$matrix) & cpptsemObject$parameterTable$row == cpptsemObject$parameterTable$col]
      if(!any(varianceLabels %in% regIndicators)){
        # only use nudging if variances are not already regularized
        regIndicators_nudged <- c(regIndicators, varianceLabels)
        targetVector_nudged <- c(targetVector, rep(controlApproxOptimizer$nudgeVariancesTarget, length(varianceLabels))) # exp(-1.609438) = .2
        names(targetVector_nudged) <- regIndicators_nudged
        adaptiveLassoWeights_nudged <- adaptiveLassoWeights
        adaptiveLassoWeights_nudged[varianceLabels] <- ifelse(lambda == 0, 0, controlApproxOptimizer$nudgeVariancesLambda/lambda)
      }else{
        regIndicators_nudged <- regIndicators
        targetVector_nudged <- targetVector
        adaptiveLassoWeights_nudged <- adaptiveLassoWeights
      }
    }else{
      regIndicators_nudged <- regIndicators
      targetVector_nudged <- targetVector
      adaptiveLassoWeights_nudged <- adaptiveLassoWeights
    }
    if(controlApproxOptimizer$package == "optimx"){
      optimized <- try(approx_cpptsemOptimx(cpptsemmodel = cpptsemObject,
                                            regM2LLCpptsem = ifelse(tolower(objective) == "ml",
                                                                    approx_RAMRegM2LLCpptsem,
                                                                    approx_KalmanRegM2LLCpptsem),
                                            gradCpptsem = approx_gradCpptsem,
                                            startingValues = startingValues,
                                            adaptiveLassoWeights = adaptiveLassoWeights_nudged,
                                            N = 1, # lambda is already adjusted
                                            lambda = lambda,
                                            regIndicators = regIndicators_nudged,
                                            targetVector = targetVector_nudged,
                                            epsilon = 10^(-8),
                                            objective = objective,
                                            testGradients = TRUE,
                                            controlOptimx = controlApproxOptimizer, silent = TRUE), silent = TRUE)
    }else if(tolower(controlApproxOptimizer$package) == "rsolnp"){
      optimized <- try(approx_cpptsemRsolnp(cpptsemmodel = cpptsemObject,
                                            regM2LLCpptsem = ifelse(tolower(objective) == "ml",
                                                                    approx_RAMRegM2LLCpptsem,
                                                                    approx_KalmanRegM2LLCpptsem),
                                            gradCpptsem = approx_gradCpptsem,
                                            startingValues = startingValues,
                                            adaptiveLassoWeights = adaptiveLassoWeights_nudged,
                                            N = 1, # lambda is already adjusted
                                            lambda = lambda,
                                            regIndicators = regIndicators_nudged,
                                            targetVector = targetVector_nudged,
                                            epsilon = 10^(-8),
                                            objective = objective,
                                            testGradients = TRUE,
                                            controlRsolnp = controlApproxOptimizer,
                                            silent = TRUE), silent = TRUE)
    }else{
      stop("No package in controlApproxOptimizer")
    }
    if(!any(class(optimized) == "try-error") && returnAs == "full"){
      return(optimized)
    }
    if(any(class(optimized) == "try-error") && returnAs == "full"){
      stop("Error in tryApproxFirst")
    }
    if(!any(class(optimized) == "try-error")){
      startingValues <- optimized$parameters
    }
    if(returnAs == "vector"){
      return(startingValues)
    }else{
      startingValues <- as.matrix(startingValues)
      return(startingValues)
    }
  }

  stop("Error in tryApproxFirst")

}
