### The following functions are used in the exact optimization. Note that the actual optimizers
### are defined in separate files (GIST and GLMNET)

#' exact_getCVFit
#'
#' computes cross-validation fit for exact optimization.
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param objective Kalman or mx
#' @param ctsemObject Fitted object of class ctsemFit
#' @param mxObject Fitted object of class MxObject
#' @param parameterLabels labels of optimized parameters
#' @param parameterValuesTable table with parameter values+
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param cvSample cross-validation sample. Has to be of type mxData
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
exact_getCVFit <- function(objective, ctsemObject, mxObject, parameterLabels,
                           parameterValuesTable, lambdas,
                           cvSample){
  fits <- c("cvM2LL")
  fitTable <- matrix(NA, nrow = length(fits), ncol = length(lambdas), dimnames = list(fits, lambdas))

  if(tolower(objective) == "kalman"){
    suppressMessages(invisible(capture.output(ctsemTemp <- ctFit(dat = cvSample, ctmodelobj = ctsemObject$ctmodelobj, objective = "Kalman", useOptimizer = FALSE))))
    cvModel <- ctsemTemp$mxobj
  }else{
    # create Model
    cvModel <- mxObject
    # replace sample
    cvModel$data <- cvSample
  }

  for(lambda in 1:length(lambdas)){
    if(any(is.na(parameterValuesTable[,lambda]))){
      next
    }

    # set parameters to iteration parameters
    cvModel <- try(OpenMx::omxSetParameters(model = cvModel, labels = parameterLabels, values = parameterValuesTable[,lambda]))

    if(!any(class(cvModel) == "try-error")){
      # compute -2LL
      fitCvModel <- OpenMx::mxRun(cvModel, useOptimizer = FALSE, silent = TRUE)
      fitTable["cvM2LL",lambda] <- fitCvModel$fitfunction$result[[1]]
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
#' @param mxObject Fitted object of class MxObject
#' @param parameterLabels labels of optimized parameters
#' @param fitAndParameters table with fit and parameter values
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param sampleSize sample size
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
exact_getFitIndices <- function(mxObject, parameterLabels, fitAndParameters, lambdas, sampleSize){
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


    #fitModel <- mxObject
    #fitModel <- OpenMx::omxSetParameters(model = fitModel, labels = parameterLabels, free = currentParameterFree, values = currentParameterValues)

    # run Model
    #fitFitModel <- OpenMx::mxRun(fitModel, useOptimizer = F, silent = T)

    #fitTable["AIC", lambda] <- AIC(fitFitModel)
    #fitTable["BIC", lambda] <- BIC(fitFitModel)
    #fitTable["estimatedParameters", lambda] <- length(OpenMx::omxGetParameters(fitFitModel))

  }

  return(fitTable)

}


#' exact_getFitIndicesWithTarget
#'
#' computes fit indices for optimization = "exact"
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param mxObject Fitted object of class MxObject
#' @param parameterLabels labels of optimized parameters
#' @param fitAndParameters table with fit and parameter values
#' @param regIndicators Vector with names of regularized parameters
#' @param targetVector named vector with values towards which the parameters are regularized
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param sampleSize sample size
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
exact_getFitIndicesWithTarget <- function(mxObject, parameterLabels, regIndicators, fitAndParameters, targetVector, lambdas, sampleSize){
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


    #fitModel <- mxObject
    #fitModel <- OpenMx::omxSetParameters(model = fitModel, labels = parameterLabels, free = currentParameterFree, values = currentParameterValues)

    # run Model
    #fitFitModel <- OpenMx::mxRun(fitModel, useOptimizer = F, silent = T)

    #fitTable["AIC", lambda] <- AIC(fitFitModel)
    #fitTable["BIC", lambda] <- BIC(fitFitModel)
    #fitTable["estimatedParameters", lambda] <- length(OpenMx::omxGetParameters(fitFitModel))

  }

  return(fitTable)

}

#' exact_getFlatStdizer
#'
#' returns the standardizer for the standardized drift parameters if standardizeDrift = TRUE. Computes T0SD_predictor/T0SD_dependent
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param T0VAR matrix with T0VAR values
#' @param thetaNames vector with parameter names
#' @export
exact_getFlatStdizer <- function(T0VAR, thetaNames){
  stdizer <- matrix(1, nrow = nrow(T0VAR), ncol = ncol(T0VAR))
  stdizer <- stdizer%*%diag(sqrt(diag(T0VAR))) # times predictor sd
  stdizer <- t(t(stdizer)%*%diag(sqrt(diag(T0VAR))^(-1))) # divided by dependent sd
  flatStdizer <- OpenMx::cvectorize(stdizer) # flatten
  rownames(flatStdizer) <- thetaNames[grep("drift", thetaNames)]
  return(flatStdizer)
}




#' exact_getPenaltyValue
#'
#' computes sum(lambda*abs(regularized Values))
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param regIndicators Names of regularized parameters
#' @param lambda Penalty value
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @export
exact_getPenaltyValue <- function(lambda, theta, regIndicators, adaptiveLassoWeights = NULL){

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
#' @param regIndicators Names of regularized parameters
#' @param lambda Penalty value
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param targetVector named vector with values towards which the parameters are regularized
#' @export
exact_getPenaltyValueWithTarget <- function(lambda, theta, regIndicators, targetVector, adaptiveLassoWeights = NULL){

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
#' @export
exact_getSubgradients <- function(theta, jacobian, regIndicators, lambda, lineSearch, adaptiveLassoWeightsMatrix){

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


#' exact_getT0VAR
#'
#' computes the updated T0VAR given the old parameter values in an mxObject and the updates to the parameters in a vector d
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param mxObject mxObject with old parameter values
#' @param d vector with updates to parameter values
#' @param stationarity set to TRUE, if stationaritiy is assumed for the T0VAR
#' @export
exact_getT0VAR <- function(mxObject, d){
  stationarity <- !("T0VAR" %in% names(mxObject$algebras))
  if(stationarity){
    stop("Not yet implemented for stationarity = 'T0VAR'. Changes in derivative necessary")

    T0Varmodel <- mxModel(
      # with free parameters
      mxObject$DRIFT,
      mxObject$DIFFUSIONbase,

      # algebras
      mxObject$T0VAR,
      mxObject$asymDIFFUSIONalg,
      mxObject$DRIFTHATCH,
      mxObject$DIFFUSION,
      mxObject$DIFFUSIONchol,
      mxObject$II)

    # set Parameters
    d_subset <- names(OpenMx::omxGetParameters(T0Varmodel))
    T0Varmodel_new <- OpenMx::omxSetParameters(T0Varmodel, labels = d_subset, values = OpenMx::omxGetParameters(T0Varmodel)+d[d_subset,])

    T0Varmodel_new.fit <- OpenMx::mxRun(T0Varmodel_new, silent = TRUE)

    return(T0Varmodel_new.fit$T0VAR$values)
  }else{
    if(!any(mxObject$T0VARbase$free)){
      return(mxObject$T0VAR$result)
    }
    d_subset <- unique(OpenMx::cvectorize(mxObject$T0VARbase$labels[!is.na(mxObject$T0VARbase$labels)])) # rownames(d)[grepl("T0var", rownames(d))]
    T0VARbaseVal <- mxObject$T0VARbase$values
    T0VARbaseLab <- mxObject$T0VARbase$labels
    for(d_sub in d_subset){
      T0VARbaseVal[T0VARbaseLab == d_sub & !is.na(T0VARbaseLab)] <- T0VARbaseVal[T0VARbaseLab == d_sub & !is.na(T0VARbaseLab)] + d[d_sub,]
    }
    T0VARchol <- OpenMx::vec2diag(exp(OpenMx::diag2vec(T0VARbaseVal))) + T0VARbaseVal - OpenMx::vec2diag(OpenMx::diag2vec(T0VARbaseVal))
    T0VAR <- T0VARchol %*% t(T0VARchol)

    return(T0VAR)

  }

}




#' exact_tryStartingValues
#'
#' tries different starting values to find a good starting point for exact optimization
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param gradientModel Model used for computing gradients (from OpenMx)
#' @param cppmodel cpptsem object to compute gradients
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param currentParameters current parameter values
#' @param sparseParameters parameters of the sparsest model
#' @param regIndicators named vector of regularized parameters
#' @param newLambda new lambda
#' @param adaptiveLassoWeights named vector with adaptive lasso weights
#' @param numStartingValues How many starting values should be tried?
#' @param numOuter number of outer iterations of optim or Rsolnp
#' @export
exact_tryStartingValues <- function(gradientModel, cppmodel,
                                    objective, currentParameters,
                                    sparseParameters, regIndicators,
                                    newLambda,
                                    adaptiveLassoWeights,
                                    numStartingValues, numOuter){
  if(is.null(sparseParameters)){
    stop()
  }
  weights <- seq(-5.2,5.6,length.out = numStartingValues)
  weights <- c(1,1/(1+exp(weights)), 0)
  parameterNames <- names(currentParameters)
  startValuesTable <- matrix(NA, nrow = length(currentParameters), ncol = length(weights))
  regM2LLs <- rep(NA, length(weights))
  rownames(startValuesTable) <- parameterNames

  for(i in seq_len(length(weights))){
    # middle ground between current parameters and sparse parameters
    #weight <- (a^weights[i] - 1)/(a-1)
    weight <- weights[i]
    parameterValues_i <- weight*currentParameters[parameterNames]+(1-weight)*sparseParameters[parameterNames]
    startValuesTable[parameterNames,i] <- parameterValues_i[parameterNames]
    # set parameters
    if(!is.null(cppmodel)){
      cppmodel$setParameterValues(parameterValues_i, names(parameterValues_i))
      if(tolower(objective) == "ml"){
        invisible(capture.output(out1 <- try(cppmodel$computeRAM(), silent = T), type = "message"))
        invisible(capture.output(out2 <- try(cppmodel$fitRAM(), silent = T), type = "message"))
      }else{
        invisible(capture.output(out1 <- try(cppmodel$computeAndFitKalman(), silent = TRUE), type = "message"))
        out2 <- NA
      }
      if(any(class(out1) == "try-error")  |
         any(class(out2) == "try-error")){
        next
      }
      m2LL_i <- cppmodel$m2LL
    }else{
      gradientModel <- OpenMx::omxSetParameters(model = gradientModel, labels = names(parameterValues_i), values = parameterValues_i)
      gradientModel <- try(mxRun(gradientModel, silent = TRUE))
      if(any(class(gradientModel) == "try-error")){
        next
      }
      m2LL_i <- gradientModel$fitfunction$result[[1]]
    }
    # get regularized Likelihood
    regM2LL_i <- m2LL_i + regCtsem::exact_getPenaltyValue(lambda = newLambda,
                                                          theta = parameterValues_i,
                                                          regIndicators = regIndicators,
                                                          adaptiveLassoWeights = adaptiveLassoWeights)
    regM2LLs[i] <- regM2LL_i

    if(!is.na(m2LL_i)){
      optimized <- try(approx_cpptsemOptim(cpptsemmodel = cppmodel,
                                           regM2LLCpptsem = ifelse(tolower(objective) == "ml",
                                                                   regCtsem::approx_RAMRegM2LLCpptsem,
                                                                   regCtsem::approx_KalmanRegM2LLCpptsem),
                                           gradCpptsem = regCtsem::approx_gradCpptsem,
                                           startingValues = parameterValues_i,
                                           adaptiveLassoWeights = adaptiveLassoWeights,
                                           N = 1, lambda = newLambda,
                                           regIndicators = regIndicators,
                                           epsilon = 10^(-8),
                                           maxit = numOuter,
                                           objective,
                                           testGradients = TRUE), silent = TRUE)

      if(any(class(optimized) == "try-error")){
        next
      }


      pars <- optimized$parameters[parameterNames]
      pars[(names(pars) %in% regIndicators) & (abs(pars) < .001)] <- 0
      regM2LLTemp <- try(ifelse(tolower(objective) == "ml",
                                regCtsem::approx_RAMRegM2LLCpptsem(parameters = pars,
                                                                   cpptsemmodel = cppmodel,
                                                                   adaptiveLassoWeights = adaptiveLassoWeights,
                                                                   N = 1, lambda = newLambda,
                                                                   regIndicators = regIndicators,
                                                                   epsilon = 10^(-8), objective = objective,
                                                                   failureReturns = NA
                                ),
                                regCtsem::approx_KalmanRegM2LLCpptsem(parameters = pars,
                                                                      cpptsemmodel = cppmodel,
                                                                      adaptiveLassoWeights = adaptiveLassoWeights,
                                                                      N = 1, lambda = newLambda,
                                                                      regIndicators = regIndicators,
                                                                      epsilon = 10^(-8), objective = objective,
                                                                      failureReturns = NA)))
      if(!any(class(regM2LLTemp) == "try-error") && regM2LLTemp < optimized$regM2LL && regM2LLTemp > -99999){
        optimized$regM2LL <- regM2LLTemp
        optimized$parameters <- pars[parameterNames]
      }

      if(!any(class(optimized) == "try-error") && regM2LLs[i] > optimized$regM2LL){
        parameterValues_i <- optimized$parameters[parameterNames]
        startValuesTable[parameterNames,i] <- parameterValues_i[parameterNames]
        regM2LLs[i] <- optimized$regM2LL
      }
      #optimized <- GIST(gradientModel = gradientModel, cppmodel = cppmodel,
      #                  startingValues = parameterValues_i, objective = objective,
      #                  lambda = newLambda, adaptiveLassoWeights = adaptiveLassoWeights, regularizedParameters = regularizedParameters,
      #                  eta = 1.5, sig = .2, initialStepsize = 1, stepsizeMin = 0, stepsizeMax = 10^30,
      #                  GISTLinesearchCriterion = "monotone", GISTNonMonotoneNBack = 5,
      #                  maxIter_out = numOuter, maxIter_in = 1000,
      #                  break_outer = .00000001, differenceApprox = "central", verbose = 0, silent = TRUE)
      #if(regM2LLs[i] > min(optimized$regM2LL, na.rm = TRUE) & min(optimized$regM2LL, na.rm = TRUE) > 0){
      #  parameterValues_i <- optimized$model$getParameterValues()
      #  startValuesTable[parameterNames,i] <- parameterValues_i[parameterNames]
      #  regM2LLs[i] <- min(optimized$regM2LL, na.rm = TRUE)
      #}
    }

  }
  if(all(is.na(regM2LLs))){stop()}
  startValues <- startValuesTable[parameterNames,which(regM2LLs == min(regM2LLs, na.rm = T))[1]]
  names(startValues) <- parameterNames
  return(startValues)
}


#' exact_initialHessian
#'
#' computes an initial positive-definite Hessian matrix
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param mxObject Fitted object of class MxObject
#' @param approximationType which Hessian should be used? Currently available are "ident" for an identity matrix and "OpenMx" for the Hessian from OpenMx Gradient Descent (recommended)
#' @param estimatedHessian estimated Hessian from OpenMx
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
exact_initialHessian <- function(mxObject,
                                 approximationType,
                                 estimatedHessian){

  numberParameters <- nrow(estimatedHessian)

  # identity matrix
  if(approximationType == "ident"){
    return(diag(numberParameters))
  }

  if(approximationType == "OpenMx"){
    #positiveDefiniteHessianModel <- suppressWarnings(try(mxOption(model = mxObject, key = "Calculate Hessian", value = "No")))
    #positiveDefiniteHessianModel <- suppressWarnings(try(mxRun(positiveDefiniteHessianModel, useOptimizer = TRUE, silent = TRUE)))

    hessian <- mxObject$output$hessian

    if(is.null(hessian)){
      # if the hessian is not returned in the mxObject:
      hessian <- regCtsem::exact_initialHessian(mxObject = mxObject,
                                                approximationType = "ident",
                                                estimatedHessian = estimatedHessian)
    }

    return(hessian)
  }

  # approximation from Nocedal, p. 200
  if(approximationType == "nocedal"){
    initialHessian <- (t(y_1)%*%s_1)/(t(y_1)%*%y_1)%*%diag(numberParameters)
    return(initialHessian)
  }

  stop("Undefined approximation type")
}


#' setStartingValuesFromApprox
#'
#' set the parameter values of mxObject to the values in approx_regModel for the specified lambda
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @param approx_regModel fitted regCtsem with optimization = "approx" and without automatic cross-validation
#' @param mxObject Fitted object of class MxObject
#' @param lambda single lambda
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
setStartingValuesFromApprox <- function(approx_regModel, mxObject, lambda){
  parameterLabels <- names(OpenMx::omxGetParameters(mxObject))
  parameterValues <- approx_regModel$fitAndParameters[parameterLabels, as.character(lambda)]

  newModel <- OpenMx::omxSetParameters(model = mxObject, labels = parameterLabels, values = parameterValues)
  return(newModel)
}


#' tryApproxFirst
#'
#' NOTE: Function located in file exact_optimization.R
#'
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
tryApproxFirst <- function(startingValues, returnAs,
                           approxFirst, numStart, approxMaxIt,
                           lambda, lambdas,
                           gradientModelcpp,
                           mxObject,
                           regIndicators, targetVector = NULL, adaptiveLassoWeights, objective, sparseParameters,
                           extraTries){
  if(approxFirst && is.null(gradientModelcpp)){
    if(any(! (regIndicators %in% mxObject$DRIFT$labels))){
      stop("Requested approxFirst for a model where no C++ model could be built and non-drift parameters are regularized. Cannot build the approximation for non-drift regularization.")
    }

    approxModel <- regCtsem::approx_initializeModel(mxObject = mxObject,
                                                    sampleSize = 1, # scaling with N is handled above
                                                    regOn = "DRIFT",
                                                    regIndicators = regIndicatorsFromNameToMatrix(mxObject = mxObject,
                                                                                                  regOn = regOn,
                                                                                                  regIndicators = regIndicators),
                                                    lambda = lambda,
                                                    adaptiveLassoWeights = adaptiveLassoWeights,
                                                    penalty = "lasso"
    )
    suppressMessages(invisible(capture.output(approxModel <- try(expr = OpenMx::mxTryHardctsem(approxModel, extraTries = extraTries), silent = TRUE))))
    if(!any(class(approxModel) == "try-error")){
      startingValues <- omxGetParameters(approxModel)
    }
    if(returnAs == "vector"){
      return(startingValues)
    }else{
      startingValues <- as.matrix(startingValues)
      return(startingValues)
    }
  }

  # define a target vector if none is provided
  if(is.null(targetVector)){
    targetVector <- rep(0, length(regIndicators))
    names(targetVector) <- regIndicators
  }

  if(approxFirst && !is.null(sparseParameters)){
    temp_startValues <- try(exact_tryStartingValues(gradientModel = gradientModel,
                                                    cppmodel = gradientModelcpp,
                                                    objective = objective,
                                                    currentParameters = startingValues,
                                                    sparseParameters = sparseParameters,
                                                    regIndicators = regIndicators,
                                                    newLambda = lambda,
                                                    adaptiveLassoWeights = adaptiveLassoWeights,
                                                    numStartingValues = numStart,
                                                    numOuter = approxMaxIt), silent = TRUE)
    if(!any(class(temp_startValues) == "try-error")){
      startingValues <- temp_startValues
    }
    if(returnAs == "vector"){
      return(startingValues)
    }else{
      startingValues <- as.matrix(startingValues)
      return(startingValues)
    }
  }
  if(approxFirst && !is.null(gradientModelcpp)){
    optimized <- try(approx_cpptsemOptim(cpptsemmodel = gradientModelcpp,
                                         regM2LLCpptsem = ifelse(tolower(objective) == "ml",
                                                                 regCtsem::approx_RAMRegM2LLCpptsem,
                                                                 regCtsem::approx_KalmanRegM2LLCpptsem),
                                         gradCpptsem = regCtsem::approx_gradCpptsem,
                                         startingValues = startingValues,
                                         adaptiveLassoWeights = adaptiveLassoWeights,
                                         N = 1, lambda = lambda,
                                         regIndicators = regIndicators,
                                         targetVector = targetVector,
                                         epsilon = 10^(-8),
                                         maxit = approxMaxIt,
                                         objective, testGradients = TRUE), silent = TRUE)
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
