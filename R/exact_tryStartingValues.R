#' exact_tryStartingValues
#'
#' tries different starting values to find a good starting point for exact optimization
#' @param gradientModel Model used for computing gradients (from OpenMx)
#' @param cppmodel cpptsem object to compute gradients
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param currentParameters current parameter values
#' @param sparseParameters parameters of the sparsest model
#' @param regIndicators named vector of regularized parameters
#' @param newLambda new regValue
#' @param adaptiveLassoWeights named vector with adaptive lasso weights
#' @param differenceApprox difference approximation for gradients (detault should be central)
#' @param numStartingValues How many starting values should be tried?
#' @param optimize 1 = optimization with optim, 2 = optimization with Rsolnp
#' @param numOuter number of outer iterations of optim or Rsolnp
#' @export
exact_tryStartingValues <- function(gradientModel, cppmodel,
                                    objective, currentParameters,
                                    sparseParameters, regIndicators,
                                    newLambda,
                                    adaptiveLassoWeights, differenceApprox,
                                    numStartingValues, optimize, numOuter){
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
    regM2LL_i <- m2LL_i + regCtsem::exact_getRegValue(lambda = newLambda,
                                                       theta = parameterValues_i,
                                                       regIndicators = regIndicators,
                                                       adaptiveLassoWeights = adaptiveLassoWeights)
    regM2LLs[i] <- regM2LL_i

    if((optimize > 0) & !is.na(m2LL_i)){
      if(optimize == 1){
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
      }else{
        optimized <- try(approx_cpptsemSolnp(cpptsemmodel = cppmodel,
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
      }

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
      if(!any(class(regM2LLTemp) == "try-error") && regM2LLTemp < optimized$regM2LL && regM2LLTemp > 0){
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




