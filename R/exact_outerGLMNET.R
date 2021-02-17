#' exact_outerGLMNET
#'
#' Performs the outer iterations of GLMNET
#'
#' @param mxObject Object of type MxModel
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param sampleSize sample size
#' @param gradientModel Object of Type MxModel which specifies how the gradients of the likelihood-function are computed (the jacobian)
#' @param gradientModelcpp cpptsem object which specifies how the gradients of the likelihood-function are computed (the jacobian)
#' @param parameterNames Vector with names of theta-parameters
#' @param initialParameters initial parameter estimates
#' @param initialGradients initial gradients of the likelihood function
#' @param initialHessian initial Hessian of the likelihood function
#' @param lambda Penalty value
#' @param regIndicators Names of regularized parameters
#' @param stepSize initial Stepsize of the outer iteration (theta_{k+1} = theta_k + stepSize \* Stepdirection) in case of lineSearch
#' @param lineSearch String indicating if Wolfe conditions (lineSearch = "Wolfe") should be used in the outer iteration. Set lineSearch = "dynamic" to only use line search if optimization without line search fails.
#' @param c1 c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
#' @param c2 c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
#' @param sig GLMNET & GIST: GLMNET: only relevant when lineSearch = 'GLMNET' | GIST: sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
#' @param gam GLMNET when lineSearch = 'GLMNET'. Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param differenceApprox Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#' @param maxIter_out Maximal number of outer iterations
#' @param maxIter_in Maximal number of inner iterations
#' @param maxIter_line Maximal number of iterations for the lineSearch procedure
#' @param eps_out Stopping criterion for outer iterations
#' @param eps_in Stopping criterion for inner iterations
#' @param eps_WW Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned
#' @param eps_numericDerivative controls the precision of the central gradient approximation. The default (1.1 * 10^(-16))^(1/3) is derived in Nocedal, J., & Wright, S. J. (2006). Numerical optimization (2nd ed), p. 197
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @import OpenMx
#' @export
exact_outerGLMNET <- function(mxObject, objective, adaptiveLassoWeights, sampleSize, gradientModel, gradientModelcpp = NULL, parameterNames, initialParameters, initialGradients, initialHessian,
                              lambda, regIndicators,
                              stepSize = 1, lineSearch = "none", c1 = .0001, c2 = .9, sig = .2, gam = 0,
                              differenceApprox = "central",
                              maxIter_out, maxIter_in, maxIter_line = 100, eps_out, eps_in, eps_WW, eps_numericDerivative = (1.1 * 10^(-16))^(1/3),
                              scaleLambdaWithN,
                              verbose = 0){
  # outer loop: optimize parameters
  iter_out <- 0

  converged <- TRUE

  # Convergence Plot
  convergencePlotValues <- NULL
  if(verbose == 1){
    convergencePlotValues <- matrix(NA, nrow = 1, ncol = maxIter_out, dimnames = list("f(theta)", 1:maxIter_out))
  }
  if(verbose == 2){
    convergencePlotValues <- matrix(NA, nrow = length(initialParameters), ncol = maxIter_out, dimnames = list(parameterNames, 1:maxIter_out))
  }

  newParameters <- initialParameters
  newGradients <- initialGradients
  newHessian <- initialHessian
  if(!is.null(gradientModelcpp)){
    newM2LL <- gradientModelcpp$m2LL
  }else{
    newM2LL <- gradientModel$fitfunction$result[[1]]
  }
  newRegM2LL <- newM2LL + exact_getRegValue(lambda = lambda,
                                            theta = newParameters,
                                            regIndicators = regIndicators,
                                            adaptiveLassoWeights = adaptiveLassoWeights)

  while(iter_out < maxIter_out){
    iter_out <- iter_out +1

    oldParameters <- newParameters
    oldGradients <- newGradients
    oldHessian <- newHessian

    oldM2LL <- newM2LL
    oldRegM2LL <- newRegM2LL

    # outer breaking condition
    if(iter_out > 1){
      # requires a direction vector d; therefore, we cannot evaluate the
      # convergence in the first iteration
      convergenceCriterion <- max(diag(diag(oldHessian))%*%d^2) < eps_out
      if(is.na(convergenceCriterion)){
        converged <- FALSE
        # save last working parameters
        newParameters <- oldParameters
        newGradients <- oldGradients
        newHessian <- oldHessian
        warning(paste("The model did NOT CONVERGE for lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda),". Returning the last working values. These values are NOT acceptable; this can happen if the function is difficult to optimize. try the approximate procedure."))
        break
      }
      if(convergenceCriterion){
        break
      }
    }


    ## inner loop: optimize directions
    d <- regCtsem::exact_innerGLMNET(
      adaptiveLassoWeights = adaptiveLassoWeights,
      thetaNames = parameterNames,
      regIndicators = regIndicators,
      lambda = lambda,
      theta_kp1 = oldParameters,
      g_kp1 = oldGradients,
      H_kp1 = oldHessian,
      maxIter_in = maxIter_in,
      eps_in = eps_in)


    # perform Line Search
    if(tolower(lineSearch) == "wolfe"){

      stepSize_k <- regCtsem::exact_weakWolfeLineSearch(gradientModel = gradientModel, objective = objective,
                                                         gradientModelcpp = gradientModelcpp,
                                                         adaptiveLassoWeights = adaptiveLassoWeights, thetaNames = parameterNames,
                                                         regIndicators = regIndicators, lambda = lambda,
                                                         theta_kp1 = oldParameters, m2LL_kp1 = oldM2LL, g_kp1 = oldGradients,
                                                         d = d, differenceApprox = differenceApprox, eps_numericDerivative = eps_numericDerivative,
                                                         stepSize= stepSize, c1 = c1, c2 = c2, maxIter_line = maxIter_line,
                                                         eps_WW = eps_WW, verbose = verbose)
    }else if(tolower(lineSearch) == "glmnet") {
      if(stepSize > 1 | stepSize < 0){
        warning("Stepsize not allowed. Setting to .9.")
        stepSize = .9
      }
      stepSize_k <- regCtsem::exact_GLMNETLineSearch(gradientModel = gradientModel, objective = objective,
                                                      gradientModelcpp = gradientModelcpp,
                                                      adaptiveLassoWeights = adaptiveLassoWeights, thetaNames = parameterNames,
                                                      regIndicators = regIndicators, lambda = lambda,
                                                      theta_kp1 = oldParameters, m2LL_kp1 = oldM2LL, g_kp1 = oldGradients, H_k = oldHessian,
                                                      d = d, differenceApprox = differenceApprox, eps_numericDerivative = eps_numericDerivative,
                                                      stepSize= stepSize, sig = sig, gam = gam, maxIter_line = maxIter_line)
    }else{
      stepSize_k <- stepSize
    }

    newParameters <- oldParameters+stepSize_k*d

    # update model: set parameter values and compute
    if(!is.null(gradientModelcpp)){
      gradientModelcpp$setParameterValues(newParameters, parameterNames)
      if(tolower(objective) == "ml"){
        gradientModelcpp$computeRAM()
        gradientModelcpp$fitRAM()

        # get fit
        newM2LL <- gradientModelcpp$m2LL
        newRegM2LL <- newM2LL + exact_getRegValue(lambda = lambda,
                                                  theta = newParameters,
                                                  regIndicators = regIndicators,
                                                  adaptiveLassoWeights = adaptiveLassoWeights)

        # extract gradients:
        newGradients <-  gradientModelcpp$approxRAMGradients(eps_numericDerivative)[parameterNames]
        newGradients <- matrix(newGradients, nrow = length(newGradients), ncol = 1)
        rownames(newGradients) <- parameterNames
      }else{
        gradientModelcpp$computeAndFitKalman()

        # get fit
        newM2LL <- gradientModelcpp$m2LL
        newRegM2LL <- newM2LL + exact_getRegValue(lambda = lambda,
                                                  theta = newParameters,
                                                  regIndicators = regIndicators,
                                                  adaptiveLassoWeights = adaptiveLassoWeights)

        # extract gradients:
        newGradients <-  gradientModelcpp$approxKalmanGradients(eps_numericDerivative)[parameterNames]
        newGradients <- matrix(newGradients, nrow = length(newGradients), ncol = 1)
        rownames(newGradients) <- parameterNames
      }

    }else{
      gradientModel <- OpenMx::omxSetParameters(model = gradientModel, labels = parameterNames, values = newParameters)
      gradientModel <- suppressWarnings(try(OpenMx::mxRun(gradientModel, silent = TRUE), silent = TRUE))
      # get fit
      newM2LL <- gradientModel$fitfunction$result[[1]]
      newRegM2LL <- newM2LL + exact_getRegValue(lambda = lambda,
                                                theta = newParameters,
                                                regIndicators = regIndicators,
                                                adaptiveLassoWeights = adaptiveLassoWeights)

      # extract gradients:
      newGradients <- gradientModel$compute$steps[[1]]$output[["gradient"]]
      newGradients <- newGradients[,differenceApprox] # use specified gradient approximation
      newGradients <- matrix(newGradients, nrow = length(newGradients), ncol = 1)
      rownames(newGradients) <- parameterNames
    }

    if(newRegM2LL - oldRegM2LL > 10){
      # the value of 10 is rather arbitrary and is only used to prevent
      # warnings near convergence, when there might be negligible steps in
      # the wrong direction
      warning("Step in wrong direction!")
    }

    # Approximate Hessian using bfgs
    newHessian <- regCtsem::exact_getBFGS(theta_k = oldParameters, g_k = oldGradients, H_k = oldHessian, theta_kp1 = newParameters, g_kp1 = newGradients)


    if(verbose == 1){
      convergencePlotValues[,iter_out] <- newRegM2LL
      plot(x=1:maxIter_out, y = convergencePlotValues, xlab = "iteration", ylab = "f(theta)", type = "l", main = "Convergence Plot")

      cat(paste0("\r",
                 "## [", sprintf("%*d", 3, iter_out),
                 "] m2LL: ", sprintf('%.3f',newM2LL),
                 " | regM2LL:  ", sprintf('%.3f',newRegM2LL),
                 " | zeroed: ", sprintf("%*d", 3, sum(newParameters[regIndicators,] == 0)),
                 " ##"
      )
      )
      flush.console()
    }

    if(verbose == 2){
      convergencePlotValues[,iter_out] <- newParameters
      matplot(x=1:maxIter_out, y = t(convergencePlotValues), xlab = "iteration", ylab = "value", type = "l", main = "Convergence Plot")
    }

  }
  # warnings
  if(iter_out == maxIter_out){
    warning(paste("For lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda), "the maximum number of iterations was reached. Try with a higher maxIter_out or with smaller lambda-steps."))
  }
  return(list("gradientModel" = gradientModel, "gradientModelcpp" = gradientModelcpp, "theta_kp1" = newParameters, "g_kp1" = newGradients, "H_kp1" = newHessian, "convergence" = converged))
}



