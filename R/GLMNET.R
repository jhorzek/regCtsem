#### Main function ####

#' exact_bfgsGLMNET
#'
#' Performs GLMNET (see Friedman, 2010 & Yuan, 2011) with bfgs approximated Hessian
#'
#' NOTE: Function located in file GLMNET.R
#'
#' @param cpptsemObject Fitted object of class cpptsem
#' @param dataset only required if objective = "Kalman". Please provide a data set in wide format compatible to ctsemOMX
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param regIndicators Vector with names of regularized parameters
#' @param lambdas Vector with lambda values that should be tried
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param sparseParameters labeled vector with parameter estimates of the most sparse model.
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = oldParameters + stepSize * Stepdirection)
#' @param lineSearch String indicating which linesearch should be used. Defaults to the one described in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Alternatively (not recommended) Wolfe conditions (lineSearch = "Wolfe") can be used in the outer iteration. Setting to "none" is also not recommended!.
#' @param c1 c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
#' @param c2 c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param initialHessianApproximation Which initial hessian approximation should be used? Possible are: 'ident' for an identity matrix and 'OpenMx' (here the hessian approxmiation from the mxObject is used). If the Hessian from 'OpenMx' is not positive definite, the negative Eigenvalues will be 'flipped' to positive Eigenvalues. This works sometimes, but not always. Alternatively, a matrix can be provided which will be used as initial Hessian
#' @param maxIter_out Maximal number of outer iterations
#' @param maxIter_in Maximal number of inner iterations
#' @param maxIter_line Maximal number of iterations for the lineSearch procedure
#' @param eps_out Stopping criterion for outer iterations
#' @param eps_in Stopping criterion for inner iterations
#' @param eps_WW Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param sampleSize sample size for scaling lambda with N
#' @param approxFirst Should approximate optimization be used first to obtain start values for exact optimization?
#' @param numStart Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param controlApproxOptimizer settings passed to optimx or Rsolnp
#' @param extraTries number of extra tries in mxTryHard for warm start
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @keywords internal
exact_bfgsGLMNET <- function(cpptsemObject,
                             dataset,
                             objective,
                             regIndicators,
                             lambdas,
                             adaptiveLassoWeights,
                             # additional settings
                             sparseParameters,# = NULL,
                             stepSize,# = 1,
                             lineSearch,# = "none",
                             c1,# = .0001,
                             c2,# = .9,
                             sig,# = 10^(-5),
                             gam,# = 0,
                             initialHessianApproximation,# = "optim",
                             maxIter_out,# = 100,
                             maxIter_in,# = 1000,
                             maxIter_line,# = 500,
                             eps_out,# = .0000000001,
                             eps_in,# = .0000000001,
                             eps_WW,# = .0001,
                             scaleLambdaWithN,# = TRUE,
                             sampleSize,
                             approxFirst,# = T,
                             numStart,# = 0,
                             controlApproxOptimizer,
                             extraTries,# = 3,
                             verbose # = 0
){

  # Setup
  # save current parameters
  parameters <- cpptsemObject$getParameterValues()

  ## compute Hessian
  if(!is.matrix(initialHessianApproximation)){
    initialHessian <- try(stats::optimHess(par = parameters, fn = fitCpptsem,
                                    cpptsemObject = cpptsemObject,
                                    objective = objective,
                                    failureReturns = NA
    ), silent = TRUE)
    if(any(class(Hessian) == "try-error")){
      message("Could not compute Hessian. Using Identity")
      initialHessian <- diag(.1, length(parameters))
    }
  }else{
    initialHessian <- initialHessianApproximation
  }

  # get parameter values
  newParameters <- as.matrix(cpptsemObject$getParameterValues())
  parameterLabels <- names(cpptsemObject$getParameterValues())
  initialParameters <- newParameters

  # get initial gradients
  newGradient <- exact_getCppGradients(cpptsemObject, objective = objective)
  initialGradients <- newGradient

  # get initial Hessian
  HessianNew <- initialHessian
  eigenDecomp <- eigen(HessianNew)
  if(any(eigenDecomp$values < 0)){
    message("Initial Hessian is not positive definite. Flipping Eigen values to obtain a positive definite initial Hessian.")
    D <- abs(diag(eigenDecomp$values))
    L <- eigenDecomp$vectors
    HessianNew <- L%*%D%*%solve(L)
  }

  # define return values
  thetas <- matrix(NA,
                   nrow = length(newParameters),
                   ncol = length(lambdas),
                   dimnames = list(parameterLabels,
                                   lambdas))
  m2LL <- rep(NA, length(lambdas))
  regM2LL <- rep(NA, length(lambdas))

  Hessians <- array(NA, dim = c(length(newParameters),length(newParameters),length(lambdas)))

  # Progress bar
  pbar <- utils::txtProgressBar(min = 0, max = length(lambdas), initial = 0, style = 3)

  # iterate over lambda values
  numLambdas <- length(lambdas)
  iteration <- 1
  retryOnce <- TRUE # will retry optimizing once with new starting values
  while(iteration <= numLambdas){
    lambda <- lambdas[iteration]
    # update progress bar
    utils::setTxtProgressBar(pbar, iteration)

    lambda = ifelse(scaleLambdaWithN, lambda*sampleSize, lambda) # set lambda*sampleSize

    # should the results first be approximated?
    startingValues <- as.vector(newParameters)
    names(startingValues) <- rownames(newParameters)
    targetVector <- rep(0, length(regIndicators))
    names(targetVector) <- regIndicators

    if(approxFirst){
      newParameters <- exact_tryStartingValues(startingValues = startingValues,
                                               returnAs = "matrix",
                                               approxFirst = approxFirst,
                                               numStart = numStart,
                                               controlApproxOptimizer = controlApproxOptimizer,
                                               lambda = lambda,
                                               cpptsemObject = cpptsemObject,
                                               regIndicators = regIndicators,
                                               targetVector = targetVector,
                                               adaptiveLassoWeights = adaptiveLassoWeights,
                                               objective = objective,
                                               sparseParameters = sparseParameters)
    }

    # outer loop: optimize parameters
    resGLMNET <- try(exact_outerGLMNET(cpptsemObject = cpptsemObject,
                                                 objective = objective,
                                                 adaptiveLassoWeights = adaptiveLassoWeights,
                                                 sampleSize = sampleSize,
                                                 parameterLabels = parameterLabels,
                                                 initialParameters = newParameters,
                                                 initialGradients = newGradient,
                                                 initialHessian = HessianNew,
                                                 lambda = lambda,
                                                 regIndicators = regIndicators,
                                                 stepSize = stepSize,
                                                 lineSearch = lineSearch,
                                                 c1 = c1,
                                                 c2 = c2,
                                                 sig = sig,
                                                 gam = gam,
                                                 maxIter_out = maxIter_out,
                                                 maxIter_in = maxIter_in,
                                                 maxIter_line = maxIter_line,
                                                 eps_out = eps_out,
                                                 eps_in = eps_in,
                                                 eps_WW = eps_WW,
                                                 scaleLambdaWithN = scaleLambdaWithN,
                                                 verbose = verbose))


    if(any(class(resGLMNET) == "try-error")){

      if(retryOnce){
        message(paste("Model for lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda), " did not converge. Retrying with different starting values.", sep = ""))

      }else{
        message(paste("Model for lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda), " did not converge.", sep = ""))
      }
      # reset theta for next iteration
      newParameters <- initialParameters
      # reset Hessian approximation
      HessianNew <- initialHessian
      # reset gradients
      newGradient <- initialGradients
      if(retryOnce){
        retryOnce <- FALSE
        next
      }
      retryOnce <- TRUE
    }else{
      newParameters <- resGLMNET$newParameters
      # run final model
      cpptsemObject <- resGLMNET$cpptsemObject
      if(tolower(objective) == "ml"){
        out1 <- try(cpptsemObject$computeRAM(), silent = TRUE)
        out2 <- try(cpptsemObject$fitRAM(), silent = TRUE)
      }else{
        out1 <- try(cpptsemObject$computeAndFitKalman(0), silent = TRUE)
        out2 <- NA
      }
      if(any(class(out1)== "try-error") | any(class(out2)== "try-error")){
        message(paste("Model for lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda), "did not converge.", sep = ""))
        # reset theta for next iteration
        newParameters <- initialParameters
        # reset Hessian approximation
        HessianNew <- initialHessian
        # reset gradients
        newGradient <- initialGradients

        if(retryOnce){
          retryOnce <- FALSE
          next
        }
        retryOnce <- TRUE

        # save results
        m2LL <- c(m2LL, NA)
        regM2LL <- c(regM2LL, NA)
      }else{
        newParameters <- resGLMNET$newParameters
        # update Hessian approximation
        HessianNew <- resGLMNET$HessianNew
        # update gradients
        newGradient <- resGLMNET$newGradient
        # save fit
        cM2LL <- ifelse(resGLMNET$convergence, cpptsemObject$m2LL, Inf)
        cRegM2LL <- ifelse(resGLMNET$convergence, cpptsemObject$m2LL +  exact_getPenaltyValue(lambda = lambda,
                                                                                                        theta = resGLMNET$newParameters,
                                                                                                        regIndicators = regIndicators,
                                                                                                        adaptiveLassoWeights = adaptiveLassoWeights), Inf)
        # save results
        thetas[,iteration] <- resGLMNET$newParameters[rownames(thetas),]
        m2LL[iteration] <- cM2LL
        regM2LL[iteration] <- cRegM2LL
        Hessians[,,iteration] <- HessianNew

      }
      iteration <- iteration + 1
    }

  }

  return(list("lambdas" = lambdas, "thetas" = thetas, "m2LL" = m2LL,"regM2LL" = regM2LL, "Hessians" = Hessians))

}


#' exact_outerGLMNET
#'
#' Performs the outer iterations of GLMNET
#'
#' NOTE: Function located in file GLMNET.R
#'
#' @param cpptsemObject Object of type cpptsem
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param sampleSize sample size
#' @param parameterLabels Vector with names of theta-parameters
#' @param initialParameters initial parameter estimates
#' @param initialGradients initial gradients of the likelihood function
#' @param initialHessian initial Hessian of the likelihood function
#' @param lambda Penalty value
#' @param regIndicators Names of regularized parameters
#' @param stepSize initial Stepsize of the outer iteration (theta_{k+1} = oldParameters + stepSize \* Stepdirection) in case of lineSearch
#' @param lineSearch String indicating which linesearch should be used. Defaults to the one described in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Alternatively (not recommended) Wolfe conditions (lineSearch = "Wolfe") can be used in the outer iteration. Setting to "none" is also not recommended!.
#' @param c1 c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
#' @param c2 c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIter_out Maximal number of outer iterations
#' @param maxIter_in Maximal number of inner iterations
#' @param maxIter_line Maximal number of iterations for the lineSearch procedure
#' @param eps_out Stopping criterion for outer iterations
#' @param eps_in Stopping criterion for inner iterations
#' @param eps_WW Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @import OpenMx
#' @keywords internal
exact_outerGLMNET <- function(cpptsemObject,
                              objective,
                              adaptiveLassoWeights,
                              sampleSize,
                              parameterLabels,
                              initialParameters,
                              initialGradients,
                              initialHessian,
                              lambda,
                              regIndicators,
                              stepSize,# = 1,
                              lineSearch,# = "none",
                              c1,# = .0001,
                              c2,# = .9,
                              sig,# = 10^(-5),
                              gam,# = 0,
                              maxIter_out,
                              maxIter_in,
                              maxIter_line,# = 500,
                              eps_out,
                              eps_in,
                              eps_WW,
                              scaleLambdaWithN,
                              verbose){
  # outer loop: optimize parameters
  iter_out <- 0

  converged <- TRUE

  # Convergence Plot
  convergencePlotValues <- NULL
  if(verbose == 1){
    convergencePlotValues <- matrix(NA, nrow = 1, ncol = maxIter_out, dimnames = list("f(theta)", 1:maxIter_out))
  }
  if(verbose == 2){
    convergencePlotValues <- matrix(NA, nrow = length(initialParameters), ncol = maxIter_out, dimnames = list(parameterLabels, 1:maxIter_out))
  }

  newParameters <- initialParameters
  newGradients <- initialGradients
  newHessian <- initialHessian
  newM2LL <- cpptsemObject$m2LL
  newRegM2LL <- newM2LL + exact_getPenaltyValue(lambda = lambda,
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
      convergenceCriterion <- try(max(diag(diag(oldHessian))%*%direction^2) < eps_out, silent = TRUE)
      if(any(class(convergenceCriterion) == "try-error") || is.na(convergenceCriterion)){
        converged <- FALSE
        # save last working parameters
        newParameters <- oldParameters
        newGradients <- oldGradients
        newHessian <- oldHessian
        message(paste("The model did NOT CONVERGE for lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda),". Returning the last working values. These values are NOT acceptable; this can happen if the function is difficult to optimize. try the approximate procedure."))
        break
      }
      if(convergenceCriterion){
        break
      }
    }


    ## inner loop: optimize directions
    direction <- exact_innerGLMNET(
      adaptiveLassoWeights = adaptiveLassoWeights,
      parameterLabels = parameterLabels,
      regIndicators = regIndicators,
      lambda = lambda,
      newParameters = oldParameters,
      newGradient = oldGradients,
      HessianNew = oldHessian,
      maxIter_in = maxIter_in,
      eps_in = eps_in
    )


    # perform Line Search
    if(tolower(lineSearch) == "glmnet") {
      if(stepSize > 1 | stepSize < 0){
        message("Stepsize not allowed. Setting to .9.")
        stepSize = .9
      }
      stepSize_k <- exact_GLMNETLineSearch(cpptsemObject = cpptsemObject,
                                                     objective = objective,
                                                     adaptiveLassoWeights = adaptiveLassoWeights,
                                                     parameterLabels = parameterLabels,
                                                     regIndicators = regIndicators,
                                                     lambda = lambda,
                                                     newParameters = oldParameters,
                                                     m2LLNew = oldM2LL,
                                                     newGradient = oldGradients,
                                                     oldHessian = oldHessian,
                                                     direction = direction,
                                                     stepSize= stepSize,
                                                     sig = sig,
                                                     gam = gam,
                                                     maxIter_line = maxIter_line)
    }else{
      stepSize_k <- stepSize
    }

    newParameters <- oldParameters+stepSize_k*direction

    # update model: set parameter values and compute
    cpptsemObject$setParameterValues(newParameters, parameterLabels)
    if(tolower(objective) == "ml"){
      cpptsemObject$computeRAM()
      cpptsemObject$fitRAM()

      # get fit
      newM2LL <- cpptsemObject$m2LL
      newRegM2LL <- newM2LL + exact_getPenaltyValue(lambda = lambda,
                                                    theta = newParameters,
                                                    regIndicators = regIndicators,
                                                    adaptiveLassoWeights = adaptiveLassoWeights)

      # extract gradients:
      invisible(utils::capture.output(newGradients <- try(exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective))[parameterLabels]))
      newGradients <- matrix(newGradients, nrow = length(newGradients), ncol = 1)
      rownames(newGradients) <- parameterLabels
    }else{
      cpptsemObject$computeAndFitKalman(0)

      # get fit
      newM2LL <- cpptsemObject$m2LL
      newRegM2LL <- newM2LL + exact_getPenaltyValue(lambda = lambda,
                                                    theta = newParameters,
                                                    regIndicators = regIndicators,
                                                    adaptiveLassoWeights = adaptiveLassoWeights)

      # extract gradients:
      invisible(utils::capture.output(newGradients <- try(exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective))[parameterLabels]))
      newGradients <- matrix(newGradients, nrow = length(newGradients), ncol = 1)
      rownames(newGradients) <- parameterLabels
    }

    if(newRegM2LL - oldRegM2LL > 10){
      # the value of 10 is rather arbitrary and is only used to prevent
      # warnings near convergence, when there might be negligible steps in
      # the wrong direction
      message("Step in wrong direction!")
    }

    # Approximate Hessian using bfgs
    newHessian <- exact_getBFGS(oldParameters = oldParameters, oldGradients = oldGradients, oldHessian = oldHessian, newParameters = newParameters, newGradient = newGradients)

    if(verbose == 1){
      convergencePlotValues[,iter_out] <- newRegM2LL
      plot(x=1:maxIter_out, y = convergencePlotValues, xlab = "iteration", ylab = "f(theta)", type = "l", main = "Convergence Plot")

      cat(paste0("\r",
                 "## [", sprintf("%*direction", 3, iter_out),
                 "] m2LL: ", sprintf('%.3f',newM2LL),
                 " | regM2LL:  ", sprintf('%.3f',newRegM2LL),
                 " | zeroed: ", sprintf("%*direction", 3, sum(newParameters[regIndicators,] == 0)),
                 " ##"
      )
      )
      utils::flush.console()
    }

    if(verbose == 2){
      convergencePlotValues[,iter_out] <- newParameters
      graphics::matplot(x=1:maxIter_out, y = t(convergencePlotValues), xlab = "iteration", ylab = "value", type = "l", main = "Convergence Plot")
    }

  }
  # warnings
  if(iter_out == maxIter_out){
    message(paste("For lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda), "the maximum number of iterations was reached. Try with a higher maxIter_out or with smaller lambda-steps."))
  }
  return(list("cpptsemObject" = cpptsemObject, "newParameters" = newParameters, "newGradient" = newGradients, "HessianNew" = newHessian, "convergence" = converged))
}


#' exact_innerGLMNET
#'
#' performs the inner optimization routine of GLMNET. exact_innerGLMNET returns a direction vector for the next step in the outer optimization routine.
#'
#' NOTE: Function located in file GLMNET.R
#'
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param parameterLabels Vector with names of theta-parameters
#' @param regIndicators Names of regularized parameters
#' @param lambda Penalty value
#' @param newParameters Theta at iteration k+1
#' @param newGradient Gradients of the likelihood function at iteration k+1
#' @param HessianNew Hessian of the likelihood function at iteration k+1
#' @param maxIter_in Maximal number of iterations of the inner optimization algorithm
#' @param eps_in Stopping criterion for the inner iterations
#' @keywords internal
exact_innerGLMNET <- function(adaptiveLassoWeights, parameterLabels, regIndicators, lambda, newParameters, newGradient, HessianNew, maxIter_in, eps_in){
  ## inner loop: optimize directions

  direction <- matrix(0,nrow = length(parameterLabels), ncol = 1) # initialize step vector (direction)
  rownames(direction) = parameterLabels

  iter_in <- 0

  while(iter_in < maxIter_in){
    iter_in <- iter_in+1

    # initialize direction z
    z <- matrix(0,nrow = length(parameterLabels), ncol = 1)

    # random order of updates
    d_order <- sample(1:length(direction),length(direction))

    # iterate over parameters
    for(d_i in d_order){

      # compute derivative elements:
      dp_k <- newGradient[d_i]+(HessianNew%*%direction)[d_i]
      d2p_k <- HessianNew[d_i,d_i]

      # if the parameter is regularized:
      if(names(direction[d_i,]) %in% regIndicators){
        theta_i <- newParameters[d_i]
        adaptiveWeight <- ifelse(is.null(adaptiveLassoWeights), 1, adaptiveLassoWeights[d_i])

        # adjust direction for regularized parameters
        if((dp_k-lambda*adaptiveWeight)>=(d2p_k*(theta_i+direction[d_i]))){
          # condition 1
          z_j <- -(dp_k-lambda*adaptiveWeight)/(d2p_k)
          z[d_i] <- z_j
          direction[d_i] <- direction[d_i] + z_j
        }else if((dp_k+lambda*adaptiveWeight)<=(d2p_k*(theta_i+direction[d_i]))){
          # condition 2
          z_j <- -(dp_k+lambda*adaptiveWeight)/(d2p_k)
          z[d_i] <- z_j
          direction[d_i] <- direction[d_i] + z_j
        }else{
          # condition 3
          z_j <- -(theta_i+direction[d_i])
          z[d_i] <- z_j
          direction[d_i] <- direction[d_i]+z_j
        }
      }else{
        # if not regularized: coordinate descent with newton direction
        z_j <- -dp_k/d2p_k
        z[d_i] <- z_j
        direction[d_i] <- direction[d_i]+z_j
      }
    }

    # check inner stopping criterion:
    HessTimesD <- diag(diag(HessianNew))%*%z^2

    breakInnner <- max(HessTimesD)<eps_in
    if(is.na(breakInnner)){
      stop("Error in inner exact_innerGLMNET: The direction z appears to go to infinity.")
    }
    if(breakInnner){
      break
    }
  }

  # return step direction
  return(direction)
}



#### Line search ####

#' exact_armijoLineSearch
#'
#' performs Armijo line search
#'
#' NOTE: Function located in file GLMNET.R
#'
#' @param gradientModel mxObject for computing the derivarive of the likelihood with respect to the parameter estimates
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param parameterLabels names of the parameter estimates
#' @param regIndicators vector with names of parameters to regularize
#' @param lambda penalty value
#' @param newParameters parameter values of iteration k plus 1
#' @param m2LLNew -2 log likelihood of iteration k plus 1
#' @param newGradient gradients of iteration k plus 1
#' @param direction vector with updates to parameter estimates
#' @param c1 tuning parameter for Armijo condition
#' @param stepSize Initial stepsize of the outer iteration (theta_{k+1} = oldParameters + Stepsize \* Stepdirection)
#' @param differenceApprox which approximation for the gradients should be used? Recommended is central
#' @param maxIter_line maximal number of iterations for line search
#' @keywords internal
exact_armijoLineSearch <- function(gradientModel,
                                   adaptiveLassoWeights,
                                   parameterLabels,
                                   regIndicators,
                                   lambda,
                                   newParameters,
                                   m2LLNew,
                                   newGradient,
                                   direction,
                                   c1,
                                   stepSize,
                                   differenceApprox,
                                   maxIter_line){
  stop("lineSearch = 'armijo' is deprecated. Use lineSearch = 'GLMNET'.")
  if(!is.null(adaptiveLassoWeights)){stop("not implemented for adaptive lasso")}
  # Adapted from Huang lslx

  # get penalized M2LL for step size 0:
  h_0 <- m2LLNew + lambda*sum(abs((newParameters)[regIndicators,]))

  # get (sub-)gradients for step size 0:
  g_0 <- exact_getSubgradients(theta = newParameters, jacobian = newGradient, regIndicators = regIndicators, lambda = lambda)

  # Inexact Line Search
  i <- 0

  while(i<maxIter_line){
    i <- i+1 # update iterations
    if(stepSize <.00001){
      return(stepSize)
    }
    # new theta
    parametersNew_td <- newParameters+stepSize*direction
    # get L(x+stepSize*direction) and L'(x+stepSize*direction)
    gradientModel <- OpenMx::omxSetParameters(model = gradientModel, labels = parameterLabels, values = parametersNew_td)
    gradientModel <- OpenMx::mxRun(gradientModel, silent = TRUE)
    m2LL_kp1_td <- gradientModel$fitfunction$result[[1]]
    g_kp1_td <- gradientModel$compute$steps[[1]]$output[["gradient"]]
    g_kp1_td <- g_kp1_td[,differenceApprox] # use specified gradient approximation
    g_kp1_td <- matrix(g_kp1_td, nrow = length(g_kp1_td), ncol = 1)
    rownames(g_kp1_td) <- parameterLabels

    # compute h(stepSize) = L(x+td) + p(x+td) - L(x) - p(x), where p(x) is the penalty function
    h_t <- m2LL_kp1_td + lambda*sum(abs((parametersNew_td)[regIndicators,]))
    # compute h'(stepSize)
    g_t <- exact_getSubgradients(theta = parametersNew_td, jacobian = g_kp1_td, regIndicators = regIndicators, lambda = lambda)

    # Check Armijo
    if(h_t-h_0 <= c1*stepSize*(t(newGradient)%*%direction+lambda*sum(abs((parametersNew_td)[regIndicators,]))-lambda*sum(abs((newParameters)[regIndicators,])))){
      return(stepSize)
    }else{
      stepSize <- stepSize^2
    }
  }
}


#' exact_GLMNETLineSearch
#'
#' performs the line search procedure described by Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421 Equation 20.
#'
#' NOTE: Function located in file GLMNET.R
#'
#' @param cpptsemObject cpptsem object which specifies how the gradients of the likelihood-function are computed
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param parameterLabels names of the parameter estimates
#' @param regIndicators vector with names of parameters to regularize
#' @param lambda penalty value
#' @param newParameters parameter values of iteration k plus 1
#' @param m2LLNew -2 log likelihood of iteration k plus 1
#' @param newGradient gradients of iteration k plus 1
#' @param oldHessian Hessian approximation
#' @param direction vector with updates to parameter estimates
#' @param stepSize Initial stepsize of the outer iteration (theta_{k+1} = oldParameters + Stepsize \* Stepdirection)
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIter_line maximal number of iterations for line search
#' @keywords internal
exact_GLMNETLineSearch <- function(cpptsemObject,
                                   objective,
                                   adaptiveLassoWeights,
                                   parameterLabels,
                                   regIndicators,
                                   lambda,
                                   newParameters,
                                   m2LLNew,
                                   newGradient,
                                   oldHessian,
                                   direction,
                                   stepSize,
                                   sig,
                                   gam,
                                   maxIter_line){

  # deep clone of cpptsemObject
  cpptsemObject_i <- rlang::duplicate(cpptsemObject, shallow =FALSE)

  if(!is.null(adaptiveLassoWeights)){
    adaptiveLassoWeightsMatrix <- diag(adaptiveLassoWeights)
  }else{
    adaptiveLassoWeightsMatrix <- diag(length(parameterLabels))
  }
  rownames(adaptiveLassoWeightsMatrix) <- parameterLabels
  colnames(adaptiveLassoWeightsMatrix) <- parameterLabels

  # get penalized M2LL for step size 0:

  f_0 <- m2LLNew + lambda*sum(abs((adaptiveLassoWeightsMatrix%*%newParameters)[regIndicators,]))
  pen_0 <- lambda*sum(abs((adaptiveLassoWeightsMatrix%*%newParameters)[regIndicators,]))

  i <- 0

  if(stepSize >= 1){
    stepSizeInit <- .9
  }else{
    stepSizeInit <- stepSize
  }
  while(TRUE){
    stepSize <- stepSizeInit^i
    parametersNew_td <- newParameters+stepSize*direction

    # compute new fitfunction value
    cpptsemObject_i$setParameterValues(parametersNew_td, parameterLabels)
    if(tolower(objective) == "ml"){
      invisible(utils::capture.output(out1 <- try(cpptsemObject_i$computeRAM(), silent = TRUE), type = "message"))
      invisible(utils::capture.output(out2 <- try(cpptsemObject_i$fitRAM(), silent = TRUE), type = "message"))
    }else{
      invisible(utils::capture.output(out1 <- try(cpptsemObject_i$computeAndFitKalman(0), silent = TRUE), type = "message"))
      out2 <- NA
    }
    if(any(class(out1)== "try-error") | any(class(out2)== "try-error") | !is.finite(cpptsemObject_i$m2LL)){
      # if the starting values are not feasible
      i <- i+1
      next
    }

    if(any(class(cpptsemObject_i)=="try-error")){
      # if the starting values are far off, the gradients can be very large and testing for the initial
      # step size might result in infeasible parameter values. In this case: try smaller step size
      i <- i+1
      next
    }

    m2LL_kp1_td <- cpptsemObject_i$m2LL

    # compute h(stepSize) = L(x+td) + p(x+td) - L(x) - p(x), where p(x) is the penalty function
    f_new <- m2LL_kp1_td +  lambda*sum(abs((adaptiveLassoWeightsMatrix%*%parametersNew_td)[regIndicators,]))
    p_new <- lambda*sum(abs((adaptiveLassoWeightsMatrix%*%parametersNew_td)[regIndicators,]))

    # test line search criterion
    lineCriterion <- f_new - f_0 <= sig*stepSize*(t(newGradient)%*%direction + gam*t(direction)%*%oldHessian%*%direction + p_new - pen_0)
    if(lineCriterion){
      break
    }
    i <- i+1
    if(i >= maxIter_line){
      #warning("Line search found no stepSize within the maximal number of line search iterations.")
      break
    }
  }
  return(stepSize)
}





#' exact_getBFGS
#'
#' computes the BFGS Hessian approximation
#'
#' NOTE: Function located in file GLMNET.R
#'
#' @param oldParameters Theta at iteration k
#' @param oldGradients Gradients of the likelihood function at iteration k
#' @param oldHessian Hessian of the likelihood function at iteration k
#' @param newParameters Theta at iteration k+1
#' @param newGradient Gradients of the likelihood function at iteration k+1
#' @param cautious boolean: should the update be skipped if it would result in a non positive definite Hessian?
#' @param hessianEps controls when the update of the Hessian approximation is skipped
#' @keywords internal
exact_getBFGS <- function(oldParameters, oldGradients, oldHessian, newParameters, newGradient, cautious = TRUE, hessianEps = .001){

  y <- newGradient-oldGradients
  d <- newParameters-oldParameters

  # test if positive definiteness is ensured
  skipUpdate <- try((t(y)%*%d < hessianEps) && cautious, silent = TRUE)
  if(any(class(skipUpdate) == "try-error") || skipUpdate || is.na(skipUpdate)){
    # Hessian might become non-positive definite. Return without update
    return(oldHessian)
  }
  if(t(y)%*%d < 0){
    message("Hessian update possibly non-positive definite.")
  }

  HessianNew <- oldHessian - (oldHessian%*%d%*%t(d)%*%oldHessian)/as.numeric(t(d)%*%oldHessian%*%d) + (y%*%t(y))/as.numeric(t(y)%*%d)

  if(anyNA(HessianNew)){
    message("Invalid Hessian. Returning previous Hessian")
    return(oldHessian)
  }
  HessianEigen <- eigen(HessianNew)
  iscomplex <- any(!Im(HessianEigen$values) == 0)
  if(iscomplex){
    eigenVectors <- Re(HessianEigen$vectors)
    eigenValues <- Re(HessianEigen$values)
    HessianNew <- eigenVectors%*%diag(eigenValues)%*%t(eigenVectors)
  }
  if(any(HessianEigen$values < 0)){
    while(any(eigen(HessianNew)$values < 0)){
      HessianNew <- HessianNew + .01*diag(nrow(HessianNew))
    }
  }

  return(HessianNew)

}






