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
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = theta_k + stepSize * Stepdirection)
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
#' @export
exact_bfgsGLMNET <- function(cpptsemObject, dataset, objective, regIndicators, lambdas, adaptiveLassoWeights,
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
    Hessian <- try(optimHess(par = parameters, fn = fitCpptsem,
                             cpptsemObject = cpptsemObject,
                             objective = objective,
                             failureReturns = NA
    ), silent = TRUE)
    if(any(class(Hessian) == "try-error")){
      warning("Could not compute Hessian. Using Identity")
      Hessian <- diag(.1, length(parameters))
    }
  }else{
    Hessian <- initialHessianApproximation
  }

  # get parameter values
  theta_kp1 <- as.matrix(cpptsemObject$getParameterValues())
  thetaNames <- names(cpptsemObject$getParameterValues())
  initialParameters <- theta_kp1

  # get initial gradients
  g_kp1 <- regCtsem::exact_getCppGradients(cpptsemObject, objective = objective)
  initialGradients <- g_kp1

  # get initial Hessian
  H_kp1 <- Hessian
  eigenDecomp <- eigen(H_kp1)
  if(any(eigenDecomp$values < 0)){
    warning("Initial Hessian is not positive definite. Flipping Eigen values to obtain a positive definite initial Hessian.")
    D <- abs(diag(eigenDecomp$values))
    L <- eigenDecomp$vectors
    H_kp1 <- L%*%D%*%solve(L)
  }

  # define return values
  thetas <- matrix(NA,
                   nrow = length(theta_kp1),
                   ncol = length(lambdas),
                   dimnames = list(thetaNames,
                                   lambdas))
  m2LL <- rep(NA, length(lambdas))
  regM2LL <- rep(NA, length(lambdas))


  # Progress bar
  pbar <- txtProgressBar(min = 0, max = length(lambdas), initial = 0, style = 3)

  # iterate over lambda values
  numLambdas <- length(lambdas)
  iteration <- 1
  retryOnce <- TRUE # will retry optimizin once with new starting values
  while(iteration <= numLambdas){
    lambda <- lambdas[iteration]
    # update progress bar
    setTxtProgressBar(pbar, iteration)

    lambda = ifelse(scaleLambdaWithN, lambda*sampleSize, lambda) # set lambda*sampleSize

    # should the results first be approximated?
    startingValues <- as.vector(theta_kp1)
    names(startingValues) <- rownames(theta_kp1)
    targetVector <- rep(0, length(regIndicators))
    names(targetVector) <- regIndicators

    if(approxFirst){
      theta_kp1 <- exact_tryStartingValues(startingValues = startingValues, returnAs = "matrix",
                                  approxFirst = approxFirst, numStart = numStart, controlApproxOptimizer = controlApproxOptimizer,
                                  lambda = lambda,
                                  cpptsemObject = cpptsemObject,
                                  regIndicators = regIndicators, targetVector = targetVector,
                                  adaptiveLassoWeights = adaptiveLassoWeights, objective = objective, sparseParameters = sparseParameters)
    }

    # outer loop: optimize parameters
    resGLMNET <- try(regCtsem::exact_outerGLMNET(cpptsemObject = cpptsemObject, objective = objective,
                                                 adaptiveLassoWeights = adaptiveLassoWeights, sampleSize = sampleSize,
                                                 parameterNames = thetaNames,
                                                 initialParameters = theta_kp1, initialGradients = g_kp1, initialHessian = H_kp1,
                                                 lambda = lambda, regIndicators = regIndicators,
                                                 stepSize = stepSize, lineSearch = lineSearch, c1 = c1, c2 = c2, sig = sig, gam = gam,
                                                 maxIter_out = maxIter_out,
                                                 maxIter_in = maxIter_in, maxIter_line = maxIter_line, eps_out = eps_out,
                                                 eps_in = eps_in, eps_WW = eps_WW, scaleLambdaWithN = scaleLambdaWithN, verbose = verbose))


    if(any(class(resGLMNET) == "try-error")){

      if(retryOnce){
        warning(paste("Model for lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda), " did not converge. Retrying with different starting values.", sep = ""))

      }else{
        warning(paste("Model for lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda), " did not converge.", sep = ""))
      }
      # reset theta for next iteration
      theta_kp1 <- initialParameters
      # reset Hessian approximation
      H_kp1 <- initialHessian
      # reset gradients
      g_kp1 <- initialGradients
      if(retryOnce){
        retryOnce <- FALSE
        next
      }
      retryOnce <- TRUE
    }else{
      theta_kp1 <- resGLMNET$theta_kp1
      # run final model
        cpptsemObject <- resGLMNET$cpptsemObject
        if(tolower(objective) == "ml"){
          out1 <- try(cpptsemObject$computeRAM(), silent = TRUE)
          out2 <- try(cpptsemObject$fitRAM(), silent = TRUE)
        }else{
          out1 <- try(cpptsemObject$computeAndFitKalman(), silent = TRUE)
          out2 <- NA
        }
        if(any(class(out1)== "try-error") | any(class(out2)== "try-error")){
          warning(paste("Model for lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda), "did not converge.", sep = ""))
          # reset theta for next iteration
          theta_kp1 <- initialParameters
          # reset Hessian approximation
          H_kp1 <- initialHessian
          # reset gradients
          g_kp1 <- initialGradients

          if(retryOnce){
            retryOnce <- FALSE
            next
          }
          retryOnce <- TRUE

          # save results
          m2LL <- c(m2LL, NA)
          regM2LL <- c(regM2LL, NA)
        }else{
          theta_kp1 <- resGLMNET$theta_kp1
          # update Hessian approximation
          H_kp1 <- resGLMNET$H_kp1
          # update gradients
          g_kp1 <- resGLMNET$g_kp1
          # save fit
          cM2LL <- ifelse(resGLMNET$convergence, cpptsemObject$m2LL, Inf)
          cRegM2LL <- ifelse(resGLMNET$convergence, cpptsemObject$m2LL +  regCtsem::exact_getPenaltyValue(lambda = lambda,
                                                                                                          theta = resGLMNET$theta_kp1,
                                                                                                          regIndicators = regIndicators,
                                                                                                          adaptiveLassoWeights = adaptiveLassoWeights), Inf)
          # save results
          thetas[,iteration] <- resGLMNET$theta_kp1[rownames(thetas),]
          m2LL[iteration] <- cM2LL
          regM2LL[iteration] <- cRegM2LL
        }
      iteration <- iteration + 1
    }

  }

  return(list("lambdas" = lambdas, "thetas" = thetas, "m2LL" = m2LL,"regM2LL" = regM2LL))

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
#' @param parameterNames Vector with names of theta-parameters
#' @param initialParameters initial parameter estimates
#' @param initialGradients initial gradients of the likelihood function
#' @param initialHessian initial Hessian of the likelihood function
#' @param lambda Penalty value
#' @param regIndicators Names of regularized parameters
#' @param stepSize initial Stepsize of the outer iteration (theta_{k+1} = theta_k + stepSize \* Stepdirection) in case of lineSearch
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
#' @export
exact_outerGLMNET <- function(cpptsemObject,
                              objective,
                              adaptiveLassoWeights,
                              sampleSize,
                              parameterNames,
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
    convergencePlotValues <- matrix(NA, nrow = length(initialParameters), ncol = maxIter_out, dimnames = list(parameterNames, 1:maxIter_out))
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
      eps_in = eps_in
    )


    # perform Line Search
    if(tolower(lineSearch) == "wolfe"){

      stepSize_k <- regCtsem::exact_weakWolfeLineSearch(cpptsemObject = cpptsemObject, objective = objective,
                                                        adaptiveLassoWeights = adaptiveLassoWeights, thetaNames = parameterNames,
                                                        regIndicators = regIndicators, lambda = lambda,
                                                        theta_kp1 = oldParameters, m2LL_kp1 = oldM2LL, g_kp1 = oldGradients,
                                                        d = d,
                                                        stepSize= stepSize, c1 = c1, c2 = c2, maxIter_line = maxIter_line,
                                                        eps_WW = eps_WW, verbose = verbose)
    }else if(tolower(lineSearch) == "glmnet") {
      if(stepSize > 1 | stepSize < 0){
        warning("Stepsize not allowed. Setting to .9.")
        stepSize = .9
      }
      stepSize_k <- regCtsem::exact_GLMNETLineSearch(cpptsemObject = cpptsemObject, objective = objective,
                                                     adaptiveLassoWeights = adaptiveLassoWeights, thetaNames = parameterNames,
                                                     regIndicators = regIndicators, lambda = lambda,
                                                     theta_kp1 = oldParameters, m2LL_kp1 = oldM2LL, g_kp1 = oldGradients, H_k = oldHessian,
                                                     d = d,
                                                     stepSize= stepSize, sig = sig, gam = gam, maxIter_line = maxIter_line)
    }else{
      stepSize_k <- stepSize
    }

    newParameters <- oldParameters+stepSize_k*d

    # update model: set parameter values and compute
    cpptsemObject$setParameterValues(newParameters, parameterNames)
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
      invisible(capture.output(newGradients <- try(exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective))[parameterNames]))
      newGradients <- matrix(newGradients, nrow = length(newGradients), ncol = 1)
      rownames(newGradients) <- parameterNames
    }else{
      cpptsemObject$computeAndFitKalman()

      # get fit
      newM2LL <- cpptsemObject$m2LL
      newRegM2LL <- newM2LL + exact_getPenaltyValue(lambda = lambda,
                                                    theta = newParameters,
                                                    regIndicators = regIndicators,
                                                    adaptiveLassoWeights = adaptiveLassoWeights)

      # extract gradients:
      invisible(capture.output(newGradients <- try(exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective))[parameterNames]))
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
  return(list("cpptsemObject" = cpptsemObject, "theta_kp1" = newParameters, "g_kp1" = newGradients, "H_kp1" = newHessian, "convergence" = converged))
}


#' exact_innerGLMNET
#'
#' performs the inner optimization routine of GLMNET. exact_innerGLMNET returns a direction vector for the next step in the outer optimization routine.
#'
#' NOTE: Function located in file GLMNET.R
#'
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param thetaNames Vector with names of theta-parameters
#' @param regIndicators Names of regularized parameters
#' @param lambda Penalty value
#' @param theta_kp1 Theta at iteration k+1
#' @param g_kp1 Gradients of the likelihood function at iteration k+1
#' @param H_kp1 Hessian of the likelihood function at iteration k+1
#' @param maxIter_in Maximal number of iterations of the inner optimization algorithm
#' @param eps_in Stopping criterion for the inner iterations
#' @export
exact_innerGLMNET <- function(adaptiveLassoWeights, thetaNames, regIndicators, lambda, theta_kp1, g_kp1, H_kp1, maxIter_in, eps_in){
  ## inner loop: optimize directions

  d <- matrix(0,nrow = length(thetaNames), ncol = 1) # initialize step vector (direction d)
  rownames(d) = thetaNames

  iter_in <- 0

  while(iter_in < maxIter_in){
    iter_in <- iter_in+1

    # initialize direction z
    z <- matrix(0,nrow = length(thetaNames), ncol = 1)

    # random order of updates
    d_order <- sample(1:length(d),length(d))

    # iterate over parameters
    for(d_i in d_order){

      # compute derivative elements:
      dp_k <- g_kp1[d_i]+(H_kp1%*%d)[d_i]
      d2p_k <- H_kp1[d_i,d_i]

      # if the parameter is regularized:
      if(names(d[d_i,]) %in% regIndicators){
        theta_i <- theta_kp1[d_i]
        adaptiveWeight <- ifelse(is.null(adaptiveLassoWeights), 1, adaptiveLassoWeights[d_i])

        # adjust d for regularized parameters
        if((dp_k-lambda*adaptiveWeight)>=(d2p_k*(theta_i+d[d_i]))){
          # condition 1
          z_j <- -(dp_k-lambda*adaptiveWeight)/(d2p_k)
          z[d_i] <- z_j
          d[d_i] <- d[d_i] + z_j
        }else if((dp_k+lambda*adaptiveWeight)<=(d2p_k*(theta_i+d[d_i]))){
          # condition 2
          z_j <- -(dp_k+lambda*adaptiveWeight)/(d2p_k)
          z[d_i] <- z_j
          d[d_i] <- d[d_i] + z_j
        }else{
          # condition 3
          z_j <- -(theta_i+d[d_i])
          z[d_i] <- z_j
          d[d_i] <- d[d_i]+z_j
        }
      }else{
        # if not regularized: coordinate descent with newton direction
        z_j <- -dp_k/d2p_k
        z[d_i] <- z_j
        d[d_i] <- d[d_i]+z_j
      }
    }

    # check inner stopping criterion:
    HessTimesD <- diag(diag(H_kp1))%*%z^2

    breakInnner <- max(HessTimesD)<eps_in
    if(is.na(breakInnner)){
      stop("Error in inner exact_innerGLMNET: The direction z appears to go to infinity.")
    }
    if(breakInnner){
      break
    }
  }

  # return step direction
  return(d)
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
#' @param thetaNames names of the parameter estimates
#' @param regIndicators vector with names of parameters to regularize
#' @param lambda penalty value
#' @param theta_kp1 parameter values of iteration k plus 1
#' @param m2LL_kp1 -2 log likelihood of iteration k plus 1
#' @param g_kp1 gradients of iteration k plus 1
#' @param d vector with updates to parameter estimates
#' @param c1 tuning parameter for Armijo condition
#' @param stepSize Initial stepsize of the outer iteration (theta_{k+1} = theta_k + Stepsize \* Stepdirection)
#' @param differenceApprox which approximation for the gradients should be used? Recommended is central
#' @param maxIter_line maximal number of iterations for line search
#'
#' @export
exact_armijoLineSearch <- function(gradientModel,
                                   adaptiveLassoWeights,
                                   thetaNames,
                                   regIndicators,
                                   lambda,
                                   theta_kp1, m2LL_kp1, g_kp1, d, c1, stepSize, differenceApprox, maxIter_line){
  stop("lineSearch = 'armijo' is deprecated. Use lineSearch = 'GLMNET'.")
  if(!is.null(adaptiveLassoWeights)){stop("not implemented for adaptive lasso")}
  # Adapted from Huang lslx

  # get penalized M2LL for step size 0:
  h_0 <- m2LL_kp1 + lambda*sum(abs((theta_kp1)[regIndicators,]))

  # get (sub-)gradients for step size 0:
  g_0 <- regCtsem::exact_getSubgradients(theta = theta_kp1, jacobian = g_kp1, regIndicators = regIndicators, lambda = lambda)

  # Inexact Line Search
  i <- 0

  while(i<maxIter_line){
    i <- i+1 # update iterations
    if(stepSize <.00001){
      return(stepSize)
    }
    # new theta
    theta_kp1_td <- theta_kp1+stepSize*d
    # get L(x+stepSize*d) and L'(x+stepSize*d)
    gradientModel <- OpenMx::omxSetParameters(model = gradientModel, labels = thetaNames, values = theta_kp1_td)
    gradientModel <- OpenMx::mxRun(gradientModel, silent = TRUE)
    m2LL_kp1_td <- gradientModel$fitfunction$result[[1]]
    g_kp1_td <- gradientModel$compute$steps[[1]]$output[["gradient"]]
    g_kp1_td <- g_kp1_td[,differenceApprox] # use specified gradient approximation
    g_kp1_td <- matrix(g_kp1_td, nrow = length(g_kp1_td), ncol = 1)
    rownames(g_kp1_td) <- thetaNames

    # compute h(stepSize) = L(x+td) + p(x+td) - L(x) - p(x), where p(x) is the penalty function
    h_t <- m2LL_kp1_td + lambda*sum(abs((theta_kp1_td)[regIndicators,]))
    # compute h'(stepSize)
    g_t <- regCtsem::exact_getSubgradients(theta = theta_kp1_td, jacobian = g_kp1_td, regIndicators = regIndicators, lambda = lambda)

    # Check Armijo
    if(h_t-h_0 <= c1*stepSize*(t(g_kp1)%*%d+lambda*sum(abs((theta_kp1_td)[regIndicators,]))-lambda*sum(abs((theta_kp1)[regIndicators,])))){
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
#' @param thetaNames names of the parameter estimates
#' @param regIndicators vector with names of parameters to regularize
#' @param lambda penalty value
#' @param theta_kp1 parameter values of iteration k plus 1
#' @param m2LL_kp1 -2 log likelihood of iteration k plus 1
#' @param g_kp1 gradients of iteration k plus 1
#' @param H_k Hessian approximation
#' @param d vector with updates to parameter estimates
#' @param stepSize Initial stepsize of the outer iteration (theta_{k+1} = theta_k + Stepsize \* Stepdirection)
#' @param sig only relevant when lineSearch = 'GLMNET'. Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421.
#' @param gam Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIter_line maximal number of iterations for line search
#' @export
exact_GLMNETLineSearch <- function(cpptsemObject,
                                   objective,
                                   adaptiveLassoWeights,
                                   thetaNames,
                                   regIndicators,
                                   lambda,
                                   theta_kp1,
                                   m2LL_kp1,
                                   g_kp1,
                                   H_k,
                                   d,
                                   stepSize,
                                   sig,
                                   gam,
                                   maxIter_line){

  # deep clone of cpptsemObject
  cpptsemObject_i <- rlang::duplicate(cpptsemObject, shallow =FALSE)

  if(!is.null(adaptiveLassoWeights)){
    adaptiveLassoWeightsMatrix <- diag(adaptiveLassoWeights)
  }else{
    adaptiveLassoWeightsMatrix <- diag(length(thetaNames))
  }
  rownames(adaptiveLassoWeightsMatrix) <- thetaNames
  colnames(adaptiveLassoWeightsMatrix) <- thetaNames

  # get penalized M2LL for step size 0:

  f_0 <- m2LL_kp1 + lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1)[regIndicators,]))
  pen_0 <- lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1)[regIndicators,]))

  i <- 0

  if(stepSize >= 1){
    stepSizeInit <- .9
  }else{
    stepSizeInit <- stepSize
  }
  while(TRUE){
    stepSize <- stepSizeInit^i
    theta_kp1_td <- theta_kp1+stepSize*d

    # compute new fitfunction value
    cpptsemObject_i$setParameterValues(theta_kp1_td, thetaNames)
    if(tolower(objective) == "ml"){
      invisible(capture.output(out1 <- try(cpptsemObject_i$computeRAM(), silent = TRUE), type = "message"))
      invisible(capture.output(out2 <- try(cpptsemObject_i$fitRAM(), silent = TRUE), type = "message"))
    }else{
      invisible(capture.output(out1 <- try(cpptsemObject_i$computeAndFitKalman(), silent = TRUE), type = "message"))
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
    f_new <- m2LL_kp1_td +  lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1_td)[regIndicators,]))
    p_new <- lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1_td)[regIndicators,]))

    # test line search criterion
    lineCriterion <- f_new - f_0 <= sig*stepSize*(t(g_kp1)%*%d + gam*t(d)%*%H_k%*%d + p_new - pen_0)
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



#' exact_weakWolfeLineSearch
#'
#'
#' NOTE: Function located in file GLMNET.R
#'
#' @param cpptsemObject Object of type cpptsem
#' @param objective Kalman or ML
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param thetaNames names of the parameter estimates
#' @param regIndicators vector with names of parameters to regularize
#' @param lambda penalty value
#' @param theta_kp1 parameter values of iteration k plus 1
#' @param m2LL_kp1 -2 log likelihood of iteration k plus 1
#' @param g_kp1 gradients of iteration k plus 1
#' @param d vector with updates to parameter estimates
#' @param c1 tuning parameter for Armijo condition
#' @param c2 c2 constant for lineSearch. This constant controls the Curvature condition
#' @param stepSize Initial stepsize of the outer iteration (theta_{k+1} = theta_k + Stepsize \* Stepdirection)
#' @param maxIter_line maximal number of iterations for line search
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @param eps_WW additional parameter to shorten the search if the upper and lower bound are already very close to each other
#' @export
exact_weakWolfeLineSearch <- function(cpptsemObject,
                                      objective,
                                      adaptiveLassoWeights,
                                      thetaNames,
                                      regIndicators,
                                      lambda,
                                      theta_kp1,
                                      m2LL_kp1,
                                      g_kp1,
                                      d,
                                      c1,
                                      c2,
                                      stepSize,
                                      maxIter_line,
                                      verbose,
                                      eps_WW){
  # Adapted from HANSO (see LEWIS, D.S. &  OVERTONNON, M. L., (2013). SMOOTH OPTIMIZATION VIA BFGS)

  # deep clone of cpptsemObject
  cpptsemObject_i <- rlang::duplicate(cpptsemObject, shallow = FALSE)

  if(!is.null(adaptiveLassoWeights)){
    adaptiveLassoWeightsMatrix <- diag(adaptiveLassoWeights)
  }else{
    adaptiveLassoWeightsMatrix <- diag(length(thetaNames))
  }
  rownames(adaptiveLassoWeightsMatrix) <- thetaNames
  colnames(adaptiveLassoWeightsMatrix) <- thetaNames

  stepSize_L <- 0 # lower bound
  stepSize_U <- Inf # upper bound

  # get penalized M2LL for step size 0:

  h_0 <- m2LL_kp1 + lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1)[regIndicators,]))

  # get (sub-)gradients for step size 0:
  g_0 <- regCtsem::exact_getSubgradients(theta = as.matrix(theta_kp1),
                                         jacobian = as.matrix(g_kp1),
                                         regIndicators = regIndicators,
                                         lambda = lambda,
                                         adaptiveLassoWeightsMatrix = adaptiveLassoWeightsMatrix)

  # Inexact Line Search
  i <- 0

  while(i<maxIter_line){
    i <- i+1 # update iterations

    if(verbose==2){
      cat(paste0("\r",
                 "## Linesearch: [Lower step size: ", sprintf('%.3f', stepSize_L, 3),
                 ", current step size: ", sprintf('%.3f', stepSize, 3),
                 ", upper step size: ", sprintf('%.3f', stepSize_U, 3),
                 "] ##"
      )
      )
      flush.console()
    }

    if(stepSize_U-stepSize_L < eps_WW){
      # The following part is adapted from https://github.com/ceres-solver/ceres-solver/blob/master/internal/ceres/line_search.cc
      # and is a workaround for when the line search algorithm does not find
      # an exact point that satisfies both Wolfe conditions and results in an upper and lower bound which are
      # approximately identical. In this case, the iterations are stopped and the current lower bound is returned
      if(verbose>1){
        warning("Could not find an exact point satisfying both Wolfe conditions.")
      }
      # an additional problem when the hessian is not well approximated is that the search direction might be wrong.
      # this would result in an increasing Likelihood which is indicated by a lower bound of 0.
      if(stepSize_L == 0){warning("Searchdirection might be problematic. Line search returned step size 0.")}
      return(stepSize_L)
    }

    # new theta
    theta_kp1_td <- rlang::duplicate(theta_kp1)+stepSize*d
    # get L(x+stepSize*d) and L'(x+stepSize*d)
    cpptsemObject_i$setParameterValues(theta_kp1_td, thetaNames)
    if(tolower(objective) == "ml"){
      invisible(capture.output(out1 <- try(cpptsemObject_i$computeRAM(), silent = TRUE), type = "message"))
      invisible(capture.output(out2 <- try(cpptsemObject_i$fitRAM(), silent = TRUE), type = "message"))
    }else{
      invisible(capture.output(out1 <- try(cpptsemObject_i$computeAndFitKalman(), silent = TRUE), type = "message"))
      out2 <- NA
    }
    if(any(class(out1)== "try-error") | any(class(out2)== "try-error") | !is.finite(cpptsemObject_i$m2LL)){
      # if the starting values are not feasible
      stepSize <- .5*stepSize
      next
    }

    if(any(class(cpptsemObject_i)=="try-error")){
      # if the starting values are far off, the gradients can be very large and testing for the initial
      # step size might result in infeasible parameter values. In this case: try smaller step size
      stepSize <- .5*stepSize
      next
    }

    m2LL_kp1_td <- cpptsemObject_i$m2LL
    invisible(capture.output(g_kp1_td <- try(exact_getCppGradients(cpptsemObject = cpptsemObject_i, objective = objective))))

    if(any(class(g_kp1_td)=="try-error")){
      # sometimes the stepsize of 1 results in NA; in this case: use smaller step size. Also, remove unfeasible starting values
      stepSize <- .9*stepSize
      next
    }
    g_kp1_td <- matrix(g_kp1_td, nrow = length(g_kp1_td), ncol = 1)
    rownames(g_kp1_td) <- thetaNames

    if(any(is.na(g_kp1_td))){
      # sometimes the stepsize of 1 results in NA; in this case: use smaller step size. Also, remove unfeasible starting values
      stepSize <- .9*stepSize
      next
    }

    # compute h(stepSize) = L(x+td) + p(x+td) - L(x) - p(x), where p(x) is the penalty function
    h_t <- m2LL_kp1_td +  lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1_td)[regIndicators,]))
    # compute h'(stepSize)
    g_t <- regCtsem::exact_getSubgradients(theta = theta_kp1_td,
                                           jacobian = g_kp1_td,
                                           regIndicators = regIndicators,
                                           lambda = lambda,
                                           adaptiveLassoWeightsMatrix = adaptiveLassoWeightsMatrix)

    # Check Armijo and Curvature Condition
    Armijo <- h_t <= h_0 + c1*stepSize*t(d)%*%g_0 # for first wolfe condition (Armijo condition)
    Curvature <- t(g_t)%*%d >= c2*t(d)%*%g_0 # for second wolfe condition (curvature condition)

    if (!Armijo){
      # if the first condition is not satisfied at stepSize (i.e., stepSize is not resulting in a sufficient decrease of the function), set the upper limit to stepSize
      stepSize_U <- stepSize
    }else if(!Curvature){
      # if stepSize results in a sufficient decrease, but the second wolfe condition is not satisfied at stepSize (i.e. the step size is too small), set the lower limit to stepSize
      stepSize_L <- stepSize
    }else{
      # if both conditions are satisfied: end line search
      return(stepSize)
    }
    if(stepSize_U < Inf){
      # if the upper limit stepSize_U has been updated, set stepSize to the mean of the interval (lower, upper)
      stepSize <- (stepSize_L+stepSize_U)/2
    }else{
      # if the upper limit has not been updated, set t to 2 times the lower limit
      stepSize <- 2*stepSize_L
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
#' @param theta_k Theta at iteration k
#' @param g_k Gradients of the likelihood function at iteration k
#' @param H_k Hessian of the likelihood function at iteration k
#' @param theta_kp1 Theta at iteration k+1
#' @param g_kp1 Gradients of the likelihood function at iteration k+1
#' @param cautious boolean: should the update be skipped if it would result in a non positive definite Hessian?
#' @param hessianEps controls when the update of the Hessian approximation is skipped
#' @export
exact_getBFGS <- function(theta_k, g_k, H_k, theta_kp1, g_kp1, cautious = TRUE, hessianEps = .001){

  y <- g_kp1-g_k
  d <- theta_kp1-theta_k

  # test if positive definiteness is ensured
  skipUpdate <- try((t(y)%*%d < hessianEps) && cautious, silent = TRUE)
  if(any(class(skipUpdate) == "try-error") || skipUpdate || is.na(skipUpdate)){
    # Hessian might become non-positive definite. Return without update
    return(H_k)
  }
  if(t(y)%*%d < 0){
    warning("Hessian update possibly non-positive definite.")
  }

  H_kp1 <- H_k - (H_k%*%d%*%t(d)%*%H_k)/as.numeric(t(d)%*%H_k%*%d) + (y%*%t(y))/as.numeric(t(y)%*%d)

  if(anyNA(H_kp1)){
    warning("Invalid Hessian. Returning previous Hessian")
    return(H_k)
  }
  HessianEigen <- eigen(H_kp1)
  iscomplex <- any(!Im(HessianEigen$values) == 0)
  if(iscomplex){
    eigenVectors <- Re(HessianEigen$vectors)
    eigenValues <- Re(HessianEigen$values)
    H_kp1 <- eigenVectors%*%diag(eigenValues)%*%t(eigenVectors)
  }
  if(any(HessianEigen$values < 0)){
    while(any(eigen(H_kp1)$values < 0)){
      H_kp1 <- H_kp1 + .01*diag(nrow(H_kp1))
    }
  }

  return(H_kp1)

}






