#### Optimization with GIST

#### Main function ####

#' exact_GIST
#'
#' General Iterative Shrinkage and Thresholding Algorithm based on Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). A General Iterative Shrinkage and Thresholding Algorithm for Non-convex Regularized Optimization Problems. In S. Dasgupta & D. McAllester (Eds.), Proceedings of Machine Learning Research (PMLR; Vol. 28, Issue 2, pp. 37--45). PMLR. http://proceedings.mlr.press
#'
#' GIST minimizes a function of form f(theta) = l(theta) + g(theta), where l is the likelihood and g is a penalty function. Various penalties are supported, however currently only lasso and adaptive lasso are implemented.
#'
#' NOTE: Function located in file GIST.R
#'
#' @param cpptsemObject Fitted object of class cpptsem
#' @param dataset only required if objective = "Kalman" and ctsemObject is of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param regIndicators Vector with names of regularized parameters
#' @param targetVector named vector with values towards which the parameters are regularized
#' @param lambdas Vector with lambda values that should be tried
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param sparseParameters labeled vector with parameter estimates of the most sparse model.
#' @param eta if the current step size fails, eta will decrease the step size. Must be > 1
#' @param sig Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 < sigma < 1
#' @param initialStepsize intial step size tested at the beginning of each optimization iteration
#' @param stepsizeMin Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param stepsizeMax Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param GISTLinesearchCriterion criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
#' @param GISTNonMonotoneNBack in case of non-monotone line search: Number of preceding regM2LL values to consider
#' @param maxIter_out maximal number of outer iterations
#' @param maxIter_in maximal number of inner iterations
#' @param break_outer Stopping criterion for outer iterations. It has to be a named value. By default (name: gradient), a relative first-order condition is checked, where the maximum absolute value of the gradients is compared to break_outer (see https://de.mathworks.com/help/optim/ug/first-order-optimality-measure.html). Alternatively, an absolute tolerance can be passed to the function (e.g., break_outer = c("gradient" = .0001)). Instead of relative gradients, the change in parameters can used as breaking criterion. To this end, use c("parameterChange" = .00001)
#' @param verbose set to 1 to print additional information and plot the convergence and 2 for further details.
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param sampleSize sample size for scaling lambda with N
#' @param approxFirst Should approximate optimization be used first to obtain start values for exact optimization?
#' @param numStart Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param controlApproxOptimizer settings passed to optimx or Rsolnp
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress. Set verbose = -1 to use a C++ implementation of GIST (not much faster which is why the easier to handle R implementation is the default)
#' @keywords internal
exact_GIST <- function(cpptsemObject, dataset, objective, regIndicators, targetVector, lambdas, adaptiveLassoWeights,
                       # additional settings
                       sparseParameters,# = NULL,
                       eta,# = 2,
                       sig,# = 10^(-5),
                       initialStepsize,# = 1,
                       stepsizeMin,# = 1/(10^30),
                       stepsizeMax,# = 10^30,
                       GISTLinesearchCriterion,# = "monotone",
                       GISTNonMonotoneNBack,# = 5,
                       maxIter_out,# = 100,
                       maxIter_in,# = 1000,
                       break_outer,# = c("parameterChange" = 10^(-5)),
                       scaleLambdaWithN,# = TRUE,
                       sampleSize,
                       approxFirst,# = F,
                       numStart,# = 3,
                       controlApproxOptimizer,
                       verbose# = 0
){
  # Setup
  # get parameter values
  initialParameters <- cpptsemObject$getParameterValues()
  thetaNames <- names(initialParameters)
  startingValues <- initialParameters

  break_crit <- names(break_outer)
  if(is.null(break_crit) || !(break_crit %in% c("parameterChange", "gradient", "fitChange"))){stop("Unknown breaking criterion. The value passed to break_crit has to be named (either 'parameterChange' or 'gradient') See ?controlGIST for details on break_outer")}
  if(is.character(break_outer)){
    if(break_crit == "parameterChange" || break_crit == "fitChange"){stop("break_outer has to be numeric if parameterChange or fitChange is used as criterion")}
    break_outer <- eval(parse(text = break_outer))
    names(break_outer) <- break_crit
    if(verbose > 0){message(paste0("Breaking criterion: max(abs(gradient)) < ", break_outer))}
  }

  thetas <- matrix(NA,
                   nrow = length(initialParameters),
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
  retryOnce <- TRUE # will retry optimizing once with new starting values
  while(iteration <= numLambdas){
    lambda <- lambdas[iteration]

    # update progress bar
    setTxtProgressBar(pbar, iteration)

    lambda = ifelse(scaleLambdaWithN, lambda*sampleSize, lambda) # set lambda*samplesize

    # should the results first be approximated?
    if(approxFirst){
      startingValues <- exact_tryStartingValues(startingValues = startingValues,
                                                returnAs = "vector",
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
    # use Gist
    if(is.null(targetVector) || all(targetVector == 0)){
      resGIST <- try(GIST(cpptsemObject = cpptsemObject,
                                    startingValues = startingValues,
                                    objective = objective,
                                    lambda = lambda,
                                    adaptiveLassoWeights = adaptiveLassoWeights,
                                    regularizedParameters = regIndicators,
                                    eta = eta,
                                    sig = sig,
                                    initialStepsize = initialStepsize,
                                    stepsizeMin = stepsizeMin,
                                    stepsizeMax = stepsizeMax,
                                    GISTLinesearchCriterion = GISTLinesearchCriterion,
                                    GISTNonMonotoneNBack = GISTNonMonotoneNBack,
                                    maxIter_out = maxIter_out,
                                    maxIter_in = maxIter_in,
                                    break_outer = break_outer,
                                    verbose = verbose,
                                    silent = FALSE))
    }else{
      resGIST <- try(GISTWithTarget(cpptsemObject = cpptsemObject, startingValues = startingValues,
                                              objective = objective, lambda = lambda, adaptiveLassoWeights = adaptiveLassoWeights,
                                              regularizedParameters = regIndicators, targetVector = targetVector,
                                              eta = eta, sig = sig, initialStepsize = initialStepsize, stepsizeMin = stepsizeMin, stepsizeMax = stepsizeMax,
                                              GISTLinesearchCriterion = GISTLinesearchCriterion, GISTNonMonotoneNBack = GISTNonMonotoneNBack,
                                              maxIter_out = maxIter_out, maxIter_in = maxIter_in,
                                              break_outer = break_outer, verbose = verbose, silent = FALSE))
    }


    if(any(class(resGIST) == "try-error") || !resGIST$convergence){

      if(retryOnce){
        warning(paste("Model for lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda), " did not converge. Retrying with different starting values.", sep = ""))

      }else{
        warning(paste("Model for lambda = ", ifelse(scaleLambdaWithN, lambda/sampleSize, lambda), " did not converge.", sep = ""))
      }
      # reset theta for next iteration
      startingValues <- initialParameters

      if(retryOnce){
        retryOnce <- FALSE
        next
      }
      retryOnce <- TRUE

    }else{
      if(resGIST$type == "cpptsem"){
        newValues <- resGIST$model$getParameterValues()
        startingValues <- newValues
        # save fit
        cM2LL <- ifelse(resGIST$convergence, resGIST$model$m2LL, Inf)
        if(is.null(targetVector)){
          cRegM2LL <- ifelse(resGIST$convergence, cM2LL +  exact_getPenaltyValue(lambda = lambda,
                                                                                           theta = newValues,
                                                                                           regIndicators = regIndicators,
                                                                                           adaptiveLassoWeights = adaptiveLassoWeights), Inf)
        }else{
          cRegM2LL <- ifelse(resGIST$convergence, cM2LL +  exact_getPenaltyValueWithTarget(lambda = lambda,
                                                                                                     theta = newValues,
                                                                                                     regIndicators = regIndicators,
                                                                                                     targetVector = targetVector,
                                                                                                     adaptiveLassoWeights = adaptiveLassoWeights), Inf)
        }
      }else{
        stop("Internal error in GIST.")
      }
    }

    # save results
    thetas[,iteration] <- newValues[rownames(thetas)]
    m2LL[iteration] <- cM2LL
    regM2LL[iteration] <- cRegM2LL

    iteration <- iteration + 1
  }

  return(list("lambdas" = lambdas, "thetas" = thetas, "m2LL" = m2LL,"regM2LL" = regM2LL))

}


#' GIST
#'
#' General Iterative Shrinkage and Thresholding Algorithm based on Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). A General Iterative Shrinkage and Thresholding Algorithm for Non-convex Regularized Optimization Problems. In S. Dasgupta & D. McAllester (Eds.), Proceedings of Machine Learning Research (PMLR; Vol. 28, Issue 2, pp. 37--45). PMLR. http://proceedings.mlr.press
#'
#' GIST minimizes a function of form f(theta) = l(theta) + g(theta), where l is the likelihood and g is a penalty function. Various penalties are supported, however currently only lasso and adaptive lasso are implemented.
#'
#' NOTE: Function located in file GIST.R
#'
#' @param cpptsemObject model of type cpptsem
#' @param startingValues named vector with starting values
#' @param objective "ML" for maximum likelihood SEM, "Kalman" for Kalman filter
#' @param lambda penalty value
#' @param adaptiveLassoWeights named vector with adaptive lasso weights
#' @param regularizedParameters named vector of regularized parameters
#' @param eta if the current step size fails, eta will decrease the step size. Must be > 1
#' @param sig Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 < sigma < 1
#' @param initialStepsize initial stepsize to be tried in the outer iteration
#' @param stepsizeMin Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param stepsizeMax Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param GISTLinesearchCriterion criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
#' @param GISTNonMonotoneNBack in case of non-monotone line search: Number of preceding regM2LL values to consider
#' @param maxIter_out maximal number of outer iterations
#' @param maxIter_in maximal number of inner iterations
#' @param break_outer Stopping criterion for outer iterations. It has to be a named value. By default (name: gradient), a relative first-order condition is checked, where the maximum absolute value of the gradients is compared to break_outer (see https://de.mathworks.com/help/optim/ug/first-order-optimality-measure.html). Alternatively, an absolute tolerance can be passed to the function (e.g., break_outer = c("gradient" = .0001)). Instead of relative gradients, the change in parameters can used as breaking criterion. To this end, use c("parameterChange" = .00001)
#' @param verbose set to 1 to print additional information and plot the convergence and 2 for further details.
#' @param silent suppress all warning messages
#' @keywords internal
GIST <- function(cpptsemObject,
                 startingValues,
                 objective,
                 lambda,
                 adaptiveLassoWeights,
                 regularizedParameters,
                 eta,# = 2,
                 sig,# = 10^(-5),
                 initialStepsize,# = 1,
                 stepsizeMin,# = 1/(10^30),
                 stepsizeMax,# = 10^30,
                 GISTLinesearchCriterion,# = "monotone",
                 GISTNonMonotoneNBack,# = 5,
                 maxIter_out,# = 100,
                 maxIter_in,# = 1000,
                 break_outer,# = c("parameterChange" = 10^(-5)),
                 verbose,# = 0,
                 silent # = FALSE
){
  break_crit <- names(break_outer)
  # iteration counter
  k_out <- 1
  convergence <- TRUE

  parameterNames <- names(startingValues)
  adaptiveLassoWeightsMatrix <- diag(adaptiveLassoWeights[parameterNames])
  rownames(adaptiveLassoWeightsMatrix) <- parameterNames
  colnames(adaptiveLassoWeightsMatrix) <- parameterNames

  # set parameters
  cpptsemObject$setParameterValues(startingValues, names(startingValues))
  if(tolower(objective) == "ml"){
    invisible(capture.output(out1 <- try(cpptsemObject$computeRAM(), silent = T), type = "message"))
    invisible(capture.output(out2 <- try(cpptsemObject$fitRAM(), silent = T), type = "message"))
    out3 <- exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective)
  }else{
    invisible(capture.output(out1 <- try(cpptsemObject$computeAndFitKalman(0), silent = TRUE), type = "message"))
    out2 <- NA
    out3 <- exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective)
  }
  if(any(class(out1) == "try-error")  |
     any(class(out2) == "try-error") |
     any(class(out3) == "try-error") |
     anyNA(out3)){
    stop("Infeasible starting values in GIST.")
  }

  parameters_km1 <- NULL
  parameters_k <- cpptsemObject$getParameterValues()
  parameterNames <- names(parameters_k)
  gradients_km1 <- NULL
  gradients_k <- out3
  m2LL_k <- cpptsemObject$m2LL

  regM2LL_k <- m2LL_k + exact_getPenaltyValue(lambda = lambda,
                                                        theta = parameters_k,
                                                        regIndicators = regularizedParameters,
                                                        adaptiveLassoWeights = adaptiveLassoWeights)
  regM2LL <- rep(NA, maxIter_out)
  regM2LL[1] <- regM2LL_k

  if(verbose == 2){
    convergencePlotValues <- matrix(NA, nrow = length(parameters_k), ncol = maxIter_out, dimnames = list(parameterNames, 1:maxIter_out))
  }

  resetStepSize <- TRUE
  #resetIteration <- -1
  while(k_out < maxIter_out){
    k <- 0
    # set initial step size
    if(is.null(parameters_km1)){
      stepsize <- initialStepsize
    }else{
      x_k <- parameters_k - parameters_km1
      y_k <- gradients_k - gradients_km1

      stepsize <- (t(x_k)%*%y_k)/(t(x_k)%*%x_k)
      # sometimes this step size is extremely large and the algorithm converges very slowly
      # we found that in these cases it can help to reset the stepsize
      # this will be done randomly here:
      if(runif(1,0,1) > .9){
        stepsize <- initialStepsize
      }

      if(is.na(stepsize)){
        if(!silent){
          warning(paste0("Outer iteration ", k_out, ": NA or infinite step size..."))}
        break
      }
      if(is.infinite(stepsize)){
        if(!silent){
          warning(paste0("Outer iteration ", k_out, ": NA or infinite step size..."))}
        break
      }
      if(stepsize < stepsizeMin){

        if(!resetStepSize){
          if(!silent){warning(paste0("Outer iteration ", k_out, ": Stepsize below specified minimum..."))}
          break
        }
        stepsize <- initialStepsize
      }
      if(stepsize > stepsizeMax){
        if(!resetStepSize){
          if(!silent){warning(paste0("Outer iteration ", k_out, ": Stepsize above specified maximum..."))}
          break
        }
        stepsize <- initialStepsize
      }

    }


    # inner iteration
    while(k < maxIter_in){
      u_k <- parameters_k - gradients_k/as.vector(stepsize)
      parameters_kp1 <- rep(NA, length(parameters_k))
      names(parameters_kp1) <- parameterNames
      for(i in seq_len(length(parameters_kp1))){
        parameterName <- parameterNames[i]
        lambda_i <- lambda*adaptiveLassoWeights[parameterName]
        if(parameterName %in% regularizedParameters){
          # update parameter i with lasso
          parameters_kp1[parameterName] <- sign(u_k[parameterName])*max(c(0,abs(u_k[parameterName]) - lambda_i/stepsize))
        }else{
          parameters_kp1[parameterName] <- u_k[parameterName]
        }
      }

      # update parameters
      cpptsemObject$setParameterValues(parameters_kp1, names(parameters_kp1))
      # compute likelihood
      if(tolower(objective) == "ml"){
        invisible(capture.output(out1 <- try(cpptsemObject$computeRAM(), silent = T), type = "message"))
        invisible(capture.output(out2 <- try(cpptsemObject$fitRAM(), silent = T), type = "message"))
      }else{
        invisible(capture.output(out1 <- try(cpptsemObject$computeAndFitKalman(0), silent = TRUE), type = "message"))
        out2 <- NA
      }
      if(any(class(out1) == "try-error")  |
         any(class(out2) == "try-error")){

        # update step size
        stepsize <- eta*stepsize

        # update iteration counter
        k <- k+1

        # skip rest
        next
      }

      m2LL_kp1 <- cpptsemObject$m2LL

      if(is.na(m2LL_kp1) | is.infinite(m2LL_kp1)){

        # update step size
        stepsize <- eta*stepsize

        # update iteration counter
        k <- k+1

        # skip rest
        next
      }
      regM2LL_kp1 <- m2LL_kp1 + exact_getPenaltyValue(lambda = lambda,
                                                                theta = parameters_kp1,
                                                                regIndicators = regularizedParameters,
                                                                adaptiveLassoWeights = adaptiveLassoWeights)

      if(verbose == 2){
        cat(paste0("\r",
                   "###### [", sprintf("%*d", 3, k_out), "-", sprintf("%*d", 3, k),
                   " | regM2LL:  ", sprintf('%.3f',regM2LL_kp1),
                   " | zeroed: ", sprintf("%*d", 3, sum(parameters_kp1[regularizedParameters] == 0)),
                   " | stepsize: ", sprintf("%.3f", stepsize),
                   " ######"
        )
        )
      }

      # break if line search condition is satisfied
      if(GISTLinesearchCriterion == "monotone"){
        breakCriterion <- regM2LL_kp1 <= regM2LL_k - (sig/2) * stepsize * sum((parameters_k - parameters_kp1)^2)
      }else if(GISTLinesearchCriterion == "non-monotone"){
        nBack <- max(1,k_out-GISTNonMonotoneNBack)
        breakCriterion <- regM2LL_kp1 <= max(regM2LL[nBack:k_out]) - (sig/2) * stepsize * sum((parameters_k - parameters_kp1)^2)
      }else{
        stop("Unknown GISTLinesearchCriterion. Possible are monotone and non-monotone.")
      }
      if(breakCriterion){
        #print("breaking inner")
        break
      }

      # update step size
      stepsize <- eta*stepsize

      # update iteration counter
      k <- k+1
    }

    if(k == maxIter_in){
      if(!silent){
        warning("Maximal number of inner iterations used by GIST. Consider increasing the number of inner iterations.")}
    }

    if(tolower(objective) == "ml"){
      out3 <- exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective)
    }else{
      out3 <- exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective)
    }
    if(any(class(out3) == "try-error") |
       any(is.na(out3))){
      if(!silent){
        convergence <- FALSE
        stop("No gradients in GIST")}
    }

    gradients_kp1 <- out3


    # update parameters for next iteration
    parameters_km1 <- parameters_k
    parameters_k <- parameters_kp1
    gradients_km1 <- gradients_k
    gradients_k <- gradients_kp1
    m2LL_k <- m2LL_kp1
    regM2LL_k <- regM2LL_kp1

    k_out <- k_out + 1
    regM2LL[k_out] <- regM2LL_k


    # break outer loop if stopping criterion is satisfied
    if(break_crit == "gradient"){

      # Gradient based break condition: Problem: gradients extremely dependent on the epsilon in the approximation
      parameterMatrix <- matrix(parameters_k[parameterNames], ncol = 1)
      rownames(parameterMatrix) <- parameterNames
      gradientMatrix <- matrix(gradients_k[parameterNames], ncol = 1)
      rownames(gradientMatrix) <- parameterNames
      subgradients <- exact_getSubgradients(theta = parameterMatrix, jacobian = gradientMatrix,
                                            regIndicators = regularizedParameters, lambda = lambda,
                                            adaptiveLassoWeightsMatrix = adaptiveLassoWeightsMatrix)

      breakOuter <- max(abs(subgradients)) < break_outer

      if(verbose == 1){
        plot(x=1:maxIter_out, y = regM2LL, xlab = "iteration", ylab = "f(theta)", type = "l", main = "Convergence Plot")
        cat(paste0("\r",
                   "## [", sprintf("%*d", 3, k_out),
                   "] m2LL: ", sprintf('%.3f',m2LL_k),
                   " | regM2LL:  ", sprintf('%.3f',regM2LL_k),
                   " | zeroed: ", sprintf("%*d", 3, sum(parameters_k[regularizedParameters] == 0)),
                   " | max(abs(grad)): ", sprintf("%.3f", max(abs(subgradients))),
                   " ##"
        )
        )
      }

      if(verbose == 2){
        convergencePlotValues[,k_out] <- parameters_k
        matplot(x=1:maxIter_out, y = t(convergencePlotValues), xlab = "iteration", ylab = "value", type = "l", main = "Convergence Plot")
      }

    }else if(break_crit == "parameterChange"){

      breakOuter <- (sqrt(sum((parameters_k - parameters_km1)^2))/sqrt(sum((parameters_km1)^2))) < break_outer

      if(verbose == 1){
        plot(x=1:maxIter_out, y = regM2LL, xlab = "iteration", ylab = "f(theta)", type = "l", main = "Convergence Plot")
        cat(paste0("\r",
                   "## [", sprintf("%*d", 3, k_out),
                   "] m2LL: ", sprintf('%.3f',m2LL_k),
                   " | regM2LL:  ", sprintf('%.3f',regM2LL_k),
                   " | zeroed: ", sprintf("%*d", 3, sum(parameters_k[regularizedParameters] == 0)),
                   " ##"
        )
        )
      }

      if(verbose == 2){
        convergencePlotValues[,k_out] <- parameters_k
        matplot(x=1:maxIter_out, y = t(convergencePlotValues), xlab = "iteration", ylab = "value", type = "l", main = "Convergence Plot")
      }

    }else if(break_crit == "fitChange"){
      breakOuter <- abs(regM2LL[k_out] - regM2LL[k_out-1]) < break_outer

      if(verbose == 1){
        plot(x=1:maxIter_out, y = regM2LL, xlab = "iteration", ylab = "f(theta)", type = "l", main = "Convergence Plot")
        cat(paste0("\r",
                   "## [", sprintf("%*d", 3, k_out),
                   "] m2LL: ", sprintf('%.3f',m2LL_k),
                   " | regM2LL:  ", sprintf('%.3f',regM2LL_k),
                   " | zeroed: ", sprintf("%*d", 3, sum(parameters_k[regularizedParameters] == 0)),
                   " | fitChange: ", sprintf("%.5f", abs(regM2LL[k_out] - regM2LL[k_out-1])),
                   " ##"
        )
        )
      }

      if(verbose == 2){
        convergencePlotValues[,k_out] <- parameters_k
        matplot(x=1:maxIter_out, y = t(convergencePlotValues), xlab = "iteration", ylab = "value", type = "l", main = "Convergence Plot")
      }


    }else{
      stop("Unknown breaking criterion. The value passed to break_crit has to be named (either 'parameterChange' or 'gradient') See ?controlGIST for details on break_outer")
    }

    if(breakOuter){
      break
    }

  }
  if(is.na(regM2LL_k) | is.infinite(regM2LL_k) |
     anyNA(parameters_k) | any(is.infinite(parameters_k)) |
     anyNA(gradients_k) | any(is.infinite(gradients_k))){
    convergence <- FALSE
  }

  if(k_out == maxIter_out){
    if(!silent){
      warning("Maximal number of outer iterations used by GIST. Consider increasing the number of outer iterations.")}
  }


  return(list("type" = "cpptsem",
              "model" = cpptsemObject,
              "regM2LL" = regM2LL,
              "convergence" = convergence))
}


#' GISTWithTarget
#'
#' General Iterative Shrinkage and Thresholding Algorithm based on Gong, P., Zhang, C., Lu, Z., Huang, J., & Ye, J. (2013). A General Iterative Shrinkage and Thresholding Algorithm for Non-convex Regularized Optimization Problems. In S. Dasgupta & D. McAllester (Eds.), Proceedings of Machine Learning Research (PMLR; Vol. 28, Issue 2, pp. 37--45). PMLR. http://proceedings.mlr.press
#'
#' GIST minimizes a function of form f(theta) = l(theta) + g(theta), where l is the likelihood and g is a penalty function. Various penalties are supported, however currently only lasso and adaptive lasso are implemented.
#'
#' NOTE: Function located in file GIST.R
#'
#' @param cpptsemObject model of type MxObject to compute the gradients
#' @param startingValues named vector with starting values
#' @param objective "ML" for maximum likelihood SEM, "Kalman" for Kalman filter
#' @param lambda penalty value
#' @param adaptiveLassoWeights named vector with adaptive lasso weights
#' @param regularizedParameters named vector of regularized parameters
#' @param targetVector named vector with values towards which the parameters are regularized
#' @param eta if the current step size fails, eta will decrease the step size. Must be > 1
#' @param sig Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 < sigma < 1
#' @param initialStepsize initial stepsize to be tried in the outer iteration
#' @param stepsizeMin Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param stepsizeMax Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param GISTLinesearchCriterion criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
#' @param GISTNonMonotoneNBack in case of non-monotone line search: Number of preceding regM2LL values to consider
#' @param maxIter_out maximal number of outer iterations
#' @param maxIter_in maximal number of inner iterations
#' @param break_outer Stopping criterion for outer iterations. It has to be a named value. By default (name: gradient), a relative first-order condition is checked, where the maximum absolute value of the gradients is compared to break_outer (see https://de.mathworks.com/help/optim/ug/first-order-optimality-measure.html). Alternatively, an absolute tolerance can be passed to the function (e.g., break_outer = c("gradient" = .0001)). Instead of relative gradients, the change in parameters can used as breaking criterion. To this end, use c("parameterChange" = .00001)
#' @param verbose set to 1 to print additional information and plot the convergence and 2 for further details.
#' @param silent suppress all warning messages
#' @keywords internal
GISTWithTarget <- function(cpptsemObject,
                           startingValues,
                           objective,
                           lambda,
                           adaptiveLassoWeights,
                           regularizedParameters,
                           targetVector,
                           eta,# = 2,
                           sig,# = 10^(-5),
                           initialStepsize,# = 1,
                           stepsizeMin,# = 1/(10^30),
                           stepsizeMax,# = 10^30,
                           GISTLinesearchCriterion,# = "monotone",
                           GISTNonMonotoneNBack,# = 5,
                           maxIter_out,# = 100,
                           maxIter_in,# = 1000,
                           break_outer,# = c("parameterChange" = 10^(-5)),
                           verbose,# = 0,
                           silent # = FALSE
){
  break_crit <- names(break_outer)
  # iteration counter
  k_out <- 1
  convergence <- TRUE

  parameterNames <- names(startingValues)
  adaptiveLassoWeightsMatrix <- diag(adaptiveLassoWeights[parameterNames])
  rownames(adaptiveLassoWeightsMatrix) <- parameterNames
  colnames(adaptiveLassoWeightsMatrix) <- parameterNames

  # set parameters
  cpptsemObject$setParameterValues(startingValues, names(startingValues))
  if(tolower(objective) == "ml"){
    invisible(capture.output(out1 <- try(cpptsemObject$computeRAM(), silent = T), type = "message"))
    invisible(capture.output(out2 <- try(cpptsemObject$fitRAM(), silent = T), type = "message"))
    out3 <- exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective)
  }else{
    invisible(capture.output(out1 <- try(cpptsemObject$computeAndFitKalman(0), silent = TRUE), type = "message"))
    out2 <- NA
    out3 <- exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective)
  }
  if(any(class(out1) == "try-error")  |
     any(class(out2) == "try-error") |
     any(class(out3) == "try-error") |
     anyNA(out3)){
    stop("Infeasible starting values in GIST.")
  }

  parameters_km1 <- NULL
  parameters_k <- cpptsemObject$getParameterValues()
  parameterNames <- names(parameters_k)
  gradients_km1 <- NULL
  gradients_k <- out3
  m2LL_k <- cpptsemObject$m2LL

  regM2LL_k <- m2LL_k + exact_getPenaltyValueWithTarget(lambda = lambda,
                                                                  theta = parameters_k,
                                                                  regIndicators = regularizedParameters,
                                                                  targetVector = targetVector,
                                                                  adaptiveLassoWeights = adaptiveLassoWeights)
  regM2LL <- rep(NA, maxIter_out)
  regM2LL[1] <- regM2LL_k

  if(verbose == 2){
    convergencePlotValues <- matrix(NA, nrow = length(parameters_k), ncol = maxIter_out, dimnames = list(parameterNames, 1:maxIter_out))
  }

  resetStepSize <- TRUE
  #resetIteration <- -1
  while(k_out < maxIter_out){
    k <- 0
    # set initial step size
    if(is.null(parameters_km1)){
      stepsize <- initialStepsize
    }else{
      x_k <- parameters_k - parameters_km1
      y_k <- gradients_k - gradients_km1

      stepsize <- (t(x_k)%*%y_k)/(t(x_k)%*%x_k)
      # sometimes this step size is extremely large and the algorithm converges very slowly
      # we found that in these cases it can help to reset the stepsize
      # this will be done randomly here:
      if(runif(1,0,1) > .9){
        stepsize <- initialStepsize
      }

      if(is.na(stepsize)){
        if(!silent){
          warning(paste0("Outer iteration ", k_out, ": NA or infinite step size..."))}
        break
      }
      if(is.infinite(stepsize)){
        if(!silent){
          warning(paste0("Outer iteration ", k_out, ": NA or infinite step size..."))}
        break
      }
      if(stepsize < stepsizeMin){

        if(!resetStepSize){
          if(!silent){warning(paste0("Outer iteration ", k_out, ": Stepsize below specified minimum..."))}
          break
        }
        stepsize <- initialStepsize
      }
      if(stepsize > stepsizeMax){
        if(!resetStepSize){
          if(!silent){warning(paste0("Outer iteration ", k_out, ": Stepsize above specified maximum..."))}
          break
        }
        stepsize <- initialStepsize
      }

    }


    # inner iteration
    while(k < maxIter_in){
      u_k <- parameters_k - gradients_k/as.vector(stepsize)
      parameters_kp1 <- rep(NA, length(parameters_k))
      names(parameters_kp1) <- parameterNames
      for(i in seq_len(length(parameters_kp1))){
        parameterName <- parameterNames[i]
        lambda_i <- lambda*adaptiveLassoWeights[parameterName]
        if(parameterName %in% regularizedParameters){
          # update parameter i with lasso
          if(u_k[parameterName] + lambda_i/stepsize < targetVector[parameterName]){
            # CASE 1
            parameters_kp1[parameterName] <- u_k[parameterName] + lambda_i/stepsize
          }else if(u_k[parameterName] - lambda_i/stepsize > targetVector[parameterName]){
            # CASE 2
            parameters_kp1[parameterName] <- u_k[parameterName] - lambda_i/stepsize
          }else{
            parameters_kp1[parameterName] <- targetVector[parameterName]
          }

        }else{
          parameters_kp1[parameterName] <- u_k[parameterName]
        }
      }

      # update parameters
      cpptsemObject$setParameterValues(parameters_kp1, names(parameters_kp1))
      # compute likelihood
      if(tolower(objective) == "ml"){
        invisible(capture.output(out1 <- try(cpptsemObject$computeRAM(), silent = T), type = "message"))
        invisible(capture.output(out2 <- try(cpptsemObject$fitRAM(), silent = T), type = "message"))
      }else{
        invisible(capture.output(out1 <- try(cpptsemObject$computeAndFitKalman(0), silent = TRUE), type = "message"))
        out2 <- NA
      }
      if(any(class(out1) == "try-error")  |
         any(class(out2) == "try-error")){

        # update step size
        stepsize <- eta*stepsize

        # update iteration counter
        k <- k+1

        # skip rest
        next
      }

      m2LL_kp1 <- cpptsemObject$m2LL

      if(is.na(m2LL_kp1) | is.infinite(m2LL_kp1)){

        # update step size
        stepsize <- eta*stepsize

        # update iteration counter
        k <- k+1

        # skip rest
        next
      }

      regM2LL_kp1 <- m2LL_kp1 + exact_getPenaltyValueWithTarget(lambda = lambda,
                                                                          theta = parameters_kp1,
                                                                          regIndicators = regularizedParameters,
                                                                          targetVector = targetVector,
                                                                          adaptiveLassoWeights = adaptiveLassoWeights)

      if(verbose == 2){
        cat(paste0("\r",
                   "###### [", sprintf("%*d", 3, k_out), "-", sprintf("%*d", 3, k),
                   " | regM2LL:  ", sprintf('%.3f',regM2LL_kp1),
                   " | on target: ", sprintf("%*d", 3, sum(parameters_k[regularizedParameters] - targetVector[regularizedParameters] == 0)),
                   " | stepsize: ", sprintf("%.3f", stepsize),
                   " ######"
        )
        )
      }

      # break if line search condition is satisfied
      if(GISTLinesearchCriterion == "monotone"){
        breakCriterion <- regM2LL_kp1 <= regM2LL_k - (sig/2) * stepsize * sum((parameters_k - parameters_kp1)^2)
      }else if(GISTLinesearchCriterion == "non-monotone"){
        nBack <- max(1,k_out-GISTNonMonotoneNBack)
        breakCriterion <- regM2LL_kp1 <= max(regM2LL[nBack:k_out]) - (sig/2) * stepsize * sum((parameters_k - parameters_kp1)^2)
      }else{
        stop("Unknown GISTLinesearchCriterion. Possible are monotone and non-monotone.")
      }
      if(breakCriterion){
        #print("breaking inner")
        break
      }

      # update step size
      stepsize <- eta*stepsize

      # update iteration counter
      k <- k+1
    }

    if(k == maxIter_in){
      if(!silent){
        warning("Maximal number of inner iterations used by GIST. Consider increasing the number of inner iterations.")}
    }

    if(tolower(objective) == "ml"){
      out3 <- exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective)
    }else{
      out3 <- exact_getCppGradients(cpptsemObject = cpptsemObject, objective = objective)
    }
    if(any(class(out3) == "try-error") |
       any(is.na(out3))){
      if(!silent){
        convergence <- FALSE
        stop("No gradients in GIST")}
    }

    gradients_kp1 <- out3


    # update parameters for next iteration
    parameters_km1 <- parameters_k
    parameters_k <- parameters_kp1
    gradients_km1 <- gradients_k
    gradients_k <- gradients_kp1
    m2LL_k <- m2LL_kp1
    regM2LL_k <- regM2LL_kp1

    k_out <- k_out + 1
    regM2LL[k_out] <- regM2LL_k


    # break outer loop if stopping criterion is satisfied
    if(break_crit == "gradient"){
      stop("Gradient based stopping not yet implemented for target optimization")
      # Gradient based break condition: Problem: gradients extremely dependent on the epsilon in the approximation
      parameterMatrix <- matrix(parameters_k[parameterNames], ncol = 1)
      rownames(parameterMatrix) <- parameterNames
      gradientMatrix <- matrix(gradients_k[parameterNames], ncol = 1)
      rownames(gradientMatrix) <- parameterNames
      subgradients <- exact_getSubgradients(theta = parameterMatrix, jacobian = gradientMatrix,
                                            regIndicators = regularizedParameters, lambda = lambda,
                                            adaptiveLassoWeightsMatrix = adaptiveLassoWeightsMatrix)

      breakOuter <- max(abs(subgradients)) < break_outer

      if(verbose == 1){
        plot(x=1:maxIter_out, y = regM2LL, xlab = "iteration", ylab = "f(theta)", type = "l", main = "Convergence Plot")
        cat(paste0("\r",
                   "## [", sprintf("%*d", 3, k_out),
                   "] m2LL: ", sprintf('%.3f',m2LL_k),
                   " | regM2LL:  ", sprintf('%.3f',regM2LL_k),
                   " | on target: ", sprintf("%*d", 3, sum(parameters_k[regularizedParameters] - targetVector[regularizedParameters] == 0)),
                   " | max(abs(grad)): ", sprintf("%.3f", max(abs(subgradients))),
                   " ##"
        )
        )
      }

      if(verbose == 2){
        convergencePlotValues[,k_out] <- parameters_k
        matplot(x=1:maxIter_out, y = t(convergencePlotValues), xlab = "iteration", ylab = "value", type = "l", main = "Convergence Plot")
      }

    }else if(break_crit == "parameterChange"){

      breakOuter <- (sqrt(sum((parameters_k - parameters_km1)^2))/sqrt(sum((parameters_km1)^2))) < break_outer

      if(verbose == 1){
        plot(x=1:maxIter_out, y = regM2LL, xlab = "iteration", ylab = "f(theta)", type = "l", main = "Convergence Plot")
        cat(paste0("\r",
                   "## [", sprintf("%*d", 3, k_out),
                   "] m2LL: ", sprintf('%.3f',m2LL_k),
                   " | regM2LL:  ", sprintf('%.3f',regM2LL_k),
                   " | on target: ", sprintf("%*d", 3, sum(parameters_k[regularizedParameters] - targetVector[regularizedParameters] == 0)),
                   " ##"
        )
        )
      }

      if(verbose == 2){
        convergencePlotValues[,k_out] <- parameters_k
        matplot(x=1:maxIter_out, y = t(convergencePlotValues), xlab = "iteration", ylab = "value", type = "l", main = "Convergence Plot")
      }

    }else if(break_crit == "fitChange"){
      breakOuter <- abs(regM2LL[k_out] - regM2LL[k_out-1]) < break_outer

      if(verbose == 1){
        plot(x=1:maxIter_out, y = regM2LL, xlab = "iteration", ylab = "f(theta)", type = "l", main = "Convergence Plot")
        cat(paste0("\r",
                   "## [", sprintf("%*d", 3, k_out),
                   "] m2LL: ", sprintf('%.3f',m2LL_k),
                   " | regM2LL:  ", sprintf('%.3f',regM2LL_k),
                   " | on target: ", sprintf("%*d", 3, sum(parameters_k[regularizedParameters] - targetVector[regularizedParameters] == 0)),
                   " | fitChange: ", sprintf("%.5f", abs(regM2LL[k_out] - regM2LL[k_out-1])),
                   " ##"
        )
        )
      }

      if(verbose == 2){
        convergencePlotValues[,k_out] <- parameters_k
        matplot(x=1:maxIter_out, y = t(convergencePlotValues), xlab = "iteration", ylab = "value", type = "l", main = "Convergence Plot")
      }


    }else{
      stop("Unknown breaking criterion. The value passed to break_crit has to be named (either 'parameterChange' or 'gradient') See ?controlGIST for details on break_outer")
    }

    if(breakOuter){
      break
    }

  }
  if(is.na(regM2LL_k) | is.infinite(regM2LL_k) |
     anyNA(parameters_k) | any(is.infinite(parameters_k)) |
     anyNA(gradients_k) | any(is.infinite(gradients_k))){
    convergence <- FALSE
  }

  if(k_out == maxIter_out){
    if(!silent){
      warning("Maximal number of outer iterations used by GIST. Consider increasing the number of outer iterations.")}
  }

  return(list("type" = "cpptsem",
              "model" = cpptsemObject,
              "regM2LL" = regM2LL,
              "convergence" = convergence))
}
