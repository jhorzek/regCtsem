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
#' @param ctsemObject if objective = "ML": Fitted object of class ctsem. If you want to use objective = "Kalman", pass an object of type ctsemInit from ctModel
#' @param mxObject Object of type MxModel
#' @param dataset only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators Vector with names of regularized parameters
#' @param lambdas Vector with lambda values that should be tried
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = theta_k + stepSize \* Stepdirection)
#' @param tryCpptsem should regCtsem try to translate the model to cpptsem? This can speed up the computation considerably but might fail for some models
#' @param forceCpptsem should cpptsem be enforced even if results differ from ctsem? Sometimes differences between cpptsem and ctsem can result from problems with numerical precision which will lead to the m,atrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "cpptsem") for more details
#' @param eta if the current step size fails, eta will decrease the step size. Must be > 1
#' @param sig GIST: sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
#' @param stepsizeMin the initial step size is an integer randomly selected between stepsizeMin and stepsizeMax. All subsequent step sizes will be computed as described by Gong et al. (2013)
#' @param stepsizeMax the initial step size is an integer randomly selected between stepsizeMin and stepsizeMax. All subsequent step sizes will be computed as described by Gong et al. (2013)
#' @param GISTLinesearchCriterion criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
#' @param GISTNonMonotoneNBack in case of non-monotone line search: Number of preceding regM2LL values to consider
#' @param maxIter_out maximal number of outer iterations
#' @param maxIter_in maximal number of inner iterations
#' @param break_outer stopping criterion for the outer iteration.
#' @param verbose set to 1 to print additional information and plot the convergence
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param sampleSize sample size for scaling lambda with N
#' @param approxFirst Should approximate optimization be used first to obtain start values for exact optimization? 1 = only for first lambda, 2 = for all lambdas
#' @param numStart Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param approxOpt Used if approxFirst = 3. Should each of the generated starting values be optimized slightly? This can substantially improve the fit of the generated starting values. 1 = optimization with optim, 2 = optimization with Rsolnp
#' @param approxMaxIt Used if approxFirst = 3 and approxOpt > 1. How many outer iterations should be given to each starting values vector? More will improve the selected starting values but slow down the computation. If approxFirst =  4, or approxFirst = 5 this will control the number of outer iteration in optim or solnp .
#' @param extraTries number of extra tries in mxTryHard for warm start
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @export
exact_GIST <- function(ctsemObject, mxObject, dataset, objective, regOn = "DRIFT", regIndicators, lambdas, adaptiveLassoWeights,
                       # additional settings
                       sparseParameters = NULL,
                       tryCpptsem, forceCpptsem = FALSE, eta = 1.5, sig = .2, initialStepsize = 1, stepsizeMin = 0, stepsizeMax = 999999999,
                       GISTLinesearchCriterion = "monotone", GISTNonMonotoneNBack = 5,
                       maxIter_out = 100, maxIter_in = 100,
                       break_outer = .00000001,
                       scaleLambdaWithN = TRUE, sampleSize, approxFirst = 0,
                       numStart = 10, approxOpt = T, approxMaxIt = 5,
                       extraTries = 3, differenceApprox = "central", verbose = 0,
                       progressBar = TRUE, parallelProgressBar = NULL){
  # Setup
  # get parameter values
  initialParameters <- OpenMx::omxGetParameters(mxObject)
  thetaNames <- names(initialParameters)
  startingValues <- initialParameters
  gradientModel <- NULL

  break_crit <- names(break_outer)
  if(is.null(break_crit) || !(break_crit %in% c("parameterChange", "gradient"))){stop("Unknown breaking criterion. The value passed to break_crit has to be named (either 'parameterChange' or 'gradient') See ?controlGIST for details on break_outer")}
  if(is.character(break_outer)){
    if(break_crit == "parameterChange"){stop("break_outer has to be numeric if parameterChange is used as criterion")}
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

  # define gradient model
  ## try cpptsem
  if(tryCpptsem){
    if(tolower(objective)  == "ml"){
      gradientModelcpp <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = dataset))
      if (any(class(gradientModelcpp) == "try-error")){
        message("Using OpenMx. This can slow down the computation considerably.")
        gradientModelcpp <- NULL
      }else{
        gradientModelcpp$computeRAM()
        gradientModelcpp$fitRAM()
        m2LLcpp <- gradientModelcpp$m2LL
        testM2LL <- round(ctsemObject$mxobj$fitfunction$result[[1]] - m2LLcpp,3) == 0
        if (!testM2LL & !forceCpptsem){
          message("Using OpenMx. This can slow down the computation considerably. If you want to use cpptsem nevertheless, set forceCpptsem = TRUE")
          gradientModelcpp <- NULL
        }
      }
    }else if (tolower(objective)  == "kalman"){
      tempCtsemObject <- ctFit(ctmodelobj = ctsemObject, dat = dataset, fit = FALSE, objective = "Kalman")
      tempCtsemObject$mxobj <- mxObject
      gradientModelcpp <- try(regCtsem::cpptsemFromCtsem(ctsemModel = tempCtsemObject, wideData = dataset))
      if (any(class(gradientModelcpp) == "try-error")){
        message("Using OpenMx. This can slow down the computation considerably.")
        gradientModelcpp <- NULL
      }else{
        gradientModelcpp$computeAndFitKalman()
        m2LLcpp <- gradientModelcpp$m2LL
        testM2LL <- round(tempCtsemObject$mxobj$fitfunction$result[[1]] - m2LLcpp,3) == 0
        if (!testM2LL & !forceCpptsem){
          message("Using OpenMx. This can slow down the computation considerably. If you want to use cpptsem nevertheless, set forceCpptsem = TRUE")
          gradientModelcpp <- NULL
        }
      }
    }
  }else{
    message("Using OpenMx. This can slow down the computation considerably.")
    gradientModelcpp <- NULL
  }
  if(is.null(gradientModelcpp)){
    gradientModel <- OpenMx::mxModel(mxObject,
                                     OpenMx::mxComputeSequence(steps=list(OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
                                                                                                        hessian = FALSE))
                                     ))
  }

  # Progress bar
  if(progressBar){
    pbar <- txtProgressBar(min = 0, max = length(lambdas), initial = 0, style = 3)}

  # iterate over lambda values
  numLambdas <- length(lambdas)
  iteration <- 1
  retryOnce <- TRUE # will retry optimizing once with new starting values
  while(iteration <= numLambdas){
    lambda <- lambdas[iteration]

    # update progress bar
    if(progressBar){
      setTxtProgressBar(pbar, iteration)
    }else if(!is.null(parallelProgressBar)){
      writeProgressError <- try(write.csv2(which(lambdas == lambda),parallelProgressBar$parallelTempFiles[[parallelProgressBar$iteration]], row.names = FALSE))
      if(!any(class(writeProgressError) == "try-error")){
        parallelProgressBar$printProgress(parallelProgressBar$parallelTempFiles, parallelProgressBar$maxItSum, parallelProgressBar$cores)
      }
    }

    lambda = ifelse(scaleLambdaWithN, lambda*sampleSize, lambda) # set lambda*samplesize

    # should the results first be approximated?
    startingValues <- tryApproxFirst(startingValues = startingValues, returnAs = "vector",
                                     approxFirst = approxFirst, numStart = numStart, approxOpt = approxOpt, approxMaxIt = approxMaxIt,
                                     lambda = lambda, lambdas = lambdas,
                                     gradientModelcpp = gradientModelcpp,
                                     mxObject = mxObject,
                                     regOn = regOn, regIndicators = regIndicators, adaptiveLassoWeights = adaptiveLassoWeights, objective = objective, sparseParameters =sparseParameters,
                                     extraTries = extraTries)
    # use Gist
    resGIST <- try(regCtsem::GIST(gradientModel = gradientModel, cppmodel = gradientModelcpp, startingValues = startingValues,
                                  objective = objective, lambda = lambda, adaptiveLassoWeights = adaptiveLassoWeights,
                                  regularizedParameters = regIndicators,
                                  eta = eta, sig = sig, initialStepsize = initialStepsize, stepsizeMin = stepsizeMin, stepsizeMax = stepsizeMax,
                                  GISTLinesearchCriterion = GISTLinesearchCriterion, GISTNonMonotoneNBack = GISTNonMonotoneNBack,
                                  maxIter_out = maxIter_out, maxIter_in = maxIter_in,
                                  break_outer = break_outer, differenceApprox = differenceApprox, verbose = verbose))



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
        cRegM2LL <- ifelse(resGIST$convergence, cM2LL +  regCtsem::exact_getLambda(lambda = lambda,
                                                                                   theta = newValues,
                                                                                   regIndicators = regIndicators,
                                                                                   adaptiveLassoWeights = adaptiveLassoWeights), Inf)
      }else{
        newValues <- OpenMx::omxGetParameters(resGIST$model)
        startingValues <- newValues
        # save fit
        cM2LL <- ifelse(resGIST$convergence, resGIST$model$fitfunction$result[[1]], Inf)
        cRegM2LL <- ifelse(resGIST$convergence, cM2LL +  regCtsem::exact_getLambda(lambda = lambda,
                                                                                   theta = newValues,
                                                                                   regIndicators = regIndicators,
                                                                                   adaptiveLassoWeights = adaptiveLassoWeights), Inf)
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
#' @param cppmodel model of type cpptsemmodel
#' @param startingValues named vector with starting values
#' @param objective "ML" for maximum likelihood SEM, "Kalman" for Kalman filter
#' @param lambda penalty value
#' @param adaptiveLassoWeights named vector with adaptive lasso weights
#' @param regularizedParameters named vector of regularized parameters
#' @param eta if the current step size fails, eta will decrease the step size. Must be > 1
#' @param sig GIST: sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
#' @param initialStepsize initial stepsize to be tried in the outer iteration
#' @param stepsizeMin Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param stepsizeMax Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param GISTLinesearchCriterion criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
#' @param GISTNonMonotoneNBack in case of non-monotone line search: Number of preceding regM2LL values to consider
#' @param maxIter_out maximal number of outer iterations
#' @param maxIter_in maximal number of inner iterations
#' @param break_outer stopping criterion for the outer iteration.
#' @param verbose set to 1 to print additional information and plot the convergence
#' @export
GIST <- function(gradientModel, cppmodel, startingValues, objective, lambda, adaptiveLassoWeights, regularizedParameters,
                 eta = 1.5, sig = .2, initialStepsize = 1, stepsizeMin = 0, stepsizeMax = 999999999,
                 GISTLinesearchCriterion = "monotone", GISTNonMonotoneNBack = 5,
                 maxIter_out = 100, maxIter_in = 100,
                 break_outer, differenceApprox, verbose = 0, silent = FALSE){
  break_crit <- names(break_outer)
  # iteration counter
  k_out <- 1
  convergence <- TRUE

  parameterNames <- names(startingValues)
  adaptiveLassoWeightsMatrix <- diag(adaptiveLassoWeights[parameterNames])
  rownames(adaptiveLassoWeightsMatrix) <- parameterNames
  colnames(adaptiveLassoWeightsMatrix) <- parameterNames

  # set parameters
  if(!is.null(cppmodel)){
    cppmodel$setParameterValues(startingValues, names(startingValues))
    if(tolower(objective) == "ml"){
      invisible(capture.output(out1 <- try(cppmodel$computeRAM(), silent = T), type = "message"))
      invisible(capture.output(out2 <- try(cppmodel$fitRAM(), silent = T), type = "message"))
      out3 <- exact_getCppGradients(cppmodel = cppmodel, objective = objective)
    }else{
      invisible(capture.output(out1 <- try(cppmodel$computeAndFitKalman(), silent = TRUE), type = "message"))
      out2 <- NA
      out3 <- exact_getCppGradients(cppmodel = cppmodel, objective = objective)
    }
    if(any(class(out1) == "try-error")  |
       any(class(out2) == "try-error") |
       any(class(out3) == "try-error") |
       anyNA(out3)){
      stop("Infeasible starting values in GIST.")
    }

    parameters_km1 <- NULL
    parameters_k <- cppmodel$getParameterValues()
    parameterNames <- names(parameters_k)
    gradients_km1 <- NULL
    gradients_k <- out3
    m2LL_k <- cppmodel$m2LL
  }else{
    gradientModel <- OpenMx::omxSetParameters(model = gradientModel, labels = names(startingValues), values = startingValues)
    gradientModel <- try(mxRun(gradientModel, silent = TRUE))
    if(any(class(gradientModel) == "try-error")){
      stop("Infeasible starting values in GIST.")
    }
    # extract gradients:
    gradients_k <- gradientModel$compute$steps[[1]]$output[["gradient"]]
    gradients_k <- gradients_k[,differenceApprox] # use specified gradient approximation

    parameters_km1 <- NULL
    parameters_k <- OpenMx::omxGetParameters(gradientModel)
    parameterNames <- names(parameters_k)
    gradients_km1 <- NULL
    m2LL_k <- gradientModel$fitfunction$result[[1]]
  }

  regM2LL_k <- m2LL_k + regCtsem::exact_getLambda(lambda = lambda,
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

      if(!is.null(cppmodel)){
        # update parameters
        cppmodel$setParameterValues(parameters_kp1, names(parameters_kp1))
        # compute likelihood
        if(tolower(objective) == "ml"){
          invisible(capture.output(out1 <- try(cppmodel$computeRAM(), silent = T), type = "message"))
          invisible(capture.output(out2 <- try(cppmodel$fitRAM(), silent = T), type = "message"))
        }else{
          invisible(capture.output(out1 <- try(cppmodel$computeAndFitKalman(), silent = TRUE), type = "message"))
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

        m2LL_kp1 <- cppmodel$m2LL

        if(is.na(m2LL_kp1) | is.infinite(m2LL_kp1)){

          # update step size
          stepsize <- eta*stepsize

          # update iteration counter
          k <- k+1

          # skip rest
          next
        }
      }else{

        # update parameters
        gradientModel <- OpenMx::omxSetParameters(model = gradientModel, labels = names(parameters_kp1), values = parameters_kp1)
        # compute likelihood
        gradientModel <- try(mxRun(gradientModel, silent = TRUE))

        if(any(class(gradientModel) == "try-error")){

          # update step size
          stepsize <- eta*stepsize

          # update iteration counter
          k <- k+1

          # skip rest
          next
        }

        # extract fit
        m2LL_kp1 <- gradientModel$fitfunction$result[[1]]

        if(is.na(m2LL_kp1) | is.infinite(m2LL_kp1)){

          # update step size
          stepsize <- eta*stepsize

          # update iteration counter
          k <- k+1

          # skip rest
          next
        }

      }
      regM2LL_kp1 <- m2LL_kp1 + regCtsem::exact_getLambda(lambda = lambda,
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

    if(!is.null(cppmodel)){
      if(tolower(objective) == "ml"){
        out3 <- exact_getCppGradients(cppmodel = cppmodel, objective = objective)
      }else{
        out3 <- exact_getCppGradients(cppmodel = cppmodel, objective = objective)
      }
      if(any(class(out3) == "try-error") |
         any(is.na(out3))){
        if(!silent){
          convergence <- FALSE
          stop("No gradients in GIST")}
      }

      gradients_kp1 <- out3
    }else{

      # extract gradients:
      gradients_kp1 <- gradientModel$compute$steps[[1]]$output[["gradient"]]
      gradientLabels <- rownames(gradients_kp1)
      gradients_kp1 <- gradients_kp1[,differenceApprox] # use specified gradient approximation
      names(gradients_kp1) <- gradientLabels
    }

    # update parameters for next iteration
    parameters_km1 <- parameters_k
    parameters_k <- parameters_kp1
    gradients_km1 <- gradients_k
    gradients_k <- gradients_kp1
    m2LL_k <- m2LL_kp1
    regM2LL_k <- regM2LL_kp1

    # break outer loop if stopping criterion is satisfied
    k_out <- k_out + 1
    regM2LL[k_out] <- regM2LL_k
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

    if(break_crit == "gradient"){

      # Gradient based break condition: Problem: gradients extremely dependent on the epsilon in the approximation
      parameterMatrix <- matrix(parameters_k[parameterNames], ncol = 1)
      rownames(parameterMatrix) <- parameterNames
      gradientMatrix <- matrix(gradients_k[parameterNames], ncol = 1)
      rownames(gradientMatrix) <- parameterNames
      subgradients <- exact_getSubgradients(theta = parameterMatrix, jacobian = gradientMatrix,
                                            regIndicators = regularizedParameters, lambda = lambda,
                                            lineSearch = NULL, adaptiveLassoWeightsMatrix = adaptiveLassoWeightsMatrix)

      breakOuter <- max(abs(subgradients)) < break_outer
    }else if(break_crit == "parameterChange"){

      breakOuter <- (sqrt(sum((parameters_k - parameters_km1)^2))/sqrt(sum((parameters_km1)^2))) < break_outer

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

  if(!is.null(cppmodel)){
    return(list("type" = "cpptsem",
                "model" = cppmodel,
                "regM2LL" = regM2LL,
                "convergence" = convergence))
  }
  return(list("type" = "mx",
              "model" = gradientModel,
              "regM2LL" = regM2LL,
              "convergence" = convergence))
}




