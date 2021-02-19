#' exact_bfgsGLMNET
#'
#' Performs GLMNET (see Friedman, 2010 & Yuan, 2011) with bfgs approximated Hessian
#'
#' @examples
#' library(regCtsem)
#' library(ctsemOMX)
#'
#' # The following example is taken directly from the examples provided in the ctFit documentation of ctsemOMX
#' ### Example from Voelkle, Oud, Davidov, and Schmidt (2012) - anomia and authoritarianism.
#' data(AnomAuth)
#' AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
#'                          Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = NULL)
#' AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel, fit = T)
#'
#' # with GLMNET
#' ## with standardization
#' reg <- exact_bfgsGLMNET(mxObject = AnomAuthfit$mxobj, regIndicators = c("drift_eta2_eta1", "drift_eta1_eta2"), regValues = rev(seq(0,1,.1)), standardizeDrift = TRUE)
#' reg$regM2LL
#' reg$thetas
#'
#' ## without standardization
#' reg2 <- exact_bfgsGLMNET(mxObject = AnomAuthfit$mxobj, regIndicators = c("drift_eta2_eta1", "drift_eta1_eta2"), regValues = seq(0,1,.1), standardizeDrift = FALSE)
#' reg2$regM2LL
#' reg2$thetas
#' @param ctsemObject if objective = "ML": Fitted object of class ctsem. If you want to use objective = "Kalman", pass an object of type ctsemInit from ctModel
#' @param mxObject Object of type MxModel
#' @param dataset only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators Vector with names of regularized parameters
#' @param regValues Vector with lambda values that should be tried
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param stepSize Initial stepSize of the outer iteration (theta_{k+1} = theta_k + stepSize \* Stepdirection)
#' @param tryCpptsem should regCtsem try to translate the model to cpptsem? This can speed up the computation considerably but might fail for some models
#' @param forceCpptsem should cpptsem be enforced even if results differ from ctsem? Sometimes differences between cpptsem and ctsem can result from problems with numerical precision which will lead to the m,atrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "cpptsem") for more details
#' @param lineSearch String indicating if Wolfe conditions (lineSearch = "Wolfe") should be used in the outer iteration. Set lineSearch = "dynamic" to only use line search if optimization without line search fails.
#' @param c1 c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
#' @param c2 c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
#' @param sig GLMNET & GIST: GLMNET: only relevant when lineSearch = 'GLMNET' | GIST: sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
#' @param gam GLMNET when lineSearch = 'GLMNET'. Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param differenceApprox Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#' @param initialHessianApproximation Which initial hessian approximation should be used? Possible are: 'ident' for an identity matrix and OpenMx (here the hessian approxmiation from the mxObject is used)
#' @param maxIter_out Maximal number of outer iterations
#' @param maxIter_in Maximal number of inner iterations
#' @param maxIter_line Maximal number of iterations for the lineSearch procedure
#' @param eps_out Stopping criterion for outer iterations
#' @param eps_in Stopping criterion for inner iterations
#' @param eps_WW Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param sampleSize sample size for scaling lambda with N
#' @param exactApproximateFirst Should approximate optimization be used first to obtain start values for exact optimization? 1 = only for first regValue, 2 = for all regValues
#' @param exactApproximateFirst3NumStartingValues Used if exactApproximateFirst = 3. regCtsem will try exactApproximateFirst3NumStartingValues+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param exactApproximateFirst3Optimize Used if exactApproximateFirst = 3. Should each of the generated starting values be optimized slightly? This can substantially improve the fit of the generated starting values. 1 = optimization with optim, 2 = optimization with Rsolnp
#' @param exactApproximateFirstMaxIter_out Used if exactApproximateFirst = 3 and exactApproximateFirst3Optimize > 1. How many outer iterations should be given to each starting values vector? More will improve the selected starting values but slow down the computation. If exactApproximateFirst =  4, or exactApproximateFirst = 5 this will control the number of outer iteration in optim or solnp .
#' @param extraTries number of extra tries in mxTryHard for warm start
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @export
exact_bfgsGLMNET <- function(ctsemObject, mxObject, dataset, objective, regOn = "DRIFT", regIndicators, regValues, adaptiveLassoWeights,
                             # additional settings
                             sparseParameters = NULL,
                             tryCpptsem, forceCpptsem = FALSE, stepSize = 1, lineSearch = "none", c1 = .0001, c2 = .9,
                             sig = .2, gam = 0,
                             differenceApprox = "central", initialHessianApproximation = "OpenMx", maxIter_out = 100, maxIter_in = 1000,
                             maxIter_line = 100, eps_out = .0000000001, eps_in = .0000000001, eps_WW = .0001,
                             scaleLambdaWithN = TRUE, sampleSize, exactApproximateFirst = 0,
                             exactApproximateFirst3NumStartingValues = 10,
                             exactApproximateFirst3Optimize = T,
                             exactApproximateFirstMaxIter_out = 5, extraTries = 3, verbose = 0, progressBar = TRUE, parallelProgressBar = NULL){
  # Setup
  ## hessianModel is a models which returns the gradients and the Hessian:
  computeHessian <- !is.matrix(initialHessianApproximation) & is.null(mxObject$output$hessian)
  hessianModel <- OpenMx::mxModel(mxObject,
                                  OpenMx::mxComputeSequence(steps=list(OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
                                                                                                     hessian = computeHessian))
                                  ))
  # get initial parameter values, gradients, and hessian:
  hessianModel <- OpenMx::mxRun(hessianModel, silent = TRUE)

  # get parameter values
  theta_kp1 <- as.matrix(OpenMx::omxGetParameters(hessianModel))
  initialParameters <- theta_kp1

  # get initial gradients
  g_kp1 <- hessianModel$compute$steps[[1]]$output[["gradient"]]
  thetaNames <- rownames(g_kp1)
  g_kp1 <- g_kp1[,differenceApprox] # use specified gradient approximation
  g_kp1 <- matrix(g_kp1, nrow = length(g_kp1), ncol = 1)
  rownames(g_kp1) <- thetaNames
  initialGradients <- g_kp1

  # get initial Hessian
  if(computeHessian){
    H_kp1 <- hessianModel$output$calculatedHessian
    eigenDecomp <- eigen(H_kp1)
    if(any(eigenDecomp$values < 0)){
      warning("Initial Hessian is not positive definite. Flipping Eigen values to obtain a positive definite initial Hessian.")
      D <- abs(diag(eigenDecomp$values))
      L <- eigenDecomp$vectors
      H_kp1 <- L%*%D%*%solve(L)
    }
  }else if(is.matrix(initialHessianApproximation)){
    H_kp1 <- initialHessianApproximation
  }else{
    H_kp1 <- mxObject$output$hessian
  }


  # stop if Hessian is not positive definite
  if(any(eigen(H_kp1)$values < 0)){
    if(is.matrix(initialHessianApproximation) || !is.null(mxObject$output$hessian)){
      H_kp1 <- "ident"
    }
    warning("Initial Hessian is not positive definite and will be approximated.")
    H_kp1 <- regCtsem::exact_initialHessian(mxObject = mxObject,
                                            approximationType = initialHessianApproximation,
                                            estimatedHessian = H_kp1)
  }
  initialHessian <- H_kp1
  # define return values
  thetas <- matrix(NA,
                   nrow = length(theta_kp1),
                   ncol = length(regValues),
                   dimnames = list(names(OpenMx::omxGetParameters(hessianModel)),
                                   regValues))
  m2LL <- rep(NA, length(regValues))
  regM2LL <- rep(NA, length(regValues))

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
  gradientModel <- OpenMx::mxModel(mxObject,
                                   OpenMx::mxComputeSequence(steps=list(OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
                                                                                                      hessian = FALSE))
                                   ))


  # Progress bar
  if(progressBar){
    pbar <- txtProgressBar(min = 0, max = length(regValues), initial = 0, style = 3)}

  # iterate over lambda values
  numRegValues <- length(regValues)
  iteration <- 1
  retryOnce <- TRUE # will retry optimizin once with new starting values
  while(iteration <= numRegValues){
    lambda <- regValues[iteration]
    # update progress bar
    if(progressBar){
      setTxtProgressBar(pbar, iteration)
    }else if(!is.null(parallelProgressBar)){
      writeProgressError <- try(write.csv2(which(regValues == lambda),parallelProgressBar$parallelTempFiles[[parallelProgressBar$iteration]], row.names = FALSE))
      if(!any(class(writeProgressError) == "try-error")){
        parallelProgressBar$printProgress(parallelProgressBar$parallelTempFiles, parallelProgressBar$maxItSum, parallelProgressBar$cores)
      }
    }

    lambda = ifelse(scaleLambdaWithN, lambda*sampleSize, lambda) # set lambda*samplesize

    # should the results first be approximated?
    if((exactApproximateFirst == 1 & lambda == regValues[1]) |
       exactApproximateFirst == 2 |
       (exactApproximateFirst == 4 & is.null(gradientModelcpp)) |
       (exactApproximateFirst == 5 & is.null(gradientModelcpp))){
      startingValues <- as.vector(theta_kp1)
      names(startingValues) <- rownames(theta_kp1)

      approxModel <- regCtsem::approx_initializeModel(mxObject = mxObject,
                                                      sampleSize = 1, # scaling with N is handled above
                                                      regOn = regOn,
                                                      regIndicators = regIndicatorsFromNameToMatrix(mxObject = mxObject, regOn = regOn, regIndicators = regIndicators),
                                                      regValue = lambda,
                                                      adaptiveLassoWeights = adaptiveLassoWeights,
                                                      penalty = "lasso"
      )
      approxModel <- omxSetParameters(approxModel, labels = names(startingValues), values = startingValues)
      suppressMessages(invisible(capture.output(approxModel <- try(expr = OpenMx::mxTryHardctsem(approxModel, extraTries = extraTries), silent = TRUE))))
      theta_kp1 <- as.matrix(omxGetParameters(approxModel))
    }

    if(exactApproximateFirst == 3){
      startingValues <- as.vector(theta_kp1)
      names(startingValues) <- rownames(theta_kp1)
      temp_startValues <- try(exact_tryStartingValues(gradientModel = gradientModel,
                                                      cppmodel = gradientModelcpp,
                                                      objective = objective,
                                                      currentParameters = startingValues,
                                                      sparseParameters = sparseParameters,
                                                      regIndicators = regIndicators,
                                                      newLambda = lambda,
                                                      adaptiveLassoWeights = adaptiveLassoWeights,
                                                      differenceApprox = differenceApprox,
                                                      numStartingValues = exactApproximateFirst3NumStartingValues, optimize = exactApproximateFirst3Optimize,
                                                      numOuter = exactApproximateFirstMaxIter_out), silent = TRUE)
      if(!any(class(temp_startValues) == "try-error")){
        theta_kp1 <- as.matrix(temp_startValues)
      }
    }
    if(exactApproximateFirst == 4 & !is.null(gradientModelcpp)){
      startingValues <- as.vector(theta_kp1)
      names(startingValues) <- rownames(theta_kp1)
      optimized <- try(approx_cpptsemOptim(cpptsemmodel = gradientModelcpp,
                                           regM2LLCpptsem = ifelse(tolower(objective) == "ml",
                                                                   regCtsem::approx_RAMRegM2LLCpptsem,
                                                                   regCtsem::approx_KalmanRegM2LLCpptsem),
                                           gradCpptsem = regCtsem::approx_gradCpptsem,
                                           startingValues = startingValues,
                                           adaptiveLassoWeights = adaptiveLassoWeights,
                                           N = 1, lambda = lambda,
                                           regIndicators = regIndicators,
                                           epsilon = 10^(-8),
                                           maxit = exactApproximateFirstMaxIter_out,
                                           objective, testGradients = TRUE), silent = TRUE)
      if(!any(class(optimized) == "try-error")){
        theta_kp1 <- as.matrix(optimized$parameters)
      }
    }
    if(exactApproximateFirst == 5 & !is.null(gradientModelcpp)){
      startingValues <- as.vector(theta_kp1)
      names(startingValues) <- rownames(theta_kp1)
      optimized <- try(approx_cpptsemSolnp(cpptsemmodel = gradientModelcpp,
                                           regM2LLCpptsem = ifelse(tolower(objective) == "ml",
                                                                   regCtsem::approx_RAMRegM2LLCpptsem,
                                                                   regCtsem::approx_KalmanRegM2LLCpptsem),
                                           gradCpptsem = regCtsem::approx_gradCpptsem,
                                           startingValues = startingValues,
                                           adaptiveLassoWeights = adaptiveLassoWeights,
                                           N = 1, lambda = lambda,
                                           regIndicators = regIndicators,
                                           epsilon = 10^(-8),
                                           maxit = exactApproximateFirstMaxIter_out,
                                           objective = objective), silent = TRUE)
      if(!any(class(optimized) == "try-error")){
        theta_kp1 <- as.matrix(optimized$parameters)
      }
    }

    # outer loop: optimize parameters
    resGLMNET <- try(regCtsem::exact_outerGLMNET(mxObject = mxObject, objective =objective,
                                                 adaptiveLassoWeights = adaptiveLassoWeights, sampleSize = sampleSize,
                                                 gradientModel = gradientModel, gradientModelcpp = gradientModelcpp, parameterNames = thetaNames,
                                                 initialParameters = theta_kp1, initialGradients = g_kp1, initialHessian = H_kp1,
                                                 lambda = lambda, regIndicators = regIndicators,
                                                 stepSize = stepSize, lineSearch = lineSearch, c1 = c1, c2 = c2,
                                                 differenceApprox = differenceApprox, maxIter_out = maxIter_out,
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
      # run final gradientModel
      if(!is.null(resGLMNET$gradientModelcpp)){
        gradientModelcpp <- resGLMNET$gradientModelcpp
        if(tolower(objective) == "ml"){
          out1 <- try(gradientModelcpp$computeRAM(), silent = TRUE)
          out2 <- try(gradientModelcpp$fitRAM(), silent = TRUE)
        }else{
          out1 <- try(gradientModelcpp$computeAndFitKalman(), silent = TRUE)
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
          cM2LL <- ifelse(resGLMNET$convergence, gradientModelcpp$m2LL, Inf)
          cRegM2LL <- ifelse(resGLMNET$convergence, gradientModelcpp$m2LL +  regCtsem::exact_getRegValue(lambda = lambda,
                                                                                                         theta = resGLMNET$theta_kp1,
                                                                                                         regIndicators = regIndicators,
                                                                                                         adaptiveLassoWeights = adaptiveLassoWeights), Inf)
          # save results
          thetas[,iteration] <- resGLMNET$theta_kp1[rownames(thetas),]
          m2LL[iteration] <- cM2LL
          regM2LL[iteration] <- cRegM2LL
        }
      }else{
        gradientModel <- try(OpenMx::mxRun(resGLMNET$gradientModel, silent = TRUE))
        if(class(gradientModel) == "try-error"){
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

        }else{
          theta_kp1 <- resGLMNET$theta_kp1
          # update gradientModel
          gradientModel <- resGLMNET$gradientModel
          # update Hessian approximation
          H_kp1 <- resGLMNET$H_kp1
          # update gradients
          g_kp1 <- resGLMNET$g_kp1
          # save fit
          cM2LL <- ifelse(resGLMNET$convergence, resGLMNET$gradientModel$fitfunction$result[[1]], Inf)
          cRegM2LL <- ifelse(resGLMNET$convergence, resGLMNET$gradientModel$fitfunction$result[[1]] +  regCtsem::exact_getRegValue(lambda = lambda,
                                                                                                                                   theta = resGLMNET$theta_kp1,
                                                                                                                                   regIndicators = regIndicators,
                                                                                                                                   adaptiveLassoWeights = adaptiveLassoWeights), Inf)
          # save results
          thetas[,iteration] <- resGLMNET$theta_kp1[rownames(thetas),]
          m2LL[iteration] <- cM2LL
          regM2LL[iteration] <- cRegM2LL
        }
      }
      iteration <- iteration + 1
    }

  }

  return(list("regValues" = regValues, "thetas" = thetas, "m2LL" = m2LL,"regM2LL" = regM2LL))

}




