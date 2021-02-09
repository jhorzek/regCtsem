#' approx_regCtsem
#'
#' creates a regCtsem object for approximate optimization
#'
#' @param ctsemObject Fitted object of class ctsemOMX
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param dataset only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param regValues vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param regValuesAutoLength if regValues == "auto", regValuesAutoLength will determine the number of regValues tested.
#' @param penalty Currently supported are lasso, adaptiveLasso, and ridge for optimization = approx and lasso and adaptiveLasso for optimization = exact
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param elastic_alpha placehoder for elastic net. NOT YET IMPLEMENTED
#' @param elastic_gamma placehoder for elastic net. NOT YET IMPLEMENTED
#' @param standardizeDrift Boolean: Should Drift parameters be standardized automatically using T0VAR?
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param cvSample cross-validation sample. Has to be of type mxData
#' @param autoCV Boolean: Should automatic cross-validation be used?
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param KalmanStartValues Optional starting values for the parameters when using Kalman filter
#' @param optimizeKalman Boolen: Should the Kalman model be optimized in OpenMx first? If you want the Kalman model to start optimizing in regCtsem from the provided KalmanStartValues and not use OpenMx to optimize the initial Kalman model, set to FALSE
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
#' @param zeroThresh threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
#' @param extraTries number of extra tries in mxTryHard
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param cores how many computer cores should be used?
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @param silent silent execution
#' @param progressBar Boolean: Should a progress bar be displayed
#' @param parallelProgressBar list: used internally to display progress when executing in parallel
#'
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
approx_regCtsem <- function(  # model
  ctsemObject = NULL,
  mxObject = NULL,
  dataset = NULL,
  # penalty settings
  regOn = "DRIFT",
  regIndicators,
  regValues,
  regValuesAutoLength = 50,
  penalty = "lasso",
  adaptiveLassoWeights = NULL,
  elastic_alpha = NULL,
  elastic_gamma = NULL,
  standardizeDrift = FALSE,
  # fit settings
  returnFitIndices = TRUE,
  cvSample = NULL,
  autoCV = FALSE,
  k = 5,
  # optimization settings
  objective = "ML",
  KalmanStartValues = NULL,
  optimizeKalman = TRUE,
  epsilon = .001,
  zeroThresh = .001,
  extraTries = 1,
  # additional settings
  scaleLambdaWithN = T,
  cores = 1,
  verbose = 0,
  silent = FALSE,
  progressBar = TRUE,
  parallelProgressBar = NULL){

  if(is.null(mxObject)){
    mxObject <- ctsemObject$mxobj
  }
  approxArgsIn <- as.list(environment())

  returnList <- list("setup" = approxArgsIn)



  parameterLabels <- names(OpenMx::omxGetParameters(mxObject))

  # set up automatic cross-validation
  if(autoCV){

    # fit model
    resultsCV <- iterateOverCVFolds(argsIn = approxArgsIn, objective = objective, optimization = "approx")
    fit <- resultsCV$fit

    fit["mean CV fit",] <- apply(fit,2,mean,na.rm =TRUE)
    if(any(is.na(fit))){
      warning("NAs in fit results. Consider rerunning with different starting values.")
    }
    returnList <- rlist::list.append(returnList, "fit" = fit)
    returnList <- c(returnList, "folds" = resultsCV$folds)
    return(returnList)
  }

  ##### no automatic cross-validation ####

  # if Kalman: fit initial model
  if(tolower(objective) == "kalman"){
    if(optimizeKalman){
      mxObject <- OpenMx::mxTryHardctsem(mxObject, silent = TRUE, extraTries = extraTries)
      approxArgsIn$KalmanStartValues <- omxGetParameters(mxObject)
      approxArgsIn$optimizeKalman <- FALSE
    }else{
      mxObject <- OpenMx::mxRun(mxObject, silent = TRUE, useOptimizer = FALSE)
    }
    approxArgsIn$mxObject <- mxObject
    sampleSize <- nrow(dataset)
  }else{
    sampleSize <- mxObject$data$numObs
  }

  # set adaptiveLassoWeights

  adaptiveLassoWeights <- getAdaptiveLassoWeights(mxObject = mxObject, penalty = penalty, adaptiveLassoWeights = adaptiveLassoWeights, standardizeDrift = standardizeDrift)
  returnList$setup$adaptiveLassoWeights <- adaptiveLassoWeights

  # set regValues if regValues == "auto"
  if(any(regValues == "auto")){
    regIndicatorsString <- mxObject[[regOn]]$labels[regIndicators==1]
    maxRegValue <- getMaxRegValue(mxObject = mxObject,
                                  regIndicators = regIndicatorsString,
                                  differenceApprox = "central",
                                  adaptiveLassoWeights = adaptiveLassoWeights)
    regValues <- seq(0,maxRegValue$maxRegValue, length.out = regValuesAutoLength)
    # if scaleLambdaWithN = TRUE, the lambda values will be multiplied with N later on
    # we need to ensure that this will not change the values:
    if(scaleLambdaWithN){
      regValues <- regValues/sampleSize
    }
    returnList$setup$regValues <- regValues
    approxArgsIn$regValues <- regValues
  }

  if(cores == 1){

    # start fitting models
    fitAndParameters <- try(regCtsem::approx_iterateOverRegValues(ctsemObject = ctsemObject, mxObject = mxObject, sampleSize = sampleSize,
                                                                   # penalty settings
                                                                   regOn = regOn, regIndicators = regIndicators, regValues = regValues,
                                                                   penalty = penalty, elastic_alpha = elastic_alpha, elastic_gamma = elastic_gamma,
                                                                   adaptiveLassoWeights = adaptiveLassoWeights,
                                                                   # fit settings
                                                                   returnFitIndices = returnFitIndices, cvSample = cvSample,
                                                                   autoCV = autoCV, k = k,
                                                                   # optimization settings
                                                                   objective = objective, epsilon = epsilon, zeroThresh = zeroThresh, "extraTries" = extraTries,
                                                                   # additional settings
                                                                   scaleLambdaWithN, cores  = cores, verbose = verbose, silent = silent,
                                                                   progressBar = progressBar, parallelProgressBar = parallelProgressBar))

    return(rlist::list.append(returnList, "fitAndParameters" = fitAndParameters))
  }

  # multi core

  results <- try(do.call(regCtsem::approx_parallelRegCtsem, approxArgsIn))
  if(autoCV){
    return(c(returnList, results))
  }
  return(rlist::list.append(returnList, "fitAndParameters" = results))

}



