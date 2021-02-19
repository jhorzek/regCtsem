#' approx_iterateOverRegValues
#'
#' loops over regValues if optimization = "approx"
#'
#' @param ctsemObject Fitted object of class ctsemOMX
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param sampleSize sample size
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param regValues vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso for optimization = exact
#' @param elastic_alpha placehoder for elastic net. NOT YET IMPLEMENTED
#' @param elastic_gamma placehoder for elastic net. NOT YET IMPLEMENTED
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param cvSample cross-validation sample. Has to be of type mxData
#' @param autoCV Boolean: Should automatic cross-validation be used?
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
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
approx_iterateOverRegValues <- function(  # model
  ctsemObject = NULL,
  mxObject = NULL,
  sampleSize = NULL,
  # penalty settings
  regOn = "DRIFT",
  regIndicators,
  regValues,
  penalty = "lasso",
  elastic_alpha = NULL,
  elastic_gamma = NULL,
  adaptiveLassoWeights = NULL,
  # fit settings
  returnFitIndices = TRUE,
  cvSample = NULL,
  autoCV = FALSE,
  k = 5,
  # optimization settings
  objective = "ML",
  epsilon = .001,
  zeroThresh = .001,
  extraTries = 1,
  # additional settings
  scaleLambdaWithN,
  cores = 1,
  verbose = 0,
  silent = FALSE,
  progressBar = TRUE,
  parallelProgressBar = NULL){

  parameterLabels <- names(OpenMx::omxGetParameters(mxObject))
  fitLabels <- c("regM2LL", "m2LL")
  if(returnFitIndices){
    fitLabels <- c(fitLabels, "AIC", "BIC", "estimatedParameters")
  }
  if(!is.null(cvSample)){
    fitLabels <- c(fitLabels, "cvM2LL")
  }

  fitAndParameters <- matrix(NA, nrow = length(parameterLabels)+length(fitLabels), ncol = length(regValues))
  rownames(fitAndParameters) <- c(fitLabels, parameterLabels)
  colnames(fitAndParameters) <- regValues

  if(progressBar){
    pbar <- txtProgressBar(min = 0, max = length(regValues), initial = 0, style = 3)
  }
  for(iteration in 1:length(regValues)){

    # create model
    regModel <- regCtsem::approx_initializeModel(ctsemObject = ctsemObject, mxObject = mxObject, sampleSize = ifelse(scaleLambdaWithN,sampleSize,1),
                                                  # penalty settings
                                                  regOn = regOn, regIndicators = regIndicators, regValue = regValues[iteration],
                                                  penalty = penalty, elastic_alpha = elastic_alpha, elastic_gamma = elastic_gamma,
                                                  adaptiveLassoWeights = adaptiveLassoWeights,
                                                  # fit settings
                                                  returnFitIndices = returnFitIndices, cvSample = cvSample,
                                                  autoCV = autoCV, k = k,
                                                  # optimization settings
                                                  epsilon = epsilon, zeroThresh = zeroThresh,
                                                  # additional settings
                                                  cores  = cores, verbose = verbose, silent = silent,
                                                  progressBar = progressBar, parallelProgressBar = parallelProgressBar)
    # run model
    if(extraTries == 1){
      regModel <- try(expr = OpenMx::mxRun(regModel, silent = TRUE), silent = silent)
    }else{
      suppressMessages(invisible(capture.output(regModel <- try(expr = OpenMx::mxTryHardctsem(regModel, extraTries = extraTries), silent = TRUE))))
    }

    # check if model returned errors:
    if(any(class(regModel) == "try-error")){next} # go to next iteration if errors were encountered

    ## if no errors were produced
    # update model
    mxObject <- regModel$submodel

    # extract parameter estimates
    fitAndParameters[parameterLabels, as.character(regValues[iteration])] <- OpenMx::cvectorize(OpenMx::omxGetParameters(regModel$submodel))
    fits <- regCtsem::approx_getFitIndices(fittedRegModel = regModel, ctsemObject = ctsemObject, sampleSize = sampleSize,
                                            # penalty settings
                                            regOn = regOn, regIndicators = regIndicators,
                                            # fit settings
                                            returnFitIndices = returnFitIndices, cvSample = cvSample,
                                            # optimization settings
                                            zeroThresh = zeroThresh, objective = objective)
    # extract fit
    fitAndParameters[fitLabels,as.character(regValues[iteration])] <- fits[fitLabels,]

    # set progress bar for single core execution
    if(progressBar){setTxtProgressBar(pbar,iteration)}
    # set progress bar for multi core execution
    if(!is.null(parallelProgressBar)){
      writeProgressError <- try(write.csv2(which(regValues == regValues[iteration]),parallelProgressBar$parallelTempFiles[[parallelProgressBar$iteration]], row.names = FALSE))
      if(!any(class(writeProgressError) == "try-error")){
        parallelProgressBar$printProgress(parallelProgressBar$parallelTempFiles, parallelProgressBar$maxItSum, parallelProgressBar$cores)
      }
    }
  }

  # cvFit

  if(!is.null(cvSample)){
    fitAndParameters <- approx_getCVFit(mxObject = mxObject, ctsemObject = ctsemObject,
                                        cvSample = cvSample, objective = objective,
                                        fitAndParameters = fitAndParameters, parameterLabels = parameterLabels,
                                        regValues = regValues)
  }


  return(fitAndParameters)

}




