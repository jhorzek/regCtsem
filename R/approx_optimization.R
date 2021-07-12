#### The following functions implement approximate optimization in OpenMx


#### Main functions ####
#' approx_initializeModel
#'
#' sets up an approximate regularized ctsem
#'
#' NOTE: Function located in file approx_optimization.R
#'
#' @param ctsemObject Fitted object of class ctsemFit
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param sampleSize sample size
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso for optimization = exact
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param cvSample cross-validation sample. Has to be of type mxData
#' @param autoCV Boolean: Should automatic cross-validation be used?
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
#' @param zeroThresh threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
#' @param cores how many computer cores should be used?
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @param silent silent execution
#' @param progressBar Boolean: Should a progress bar be displayed
#' @param parallelProgressBar list: used internally to display progress when executing in parallel
#'
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
approx_initializeModel <- function(  # model
  ctsemObject = NULL,
  mxObject = NULL,
  sampleSize = NULL,
  # penalty settings
  regOn = "DRIFT",
  regIndicators,
  lambda,
  penalty = "lasso",
  adaptiveLassoWeights = NULL,
  # fit settings
  returnFitIndices = TRUE,
  cvSample = NULL,
  autoCV = FALSE,
  k = 5,
  # optimization settings
  epsilon = .001,
  zeroThresh = .001,
  # additional settings
  cores = 1,
  verbose = 0,
  silent = FALSE,
  progressBar = TRUE,
  parallelProgressBar = NULL){


  mxNumObs <- OpenMx::mxMatrix(type= "Full", nrow= 1, ncol = 1, free = FALSE, values = sampleSize,name = "numObs") # define numObs as mxMatrix

  # Define provided mxObject as submodel
  submodel <- OpenMx::mxModel(mxObject, name = "submodel") # has all the parameters and the base fit function (FML or FIML)

  # create fitFunString
  fitFunString <- approx_createFitFunString(penalty = penalty, epsilon = epsilon, silent = silent)

  # define fitfunction:
  regFitAlgebra <- OpenMx::mxAlgebraFromString(fitFunString, name = "regFitAlgebra")
  regFitFunction <- OpenMx::mxFitFunctionAlgebra("regFitAlgebra")

  # Add lambda and regIndicator
  selectedDRIFTValues <- OpenMx::mxMatrix(type = "Full", values = regIndicators, free = F, name = "selectedDRIFTValues")
  lambda <- OpenMx::mxMatrix(type = "Full", values = lambda, nrow = 1, ncol = 1, free = F, name = "lambda")

  # Add Drift weights
  driftWeights <- adaptiveLassoWeights[grepl("drift", names(adaptiveLassoWeights))]
  rowVecDriftWeights <- matrix(driftWeights, nrow = 1)
  rowVecDriftWeights <- OpenMx::mxMatrix(type = "Full", values = rowVecDriftWeights, labels = paste0("stdizer_",names(driftWeights)), free = F, name = "rowVecDriftWeights")

  # Combine everything in a new model:

  internalModel <- OpenMx::mxModel(model= "internalModel", # model that is returned
                                   submodel, # BaseModel is a submodel of internalModel The elements of BaseModel can be accessed by the internalModel
                                   mxNumObs,
                                   selectedDRIFTValues,
                                   rowVecDriftWeights,
                                   lambda,
                                   regFitAlgebra,
                                   regFitFunction)

  internalModel <- OpenMx::mxOption(internalModel, "Calculate Hessian", "No")
  internalModel <- OpenMx::mxOption(internalModel, "Standard Errors", "No")

  return(internalModel)
}

#' approx_iterateOverLambdas
#'
#' loops over lambdas if optimization = "approx"
#'
#' NOTE: Function located in file approx_optimization.R
#'
#' @param ctsemObject Fitted object of class ctsemFit
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param sampleSize sample size
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso for optimization = exact
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
approx_iterateOverLambdas <- function(  # model
  ctsemObject = NULL,
  mxObject = NULL,
  sampleSize = NULL,
  # penalty settings
  regOn = "DRIFT",
  regIndicators,
  lambdas,
  penalty = "lasso",
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

  fitAndParameters <- matrix(NA, nrow = length(parameterLabels)+length(fitLabels), ncol = length(lambdas))
  rownames(fitAndParameters) <- c(fitLabels, parameterLabels)
  colnames(fitAndParameters) <- lambdas

  if(progressBar){
    pbar <- txtProgressBar(min = 0, max = length(lambdas), initial = 0, style = 3)
  }
  for(iteration in 1:length(lambdas)){

    # create model
    regModel <- regCtsem::approx_initializeModel(ctsemObject = ctsemObject, mxObject = mxObject, sampleSize = ifelse(scaleLambdaWithN,sampleSize,1),
                                                 # penalty settings
                                                 regOn = regOn, regIndicators = regIndicators, lambda = lambdas[iteration],
                                                 penalty = penalty, adaptiveLassoWeights = adaptiveLassoWeights,
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
    fitAndParameters[parameterLabels, as.character(lambdas[iteration])] <- OpenMx::cvectorize(OpenMx::omxGetParameters(regModel$submodel))
    fits <- regCtsem::approx_getFitIndices(fittedRegModel = regModel, sampleSize = sampleSize,
                                           # penalty settings
                                           regOn = regOn, regIndicators = regIndicators,
                                           # fit settings
                                           returnFitIndices = returnFitIndices, cvSample = cvSample,
                                           # optimization settings
                                           zeroThresh = zeroThresh)
    # extract fit
    fitAndParameters[fitLabels,as.character(lambdas[iteration])] <- fits[fitLabels,]

    # set progress bar for single core execution
    if(progressBar){setTxtProgressBar(pbar,iteration)}
    # set progress bar for multi core execution
    if(!is.null(parallelProgressBar)){
      writeProgressError <- try(write.csv2(which(lambdas == lambdas[iteration]),parallelProgressBar$parallelTempFiles[[parallelProgressBar$iteration]], row.names = FALSE))
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
                                        lambdas = lambdas)
  }


  return(fitAndParameters)

}

#### Define new fitting function ####

#' approx_createFitFunString
#'
#' Creates the new approx_createFitFunString for specified regularization
#'
#' NOTE: Function located in file approx_optimization.R
#'
#' @param penalty type of penalty. Currently supported are lasso and ridge
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentible one
#' @param silent Boolean indicating if print statements will be suppressed
#' @author Jannik Orzek
#' @keywords internal
#' @export
approx_createFitFunString = function(penalty, epsilon = .00001, silent = FALSE){
  ## unregularized fitfunction:
  fitFunString <- "submodel.fitfunction"

  ## add penalty
  if(tolower(penalty) == "lasso" | tolower(penalty) == "adaptivelasso"){
    addPenalty <- paste("numObs*(lambda*(omxSelectCols(rowVecDriftWeights,t(cvectorize(selectedDRIFTValues)))%*%omxSelectRows(cvectorize((submodel.DRIFT^2 + ",epsilon,")^(1/2)), cvectorize(selectedDRIFTValues))))", sep = "")
  }else if(penalty == "ridge"){
    addPenalty <- "numObs*(lambda*(omxSelectCols(rowVecDriftWeights,t(cvectorize(selectedDRIFTValues)))%*%omxSelectRows(cvectorize((submodel.DRIFT)^2), cvectorize(selectedDRIFTValues))))"

  }else if(penalty == "elastic"){
    stop("not implemented")
  }else{
    stop(paste("could not find pealty function ", penalty, ". Currently supported are lasso, ridge, and adaptiveLasso.", sep = ""))
  }
  fitFunString <- paste(fitFunString, addPenalty, sep = "+")
  return(fitFunString)
}


#### Fit ####

#' approx_getCVFit
#'
#' computes cv fit for approximate optimization
#'
#' NOTE: Function located in file approx_optimization.R
#'
#' @param ctsemObject Fitted object of class ctsemFit
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param cvSample only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param fitAndParameters table with fit results and parameter values
#' @param parameterLabels vector with parameter labels
#' @param lambdas vector with lambdas
#' @author Jannik Orzek
#' @export
approx_getCVFit <- function(mxObject, ctsemObject, cvSample, objective, fitAndParameters, parameterLabels, lambdas){
  if(tolower(objective) == "kalman"){
    # create Model
    suppressMessages(invisible(capture.output(ctsemTemp <- ctFit(dat = cvSample, ctmodelobj = ctsemObject$ctmodelobj, objective = "Kalman", useOptimizer = FALSE))))
    cvModel <- ctsemTemp$mxobj

  }else{
    # create Model
    cvModel <- mxObject
    # replace sample
    cvModel$data <- cvSample
  }
  for(iteration in 1:length(lambdas)){

    currentParameters <- fitAndParameters[parameterLabels,iteration]
    # set parameters to iteration parameters
    cvModel <- try(OpenMx::omxSetParameters(model = cvModel, labels = parameterLabels, values = currentParameters))

    # compute -2LL
    if(any(class(cvModel) == "try-error")){next}
    fitCvModel <- OpenMx::mxRun(cvModel, useOptimizer = F, silent = T)
    fitAndParameters["cvM2LL", iteration] <- fitCvModel$fitfunction$result[[1]]

  }

  return(fitAndParameters)
}


#' approx_getFitIndices
#'
#' computes fit indices
#'
#' NOTE: Function located in file approx_optimization.R
#'
#' @param fittedRegModel fitted regCtsem object
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param cvSample cross-validation sample. Has to be of type mxData
#' @param zeroThresh threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
#' @param sampleSize sample size
#' @author Jannik Orzek
#' @export
approx_getFitIndices <- function(fittedRegModel,
                                 sampleSize,
                                 # penalty settings
                                 regOn = "DRIFT",
                                 regIndicators,
                                 # fit settings
                                 returnFitIndices = TRUE,
                                 cvSample = NULL,
                                 # optimization settings
                                 zeroThresh = .001){
  fitLabels <- c("regM2LL", "m2LL")
  if(returnFitIndices){
    fitLabels <- c(fitLabels, "AIC", "BIC", "estimatedParameters")
  }
  if(!is.null(cvSample)){
    fitLabels <- c(fitLabels, "cvM2LL")
  }

  fitTable <- matrix(NA, nrow = length(fitLabels), ncol = 1)
  colnames(fitTable) <- fittedRegModel$lambda$values
  rownames(fitTable) <- fitLabels

  fitTable["regM2LL",] <- fittedRegModel$fitfunction$result[[1]]
  fitTable["m2LL",] <- fittedRegModel$submodel$fitfunction$result[[1]]

  if(returnFitIndices){
    currentParameters <- OpenMx::omxGetParameters(fittedRegModel)

    driftLabels <- fittedRegModel$submodel[[regOn]]$labels

    regularizedDrifts <- driftLabels[regIndicators == 1]
    setZero <- abs(currentParameters[regularizedDrifts]) <= zeroThresh
    estimatedParameters <- length(currentParameters)-sum(setZero)

    fitTable["AIC", as.character(fittedRegModel$lambda$values)] <- fitTable["m2LL",] + 2*estimatedParameters
    fitTable["BIC", as.character(fittedRegModel$lambda$values)] <- fitTable["m2LL",] + log(sampleSize)*estimatedParameters
    fitTable["estimatedParameters", as.character(fittedRegModel$lambda$values)] <- estimatedParameters
  }

  return(fitTable)

}



