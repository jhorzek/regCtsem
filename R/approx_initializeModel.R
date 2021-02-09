#' approx_initializeModel
#'
#' sets up an approximate regularized ctsem
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
  regValue,
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
  fitFunString <- approx_createFitFunString(penalty = penalty, elastic_alpha = elastic_alpha, elastic_gamma = elastic_gamma, epsilon = epsilon, silent = silent)

  # define fitfunction:
  regFitAlgebra <- OpenMx::mxAlgebraFromString(fitFunString, name = "regFitAlgebra")
  regFitFunction <- OpenMx::mxFitFunctionAlgebra("regFitAlgebra")

  # Add regValue and regIndicator
  selectedDRIFTValues <- OpenMx::mxMatrix(type = "Full", values = regIndicators, free = F, name = "selectedDRIFTValues")
  regValue <- OpenMx::mxMatrix(type = "Full", values = regValue, nrow = 1, ncol = 1, free = F, name = "regValue")

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
                                   regValue,
                                   regFitAlgebra,
                                   regFitFunction)

  internalModel <- OpenMx::mxOption(internalModel, "Calculate Hessian", "No")
  internalModel <- OpenMx::mxOption(internalModel, "Standard Errors", "No")

  return(internalModel)
}



