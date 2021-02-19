#' exact_parallelRegCtsem
#'
#' parallel processing for regCtsem with optimization = "approx"
#'
#' @param ctsemObject Fitted object of class ctsemOMX
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param regValues vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param regValuesAutoLength if regValues == "auto", regValuesAutoLength will determine the number of regValues tested.
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso for optimization = exact
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
#' @keywords internal
#' @import ctsemOMX
#' @importFrom foreach %dopar%
#' @export

# define operators
`%dopar%` <- foreach::`%dopar%`

approx_parallelRegCtsem <- function(
  # model
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
  parallelProgressBar = NULL
){
  parameterLabels <- names(OpenMx::omxGetParameters(mxObject))
  fitLabels <- c("m2LL", "regM2LL", "AIC", "BIC", "cvM2LL", "estimatedParameters")


  # set up multicore execution

  # returns a global cluster variable
  regCtsemClusterVariable <<- snow::makeSOCKcluster(cores, outfile="")
  doSNOW::registerDoSNOW(regCtsemClusterVariable)


  # Split the regValues in equal regValueBins:

  regValueBins <- split(regValues,
                        cut(x = 1:length(regValues),
                            breaks = cores,
                            labels = paste0("core", 1:cores)))

  # for progress printing
  parallelTempFiles <- vector(mode = "list", length = cores)
  names(parallelTempFiles) <- paste0("core", 1:cores)

  for(core in 1:cores){
    parallelTempFiles[core] <- tempfile(pattern = paste0("core", core), fileext = ".txt")
    write.csv2(0,parallelTempFiles[[core]], row.names = FALSE)
  }

  maxItSum <- length(regValues)

  printProgress <- function(parallelTempFiles, maxItSum, cores){
    itSum <- 0
    for(core in 1:cores){
      # read number of iterations of this core
      coreIts <- try(read.csv2(parallelTempFiles[[core]])[[1]])
      if(!any(class(coreIts) == "try-error")){
        # sum number of iterations
        itSum <- itSum + coreIts
      }
    }
    # print number of total iterations
    cat(paste0("Iteration ", itSum, " of ", maxItSum),"\r")
    flush.console()
  }


  fitAndParameters_combined <- foreach::foreach(iteration = 1:length(regValueBins), .combine = cbind,
                                                .multicombine = TRUE,.packages = c("OpenMx","regCtsem"),
                                                .inorder = TRUE,
                                                .errorhandling = "remove",
                                                .verbose = FALSE
  ) %dopar% {
    # regValues in this core
    regValueBin <- regValueBins[[iteration]]

    # set up model
    iterationParallelProgressBar = list("printProgress" = printProgress,
                                        "parallelTempFiles" = parallelTempFiles,
                                        "maxItSum" = maxItSum,
                                        "cores" = cores,
                                        "iteration" = iteration)

    iterationResult <- approx_regCtsem(ctsemObject = ctsemObject, mxObject = mxObject, dataset = dataset,
                                       # penalty settings
                                       regOn = regOn, regIndicators = regIndicators, regValues = regValueBin,
                                       penalty = penalty, elastic_alpha = elastic_alpha, elastic_gamma = elastic_gamma,
                                       standardizeDrift = standardizeDrift,
                                       adaptiveLassoWeights = adaptiveLassoWeights,
                                       # fit settings
                                       returnFitIndices = returnFitIndices, cvSample = cvSample,
                                       autoCV = autoCV, k = k, extraTries = extraTries,
                                       # optimization settings
                                       objective = objective, epsilon = epsilon, zeroThresh = zeroThresh,
                                       # additional settings
                                       scaleLambdaWithN = scaleLambdaWithN, cores  = 1, verbose = verbose, silent = silent,
                                       progressBar = FALSE, parallelProgressBar = iterationParallelProgressBar)

    return(iterationResult$fitAndParameters)
  }

  # Stop cluster
  snow::stopCluster(regCtsemClusterVariable)
  rm(regCtsemClusterVariable, pos = ".GlobalEnv")

  return(fitAndParameters_combined)
}




