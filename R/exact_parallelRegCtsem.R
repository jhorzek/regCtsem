#' exact_parallelRegCtsem
#'
#' parallel processing for regCtsem with optimization = "exact"
#'
#' @param ctsemObject if objective = "ML": Fitted object of class ctsem. If you want to use objective = "Kalman", pass an object of type ctsemInit from ctModel
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
#' @param tryCpptsem should regCtsem try to translate the model to cpptsem? This can speed up the computation considerably but might fail for some models
#' @param forceCpptsem should cpptsem be enforced even if results differ from ctsem? Sometimes differences between cpptsem and ctsem can result from problems with numerical precision which will lead to the m,atrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "cpptsem") for more details
#' @param stepSize GLMNET & GIST: initial step size of the outer iteration
#' @param lineSearch GLMNET: String indicating if Wolfe conditions (lineSearch = "Wolfe") should be used in the outer iteration
#' @param c1 GLMNET: c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
#' @param c2 GLMNET: c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
#' @param sig GLMNET & GIST: GLMNET: only relevant when lineSearch = 'GLMNET' | GIST: sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
#' @param gam GLMNET when lineSearch = 'GLMNET'. Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param differenceApprox GLMNET: Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#' @param initialHessianApproximation GLMNET: Which initial hessian approximation should be used? Possible are: 'ident' for an identity matrix and OpenMx (here the hessian approxmiation from the mxObject is used)
#' @param maxIter_out GLMNET & GIST: Maximal number of outer iterations
#' @param maxIter_in GLMNET & GIST: Maximal number of inner iterations
#' @param maxIter_line GLMNET: Maximal number of iterations for the lineSearch procedure
#' @param eps_out GLMNET: Stopping criterion for outer iterations
#' @param eps_in GLMNET: Stopping criterion for inner iterations
#' @param eps_WW GLMNET: Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned
#' @param eta GIST: if the current step size fails, eta will decrease the step size. Must be > 1
#' @param stepsizeMin GIST: Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param stepsizeMax GIST: Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param GISTLinesearchCriterion criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
#' @param GISTNonMonotoneNBack in case of non-monotone line search: Number of preceding regM2LL values to consider
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param exactApproximateFirst Should approximate optimization be used first to obtain start values for exact optimization? 1 = only for first regValue, 2 = for all regValues
#' @param exactApproximateFirst3NumStartingValues Used if exactApproximateFirst = 3. regCtsem will try exactApproximateFirst3NumStartingValues+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param exactApproximateFirst3Optimize Used if exactApproximateFirst = 3. Should each of the generated starting values be optimized slightly? This can substantially improve the fit of the generated starting values. 1 = optimization with optim, 2 = optimization with Rsolnp
#' @param exactApproximateFirstMaxIter_out Used if exactApproximateFirst = 3 and exactApproximateFirst3Optimize > 1. How many outer iterations should be given to each starting values vector? More will improve the selected starting values but slow down the computation. If exactApproximateFirst =  4, or exactApproximateFirst = 5 this will control the number of outer iteration in optim or solnp .
#' @param extraTries number of extra tries in mxTryHard for warm start
#' @param cores how many computer cores should be used?
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @param silent silent execution
#' @param progressBar Boolean: Should a progress bar be displayed
#' @param parallelProgressBar list: used internally to display progress when executing in parallel
#'
#' @author Jannik Orzek
#' @import ctsemOMX doSNOW rlist
#' @importFrom foreach %dopar%
#' @export

# define operators
`%dopar%` <- foreach::`%dopar%`

exact_parallelRegCtsem <- function(# model
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
  optimizer = "GLMNET",
  objective = "ML",
  KalmanStartValues = NULL,
  optimizeKalman = TRUE,
  sparseParameters = NULL,
  tryCpptsem = TRUE,
  forceCpptsem = FALSE,
  # settings for GLMNET
  stepSize = 1,
  lineSearch = "none",
  c1 = .0001,
  c2 = .9,
  sig = .5,
  gam = 0,
  differenceApprox = "central",
  initialHessianApproximation = "OpenMx",
  maxIter_out = 100,
  maxIter_in = 1000,
  maxIter_line = 100,
  eps_out = .0000000001,
  eps_in = .0000000001,
  eps_WW = .0001,
  # settings for GIST
  eta = 1.5,
  stepsizeMin = 0,
  stepsizeMax = 999999999,
  GISTLinesearchCriterion = "monotone",
  GISTNonMonotoneNBack = 5,
  break_outer = .0001,
  # general
  scaleLambdaWithN = TRUE,
  exactApproximateFirst = 3,
  exactApproximateFirst3NumStartingValues = 10,
  exactApproximateFirst3Optimize = T,
  exactApproximateFirstMaxIter_out = 5,
  extraTries = 3,
  # additional settings
  cores = 1,
  verbose = 0,
  silent = FALSE,
  progressBar = TRUE,
  parallelProgressBar = NULL){

  # catch function call
  fun_call <- as.list(environment())

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

  # for progess printing
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

    # expected results
    fitTable <- matrix(NA, nrow = length(parameterLabels)+length(fitLabels), ncol = length(regValueBin))
    rownames(fitTable) <- c(parameterLabels, fitLabels)
    colnames(fitTable) <- regValueBin

    initialModel <- mxObject

    iterationExactArgsIn <- fun_call
    iterationExactArgsIn$mxObject <- initialModel
    iterationExactArgsIn$regValues <- regValueBin
    iterationExactArgsIn$cores <- 1
    iterationExactArgsIn$progressBar <- FALSE
    iterationExactArgsIn$parallelProgressBar = list("printProgress" = printProgress,
                                                    "parallelTempFiles" = parallelTempFiles,
                                                    "maxItSum" = maxItSum,
                                                    "cores" = cores,
                                                    "iteration" = iteration)

    regModelIteration <- try(do.call(regCtsem::regCtsem, iterationExactArgsIn))


    if(!any(class(regModelIteration) == "try-error")){
      # extract parameter estimates and fit
      fitAndParameters <- regModelIteration$fitAndParameters
    }

    return(fitAndParameters)
  }

  # Stop cluster
  snow::stopCluster(regCtsemClusterVariable)
  rm(regCtsemClusterVariable, pos = ".GlobalEnv")

  return(list("fitAndParameters" = fitAndParameters_combined))
}


