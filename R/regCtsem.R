#' regCtsem
#'
#' main function: performs regularized ctsem
#'
#' @param ctsemObject if objective = "ML": Fitted object of class ctsem. If you want to use objective = "Kalman", pass an object of type ctsemInit from ctModel
#' @param dataset only required if objective = "Kalman" and ctsemObject is of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param regValues vector of penalty values (tuning parameter). E.g., seq(0,1,.01). Alternatively, regValues can be set to "auto"; however, this is still experimental.
#' @param regValuesAutoLength if regValues == "auto", regValuesAutoLength will determine the number of regValues tested.
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso and adaptiveLasso for optimization = exact
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the unregularized parameter estimates.
#' @param elastic_alpha placehoder for elastic net. NOT YET IMPLEMENTED
#' @param elastic_gamma placehoder for elastic net. NOT YET IMPLEMENTED
#' @param standardizeDrift Boolean: Should Drift parameters be standardized automatically using T0VAR?
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param cvSample cross-validation sample. Has to be of type mxData
#' @param autoCV Boolean: Should automatic cross-validation be used?
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @param optimization which optimization procedure should be used. Possible are  "exact" or "approx".
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param KalmanStartValues Optional starting values for the parameters when using Kalman filter
#' @param optimizeKalman Boolen: Should the Kalman model be optimized in OpenMx first? If you want the Kalman model to start optimizing in regCtsem from the provided KalmanStartValues and not use OpenMx to optimize the initial Kalman model, set to FALSE
#' @param sparseParameters labeled vector with parameter estimates of the most sparse model. Required for exactApproximateFirst = 3
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
#' @param zeroThresh threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
#' @param extraTries number of extra tries in mxTryHard
#' @param tryCpptsem should regCtsem try to translate the model to cpptsem? This can speed up the computation considerably but might fail for some models
#' @param forceCpptsem should cpptsem be enforced even if results differ from ctsem? Sometimes differences between cpptsem and ctsem can result from problems with numerical precision which will lead to the m,atrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "cpptsem") for more details
#' @param stepSize GLMNET & GIST: initial step size of the outer iteration
#' @param lineSearch GLMNET: String indicating which linesearch should be used. Defaults to the one described in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Alternatively (not recommended) Wolfe conditions (lineSearch = "Wolfe") can be used in the outer iteration. Setting to "none" is also not recommended!.
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
#' @param exactApproximateFirst Should approximate optimization be used first to obtain start values for exact optimization? 1 and 2 are using OpenMx with 1 = optimization only for first regValue, 2 = optimization for all regValues. 3 ensures that the fit will not be worse than in the sparse model if regValues = "auto" or sparseParameters are provided. To this end, 10 models between the current parameter estimates and the sparse parameter estimates are tested and the one with the lowest regM2LL is used for starting values. 4 = optimizing using optim or Opemx if cpptsem is not available, 5 = optimizing using optim or Opemx if cpptsem is not available (requires installation of Rsolnp)
#' @param exactApproximateFirst3NumStartingValues Used if exactApproximateFirst = 3. regCtsem will try exactApproximateFirst3NumStartingValues+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param exactApproximateFirst3Optimize Used if exactApproximateFirst = 3. Should each of the generated starting values be optimized slightly? This can substantially improve the fit of the generated starting values. 1 = optimization with optim, 2 = optimization with Rsolnp
#' @param exactApproximateFirstMaxIter_out Used if exactApproximateFirst = 3 and exactApproximateFirst3Optimize > 1. How many outer iterations should be given to each starting values vector? More will improve the selected starting values but slow down the computation. If exactApproximateFirst =  4, or exactApproximateFirst = 5 this will control the number of outer iteration in optim or solnp .
#' @param cores how many computer cores should be used?
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @param silent silent execution
#' @param progressBar Boolean: Should a progress bar be displayed
#' @param parallelProgressBar list: used internally to display progress when executing in parallel
#' @param calledInternally Boolean: used internally for skipping checks
#'
#' @examples
#' set.seed(17046)
#'
#' library(regCtsem)
#'
#' ## define the population model:
#'
#' # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
#' ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)
#'
#' generatingModel<-ctsem::ctModel(Tpoints=10,n.latent=2,n.TDpred=0,
#'                                 n.TIpred=0,n.manifest=2,
#'                                 MANIFESTVAR=diag(0,2),
#'                                 LAMBDA=diag(1,2),
#'                                 DRIFT=ct_drift,
#'                                 DIFFUSION=matrix(c(.5,0,0,.5),2),
#'                                 CINT=matrix(c(0,0),nrow=2),
#'                                 T0MEANS=matrix(0,ncol=1,nrow=2),
#'                                 T0VAR=diag(1,2), type = "omx")
#'
#' # simulate a training data set
#' dat <- ctsem::ctGenerate(generatingModel,n.subjects = 100, wide = TRUE)
#'
#' ## Build the analysis model. Note that drift eta1_eta2 is freely estimated
#' # although it is 0 in the population.
#' myModel <- ctsem::ctModel(Tpoints=10,n.latent=2,n.TDpred=0,
#'                           n.TIpred=0,n.manifest=2,
#'                           LAMBDA=diag(1,2),
#'                           MANIFESTVAR=diag(0,2),
#'                           CINT=matrix(c(0,0),nrow=2),
#'                           DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
#'                           T0MEANS=matrix(0,ncol=1,nrow=2),
#'                           T0VAR="auto", type = "omx")
#'
#' # fit the model using ctsemOMX:
#' fit_myModel <- ctsemOMX::ctFit(dat, myModel)
#'
#' # select DRIFT values:
#' regOn = "DRIFT"
#' regIndicators = matrix(c(0,1,1,0), byrow = T, ncol = 2)
#'
#' # perform regularization
#' regModel <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
#'                                     regOn = regOn,
#'                                     regIndicators = regIndicators,
#'                                     regValues = "auto",
#'                                     regValuesAutoLength = 15,
#'                                     penalty = "lasso",
#'                                     standardizeDrift = TRUE,
#'                                     autoCV = F,
#'                                     optimization = "exact",
#'                                     returnFitIndices = TRUE,
#'                                     cores = 1))
#'
#' summary(regModel, criterion = "BIC")
#' plot(regModel, what = "drift")
#' plot(regModel, what = "fit", criterion = c("AIC", "BIC", "m2LL"))
#'
#' # The same regularization can be performed with the approximate optimization:
#' # Note that we are using extraTries to get better parameter estimates
#' regModelApprox <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
#'                                           regOn = regOn,
#'                                           regIndicators = regIndicators,
#'                                           regValues = "auto",
#'                                           regValuesAutoLength = 15,
#'                                           penalty = "lasso",
#'                                           standardizeDrift = T,
#'                                           autoCV = F,
#'                                           optimization = "approx",
#'                                           returnFitIndices = TRUE,
#'                                           cores = 1,
#'                                           epsilon = .001,
#'                                           zeroThresh = .04, # value below which parameters are evaluated as zero
#'                                           extraTries = 5))
#'
#' # Comparison of parameter estimates:
#' round(regModel$fitAndParameters - regModelApprox$fitAndParameters,4)
#'
#' #### regCtsem with Kalman Filter and N>1
#'
#' set.seed(175446)
#'
#' ## define the population model:
#'
#' # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
#' ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)
#'
#' generatingModel<-ctsem::ctModel(Tpoints=100,n.latent=2,
#'                                 n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                                 MANIFESTVAR=diag(0,2),
#'                                 LAMBDA=diag(1,2),
#'                                 DRIFT=ct_drift,
#'                                 DIFFUSION=matrix(c(.5,0,0,.5),2),
#'                                 CINT=matrix(c(0,0),nrow=2),
#'                                 T0MEANS=matrix(0,ncol=1,nrow=2),
#'                                 T0VAR=diag(1,2), type = "omx")
#'
#' # simulate a training data and testing data set
#' traindata <- ctsem::ctGenerate(generatingModel,n.subjects = 10, wide = TRUE)
#' testdata <- ctsem::ctGenerate(generatingModel,n.subjects = 10, wide = TRUE)
#'
#' ## Build the analysis model. Note that drift eta1_eta2 is freely estimated
#' # although it is 0 in the population.
#' myModel <- ctsem::ctModel(Tpoints=100,n.latent=2,n.TDpred=0,
#'                           n.TIpred=0,n.manifest=2,
#'                           LAMBDA=diag(1,2),
#'                           MANIFESTVAR=diag(0,2),
#'                           CINT=matrix(c(0,0),nrow=2),
#'                           DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
#'                           T0MEANS=matrix(0,ncol=1,nrow=2),
#'                           T0VAR="auto", type = "omx")
#'
#' # select DRIFT values:
#' regOn = "DRIFT"
#' regIndicators = matrix(c(0,1,1,0), byrow = T, ncol = 2)
#'
#' regModel <- try(regCtsem::regCtsem(ctsemObject = myModel,
#'                                     dataset = traindata,
#'                                     regOn = regOn,
#'                                     regIndicators = regIndicators,
#'                                     regValues = "auto",
#'                                     regValuesAutoLength = 15,
#'                                     penalty = "lasso",
#'                                     standardizeDrift = TRUE,
#'                                     autoCV = F,
#'                                     cvSample = testdata,
#'                                     optimization = "exact",
#'                                     returnFitIndices = TRUE,
#'                                     cores = 2,
#'                                     objective = "Kalman"))
#'
#' summary(regModel, criterion = "cvM2LL")
#' plot(regModel, what = "fit", criterion = "cvM2LL")
#'
#' @author Jannik Orzek
#' @import ctsemOMX doSNOW rlist
#' @export
regCtsem <- function(
  # model
  ctsemObject = NULL,  mxObject = NULL, dataset = NULL,
  # penalty settings
  regOn = "DRIFT",
  regIndicators,
  regValues = "auto",
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
  optimizer = "GIST",
  objective = "ML",
  KalmanStartValues = NULL,
  optimizeKalman = TRUE,
  sparseParameters = NULL,
  optimization = "exact",
  scaleLambdaWithN = TRUE,
  # for approximate estimation:
  epsilon = .001,
  zeroThresh = .04,
  extraTries = 1,
  # for exact optimization:
  tryCpptsem = TRUE,
  forceCpptsem = FALSE,
  stepSize = 1,
  lineSearch = "GLMNET",
  c1 = .0001,
  c2 = .9,
  sig = 10^(-5),
  gam = 0,
  differenceApprox = "central",
  initialHessianApproximation = "OpenMx",
  maxIter_out = 100,  maxIter_in = 1000,  maxIter_line = 100,
  eps_out = .0000000001,  eps_in = .0000000001,  eps_WW = .0001,
  # settings for GIST
  eta = 2,
  stepsizeMin = 1/(10^30),
  stepsizeMax = 10^30,
  GISTLinesearchCriterion = "monotone",
  GISTNonMonotoneNBack = 5,
  break_outer = 10^(-5),
  exactApproximateFirst = "auto",
  exactApproximateFirst3NumStartingValues = 0,
  exactApproximateFirst3Optimize = 1,
  exactApproximateFirstMaxIter_out = "auto",
  # additional settings
  cores = 1, verbose = 0,  silent = FALSE,  progressBar = TRUE,  parallelProgressBar = NULL, calledInternally = FALSE
){
  # save input in list
  argsIn <- as.list(environment())

  # check setup
  checkSetup(argsIn)
  if(exactApproximateFirst == "auto"){
    if((any(regValues == "auto")) | !is.null(sparseParameters)){
      argsIn$exactApproximateFirst <- 3
    }else{
      argsIn$exactApproximateFirst <- 4
    }
  }

  if(exactApproximateFirst == 3){
    if((!any(regValues == "auto")) & any(is.null(sparseParameters))){
      warning("exactApproximateFirst = 3 requested. This requires either regValues = 'auto' or the sparseParameters being passed to regCtsem. Setting exactApproximateFirst = 4")
      argsIn$exactApproximateFirst <- 4
    }
  }

  if(exactApproximateFirstMaxIter_out == "auto"){
    if(exactApproximateFirst == 3){
      argsIn$exactApproximateFirstMaxIter_out <- 10
    }else{
      argsIn$exactApproximateFirstMaxIter_out <- 200
    }
  }

  if(is.null(mxObject)){
    mxObject <- ctsemObject$mxobj
    argsIn$mxObject <- mxObject
  }

  if(optimization == "exact" & is.numeric(regIndicators)){
    if(tolower(objective) == "kalman"){
      regIndicators <- ctsemObject[[regOn]][regIndicators == 1]
      argsIn$regIndicators <- regIndicators
    }else{
      regIndicators <- mxObject[[regOn]]$labels[regIndicators == 1]
      argsIn$regIndicators <- regIndicators}
  }
  if(optimization == "approx" & is.character(regIndicators)){
    if(tolower(objective) == "kalman"){stop("Please provide regIndicators as matrix")}
    regIndicators <- regIndicatorsFromNameToMatrix(mxObject = mxObject, regOn = regOn, regIndicators = regIndicators)
    argsIn$regIndicators <- regIndicators
  }

  if(tolower(objective) == "kalman"){
    # create individual models
    # Note that we assume all persons to have the same parameter values
    if(autoCV){
      # if autoCV, we only need the model to extract the parameter labels. No need to use the entire data set
      argsIn$mxObject <- createKalmanMultiSubjectModel(ctsemObject = ctsemObject, dataset = t(as.matrix(dataset[1,])), useOptimizer = FALSE, silent = TRUE)
    }else{
      argsIn$mxObject <- createKalmanMultiSubjectModel(ctsemObject = ctsemObject, dataset = dataset, useOptimizer = FALSE, KalmanStartValues = KalmanStartValues)
    }
  }

  if(optimization == "exact"){
    regCtsemObject <- do.call(regCtsem::exact_regCtsem, rlist::list.remove(argsIn, c("optimization",
                                                                                      "epsilon",
                                                                                      "zeroThresh",
                                                                                      "elastic_alpha",
                                                                                      "elastic_gamma",
                                                                                      "calledInternally")))


    if(!autoCV){
      fitAndParametersSeparated <- try(separateFitAndParameters(regCtsemObject))
      if(!any(class(fitAndParametersSeparated) == "try-error")){
        regCtsemObject$fit <- fitAndParametersSeparated$fit
        regCtsemObject$parameterEstimatesRaw <- fitAndParametersSeparated$parameterEstimates
        regCtsemObject$parameters <- try(getVariancesInParameterEstimates(mxObject = regCtsemObject$setup$mxObject,
                                                                      parameterEstimates = fitAndParametersSeparated$parameterEstimates), silent = T)
        if(any(class(regCtsemObject$parameters) == "try-error")){
          warning("Could not compute the variances and covariances from the DIFFUSIONbase and T0VARbase. This is a bug and will be resolved later on.")
        }
      }
    }

    class(regCtsemObject) <- "regCtsem"
    return(regCtsemObject)
  }

  regCtsemObject <- do.call(regCtsem::approx_regCtsem, rlist::list.remove(argsIn, c("optimization", "sparseParameters", "optimizer", "tryCpptsem", "forceCpptsem", "stepSize", "lineSearch", "c1",
                                                                                     "c2", "differenceApprox", "initialHessianApproximation",
                                                                                     "maxIter_out", "maxIter_in",
                                                                                     "maxIter_line", "eps_out", "eps_in",
                                                                                     "eps_WW", "eta", "sig", "gam", "stepsizeMin", "stepsizeMax",
                                                                                     "maxIter_out", "maxIter_in", "GISTNonMonotoneNBack", "GISTLinesearchCriterion",
                                                                                     "break_outer",
                                                                                     "exactApproximateFirst", "exactApproximateFirst3NumStartingValues", "exactApproximateFirst3Optimize", "exactApproximateFirstMaxIter_out", "calledInternally")))
  if(!autoCV){
    fitAndParametersSeparated <- try(separateFitAndParameters(regCtsemObject))
    if(!any(class(fitAndParametersSeparated) == "try-error")){
      regCtsemObject$fit <- fitAndParametersSeparated$fit
      regCtsemObject$parameterEstimatesRaw <- fitAndParametersSeparated$parameterEstimates
      regCtsemObject$parameters <- getVariancesInParameterEstimates(mxObject = regCtsemObject$setup$mxObject, parameterEstimates = fitAndParametersSeparated$parameterEstimates)
    }
  }

  class(regCtsemObject) <- "regCtsem"
  cat("\n")
  return(regCtsemObject)
}



