#' regCtsem
#'
#' main function: performs regularized ctsem
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param ctsemObject if objective = "ML": Fitted object of class ctsem. If you want to use objective = "Kalman", pass an object of type ctsemInit from ctModel
#' @param dataset only required if objective = "Kalman" and ctsemObject is of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01). Alternatively, lambdas can be set to "auto"; however, this is still experimental.
#' @param lambdasAutoLength if lambdas == "auto", lambdasAutoLength will determine the number of lambdas tested.
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso and adaptiveLasso for optimization = exact
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the unregularized parameter estimates.
#' @param elastic_alpha placehoder for elastic net. NOT YET IMPLEMENTED
#' @param elastic_gamma placehoder for elastic net. NOT YET IMPLEMENTED
#' @param cvSample cross-validation sample. Has to be of type mxData
#' @param autoCV Boolean: Should automatic cross-validation be used?
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param KalmanStartValues Optional starting values for the parameters when using Kalman filter
#' @param optimizeKalman Boolen: Should the Kalman model be optimized in OpenMx first? If you want the Kalman model to start optimizing in regCtsem from the provided KalmanStartValues and not use OpenMx to optimize the initial Kalman model, set to FALSE
#' @param sparseParameters labeled vector with parameter estimates of the most sparse model. Required for approxFirst = 3
#' @param standardizeDrift Boolean: Should Drift parameters be standardized automatically using T0VAR?
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param optimization which optimization procedure should be used. Possible are  "exact" or "approx".
#' @param optimizer for exact optimization: Either GIST or GLMNET
#' @param control List with control arguments for the optimizer. See ?controlGIST, ?controlGLMNET and ?controlApprox for the respective parameters
#' @param extraTries number of extra tries in mxTryHard
#' @param cores how many computer cores should be used?
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @param silent silent execution
#' @param progressBar Boolean: Should a progress bar be displayed
#' @param parallelProgressBar list: used internally to display progress when executing in parallel. Do not pass values to parallelProgressBar
#' @param calledInternally Boolean: used internally for skipping checks
#'
#' @examples
#' set.seed(17046)
#'
#' library(regCtsem)
#'
#' #### Example 1 ####
#' ## Regularization with FIML objective function
#'
#' ## Population model:
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
#' dat <- ctsem::ctGenerate(generatingModel, n.subjects = 100, wide = TRUE)
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
#' # select DRIFT values for regularization:
#' regOn = "DRIFT"
#' regIndicators = matrix(c(0,1,
#'                          1,0),
#'                        byrow = T, nrow = 2)
#'
#' # Optimize model using GIST with lasso penalty
#' regModel <- regCtsem::regCtsem(ctsemObject = fit_myModel,
#'                                regOn = regOn,
#'                                regIndicators = regIndicators,
#'                                lambdas = "auto",
#'                                lambdasAutoLength = 15)
#'
#' summary(regModel, criterion = "BIC")
#' plot(regModel, what = "drift")
#' plot(regModel, what = "fit", criterion = c("AIC", "BIC", "m2LL"))
#'
#' # Using shinyfy can also give some insights into the model:
#' # shinyfy(regCtsemObject = regModel)
#'
#' # The best parameter estimates and the final model as mxObject can be extracted with:
#' # getFinalParameters(regCtsemObject = regModel, criterion = "BIC")
#' # bestModel <- getFinalModel(regCtsemObject = regModel, criterion = "BIC")
#'
#' # Optimize model using GLMNET with lasso penalty
#' regModel <- regCtsem::regCtsem(ctsemObject = fit_myModel,
#'                                regOn = regOn,
#'                                regIndicators = regIndicators,
#'                                lambdas = "auto",
#'                                lambdasAutoLength = 15,
#'                                optimizer = "GLMNET")
#'
#' summary(regModel, criterion = "BIC")
#' plot(regModel, what = "drift")
#' plot(regModel, what = "fit", criterion = c("AIC", "BIC", "m2LL"))
#'
#' # The same regularization can be performed with the approximate optimization:
#' # Note that we are using extraTries to get better parameter estimates
#' regModelApprox <- regCtsem::regCtsem(ctsemObject = fit_myModel,
#'                                      regOn = regOn,
#'                                      regIndicators = regIndicators,
#'                                      lambdas = "auto",
#'                                      lambdasAutoLength = 15,
#'                                      optimization = "approx",
#'                                      control = list(
#'                                        epsilon = .001, # epsilon is used to transform the non-differentiable
#'                                        #lasso penalty to a differentiable one if optimization = approx
#'                                        zeroThresh = .04 # threshold below which parameters will be evaluated as == 0
#'                                      ),
#'                                      extraTries = 5)
#'
#' # Comparison of parameter estimates:
#' round(regModel$fitAndParameters - regModelApprox$fitAndParameters,4)
#'
#' #### Example 2 ####
#' ## Regularization with Kalman objective function
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
#' ## Optimization with GIST:
#' regModel <- regCtsem::regCtsem(# important: We now have to pass the ctModel object, not the fitted model!
#'   ctsemObject = myModel,
#'   # Furthermore, the data has to be passed to regCtsem
#'   dataset = traindata,
#'   regOn = regOn,
#'   regIndicators = regIndicators,
#'   lambdas = "auto",
#'   lambdasAutoLength = 15,
#'   cvSample = testdata,
#'   objective = "Kalman",
#'   cores = 2
#' )
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
  lambdas = "auto",
  lambdasAutoLength = 50,
  penalty = "lasso",
  adaptiveLassoWeights = NULL,
  elastic_alpha = NULL,
  elastic_gamma = NULL,
  cvSample = NULL,
  autoCV = FALSE,
  k = 5,
  objective = "ML",
  KalmanStartValues = NULL,
  optimizeKalman = TRUE,
  sparseParameters = NULL,

  standardizeDrift = FALSE,
  scaleLambdaWithN = TRUE,
  returnFitIndices = TRUE,

  optimization = "exact",
  optimizer = "GIST",
  control = list(),
  extraTries = 3,
  # additional settings
  cores = 1,
  verbose = 0,
  silent = FALSE,
  progressBar = TRUE,
  parallelProgressBar = NULL,
  calledInternally = FALSE
){
  # save input in list
  argsIn <- as.list(environment())

  ## Defaults for optimizer
  if(optimization == "approx"){
    controlTemp <- controlApprox()
  } else if(optimization == "exact" && optimizer == "GIST"){
    controlTemp <- controlGIST()
  }else if(optimization == "exact" && optimizer == "GLMNET"){
    controlTemp <- controlGLMNET()
  }else{
    stop("Unknown optimization or optimizer specified. Possible are exact or approx for optimization and GIST or GLMNET for optimizer")
  }

  controlExpectedNames <- names(controlTemp)
  controlReceivedNames <- names(control)
  if(!all(controlReceivedNames %in% controlExpectedNames)){
    stop(paste0("Unknown parameters passed to control. Expected any of the following arguments: ", paste0(controlExpectedNames, collapse = ", ")))
  }
  controlTemp[controlReceivedNames] <- control[controlReceivedNames]

  argsIn <- c(argsIn, controlTemp)

  # check setup
  checkSetup(argsIn)
  if(argsIn$optimization == "exact"){
    if(argsIn$approxFirst == "auto"){
      if((any(lambdas == "auto")) | !is.null(argsIn$sparseParameters)){
        argsIn$approxFirst <- 3
      }else{
        argsIn$approxFirst <- 4
      }
    }

    if(argsIn$approxFirst == 3){
      if((!any(argsIn$lambdas == "auto")) & any(is.null(argsIn$sparseParameters))){
        warning("approxFirst = 3 requested. This requires either lambdas = 'auto' or the sparseParameters being passed to regCtsem. Setting approxFirst = 4")
        argsIn$approxFirst <- 4
      }
    }

    if(argsIn$approxMaxIt == "auto"){
      if(argsIn$approxFirst == 3){
        argsIn$approxMaxIt <- 10
      }else{
        argsIn$approxMaxIt <- 200
      }
    }
  }

  if(is.null(argsIn$mxObject)){
    argsIn$mxObject <- argsIn$ctsemObject$mxobj
  }

  if(is.logical(argsIn$regIndicators)){
    warning("regIndicators is a logical matrix. Interpreting TRUE as 1 and FALSE as 0.")
    argsIn$regIndicators <- matrix(as.numeric(argsIn$regIndicators),
                                   nrow = nrow(argsIn$regIndicators),
                                   ncol = ncol(argsIn$regIndicators))
  }

  if(argsIn$optimization == "exact" && is.numeric(argsIn$regIndicators)){
    if(tolower(argsIn$objective) == "kalman"){
      argsIn$regIndicators <- argsIn$ctsemObject[[argsIn$regOn]][argsIn$regIndicators == 1]
    }else{
      argsIn$regIndicators <- argsIn$mxObject[[argsIn$regOn]]$labels[argsIn$regIndicators == 1]
    }
  }
  if(argsIn$optimization == "approx" & is.character(argsIn$regIndicators)){
    if(tolower(argsIn$objective) == "kalman"){stop("Please provide regIndicators as matrix")}
    argsIn$regIndicators <- regIndicatorsFromNameToMatrix(mxObject = argsIn$mxObject, regOn = argsIn$regOn,
                                                          regIndicators = argsIn$regIndicators)
  }

  if(tolower(argsIn$objective) == "kalman"){
    # create individual models
    # Note that we assume all persons to have the same parameter values
    if(argsIn$autoCV){
      # if autoCV, we only need the model to extract the parameter labels. No need to use the entire data set
      argsIn$mxObject <- createKalmanMultiSubjectModel(ctsemObject = argsIn$ctsemObject, dataset = t(as.matrix(argsIn$dataset[1,])),
                                                       useOptimizer = FALSE, silent = TRUE)
    }else{
      argsIn$mxObject <- createKalmanMultiSubjectModel(ctsemObject = argsIn$ctsemObject, dataset = argsIn$dataset,
                                                       useOptimizer = FALSE, KalmanStartValues = argsIn$KalmanStartValues)
    }
  }

  if(argsIn$optimization == "exact"){
    regCtsemObject <- do.call(regCtsem::exact_regCtsem, rlist::list.remove(argsIn, c("optimization",
                                                                                     "control",
                                                                                     "elastic_alpha",
                                                                                     "elastic_gamma",
                                                                                     "calledInternally")))


    if(!argsIn$autoCV){
      fitAndParametersSeparated <- try(separateFitAndParameters(regCtsemObject))
      if(!any(class(fitAndParametersSeparated) == "try-error")){
        regCtsemObject$fit <- fitAndParametersSeparated$fit
        regCtsemObject$parameterEstimatesRaw <- fitAndParametersSeparated$parameterEstimates
        regCtsemObject$parameters <- try(getVariancesInParameterEstimates(mxObject = regCtsemObject$setup$mxObject,
                                                                          parameterEstimates = fitAndParametersSeparated$parameterEstimates),
                                         silent = T)
        if(any(class(regCtsemObject$parameters) == "try-error")){
          warning("Could not compute the variances and covariances from the DIFFUSIONbase and T0VARbase. This is a bug and will be resolved later on.")
        }
      }
    }

    class(regCtsemObject) <- "regCtsem"
    return(regCtsemObject)
  }

  regCtsemObject <- do.call(regCtsem::approx_regCtsem, rlist::list.remove(argsIn, c("optimization", "sparseParameters", "optimizer", "control", "calledInternally")))
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


#' exact_regCtsem
#'
#' creates a regCtsem object for exact optimization
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param ctsemObject if objective = "ML": Fitted object of class ctsem. If you want to use objective = "Kalman", pass an object of type ctsemInit from ctModel
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param dataset only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param lambdasAutoLength if lambdas == "auto", lambdasAutoLength will determine the number of lambdas tested.
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
#' @param sparseParameters labeled vector with parameter estimates of the most sparse model. Required for approxFirst = 3
#' @param tryCpptsem should regCtsem try to translate the model to cpptsem? This can speed up the computation considerably but might fail for some models
#' @param forceCpptsem should cpptsem be enforced even if results differ from ctsem? Sometimes differences between cpptsem and ctsem can result from problems with numerical precision which will lead to the m,atrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "cpptsem") for more details
#' @param stepSize GLMNET & GIST: initial step size of the outer iteration
#' @param lineSearch GLMNET: String indicating if Wolfe conditions (lineSearch = "Wolfe") should be used in the outer iteration. Alternatively, set to "none" (not recommended!).
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
#' @param approxFirst Should approximate optimization be used first to obtain start values for exact optimization? 1 = only for first lambda, 2 = for all lambdas, 3 = ensures that the fit will not be worse than in the sparse model if lambdas = "auto" or sparseParameters are provided. To this end, 10 models between the current parameter estimates and the sparse parameter estimates are tested and the one with the lowest regM2LL is used for starting values.
#' @param numStart Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param approxOpt Used if approxFirst = 3. Should each of the generated starting values be optimized slightly? This can substantially improve the fit of the generated starting values. 1 = optimization with optim, 2 = optimization with Rsolnp
#' @param approxMaxIt Used if approxFirst = 3 and approxOpt > 1. How many outer iterations should be given to each starting values vector? More will improve the selected starting values but slow down the computation. If approxFirst =  4, or approxFirst = 5 this will control the number of outer iteration in optim or solnp .
#' @param extraTries number of extra tries in mxTryHard for warm start
#' @param cores how many computer cores should be used?
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @param silent silent execution
#' @param progressBar Boolean: Should a progress bar be displayed
#' @param parallelProgressBar list: used internally to display progress when executing in parallel
#'
#' @author Jannik Orzek
#' @import ctsemOMX rlist
#' @export
exact_regCtsem <- function(  # model
  ctsemObject = NULL,
  mxObject = NULL,
  dataset = NULL,
  # penalty settings
  regOn = "DRIFT",
  regIndicators,
  lambdas,
  lambdasAutoLength = 50,
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
  sig = .2,
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
  break_outer = .0000000001,
  # general
  scaleLambdaWithN = TRUE,
  approxFirst = 0,
  numStart = 10,
  approxOpt = T,
  approxMaxIt = 5,
  extraTries = 3,
  # additional settings
  cores = 1,
  verbose = 0,
  silent = FALSE,
  progressBar = TRUE,
  parallelProgressBar = NULL){

  if(is.null(mxObject)){
    mxObject <- ctsemObject$mxobj
  }

  exactArgsIn <- as.list(environment())

  returnList <- list("setup" = exactArgsIn)

  parameterLabels <- names(OpenMx::omxGetParameters(mxObject))

  returnList$setup$parameterLabels <- parameterLabels

  # set up automatic cross-validation
  if(autoCV){

    # fit model
    resultsCV <- iterateOverCVFolds(exactArgsIn, objective = objective, optimization = "exact")
    fit <- resultsCV$fit
    returnList$setup$lambdas <- resultsCV$folds$models$fold1$setup$lambdas

    fit["mean CV fit",] <- apply(fit,2,mean,na.rm =TRUE)
    if(any(is.na(fit))){
      warning("NAs in fit results. Consider rerunning with different starting values.")
    }

    returnList <- rlist::list.append(returnList, "fit" = fit, "folds" = resultsCV$folds)
    return(returnList)
  }

  ##### no automatic cross-validation ####

  # if Kalman: fit initial model
  if(tolower(objective) == "kalman"){
    if(optimizeKalman){
      mxObject <- OpenMx::mxTryHardctsem(mxObject, silent = TRUE, extraTries = extraTries)
      exactArgsIn$KalmanStartValues <- omxGetParameters(mxObject)
      exactArgsIn$optimizeKalman <- FALSE
    }else{
      mxObject <- OpenMx::mxRun(mxObject, silent = TRUE, useOptimizer = FALSE)
    }
    exactArgsIn$mxObject <- mxObject
    sampleSize <- nrow(dataset)
  }else{
    sampleSize <- mxObject$data$numObs
  }

  # set adaptiveLassoWeights

  adaptiveLassoWeights <- getAdaptiveLassoWeights(mxObject = mxObject, penalty = penalty, adaptiveLassoWeights = adaptiveLassoWeights, standardizeDrift = standardizeDrift)
  returnList$setup$adaptiveLassoWeights <- adaptiveLassoWeights

  # set lambdas if lambdas == "auto"
  if(any(lambdas == "auto")){
    maxLambda <- getMaxLambda(mxObject = mxObject,
                                  regIndicators = regIndicators,
                                  differenceApprox = differenceApprox,
                                  adaptiveLassoWeights = adaptiveLassoWeights)
    sparseParameters <- maxLambda$sparseParameters
    maxLambda <- maxLambda$maxLambda
    lambdas <- seq(0,maxLambda, length.out = lambdasAutoLength)
    # if scaleLambdaWithN = TRUE, the lambda values will be multiplied with N later on
    # we need to ensure that this will not change the values:
    if(scaleLambdaWithN){
      lambdas <- lambdas/sampleSize
    }
    returnList$setup$lambdas <- lambdas
    exactArgsIn$lambdas <- lambdas
    returnList$setup$sparseParameters <- sparseParameters
    exactArgsIn$sparseParameters <- sparseParameters
  }

  if(cores == 1){
    fits <- c("regM2LL", "m2LL")
    if(returnFitIndices){
      fits <- c(fits, "AIC", "BIC", "estimatedParameters")
    }
    if(!is.null(cvSample)){
      fits <- c(fits, "cvM2LL")
    }

    # prepare table with results
    fitAndParameters <- matrix(NA,
                               nrow = length(c(fits, parameterLabels)),
                               ncol = length(lambdas),
                               dimnames = list(c(fits, parameterLabels),
                                               lambdas))

    # start fitting models
    if(tolower(optimizer) == "glmnet"){
      regModel <- try(regCtsem::exact_bfgsGLMNET(ctsemObject = ctsemObject, dataset = dataset, objective = objective,
                                                 mxObject = mxObject, regOn = regOn, regIndicators = regIndicators, lambdas = lambdas,
                                                 adaptiveLassoWeights = adaptiveLassoWeights,
                                                 # additional settings
                                                 sparseParameters = sparseParameters,
                                                 tryCpptsem = tryCpptsem, forceCpptsem = forceCpptsem, stepSize = stepSize, lineSearch = lineSearch, c1 = c1, c2 = c2,
                                                 sig = sig, gam = gam,
                                                 differenceApprox = differenceApprox,
                                                 initialHessianApproximation = initialHessianApproximation,
                                                 maxIter_out = maxIter_out,
                                                 maxIter_in = maxIter_in, maxIter_line = maxIter_line,
                                                 eps_out = eps_out, eps_in = eps_in, eps_WW = eps_WW,
                                                 scaleLambdaWithN = scaleLambdaWithN, sampleSize = sampleSize,
                                                 approxFirst = approxFirst,
                                                 numStart = numStart, approxOpt = approxOpt, approxMaxIt = approxMaxIt,
                                                 extraTries = extraTries,
                                                 verbose = verbose, progressBar = progressBar, parallelProgressBar = parallelProgressBar))
    }
    if(tolower(optimizer) == "gist"){
      regModel <- try(regCtsem::exact_GIST(ctsemObject = ctsemObject, dataset = dataset, objective = objective,
                                           mxObject = mxObject, regOn = regOn, regIndicators = regIndicators, lambdas = lambdas,
                                           adaptiveLassoWeights = adaptiveLassoWeights,
                                           # additional settings
                                           sparseParameters = sparseParameters,
                                           tryCpptsem = tryCpptsem, forceCpptsem = forceCpptsem, eta = eta, sig = sig,
                                           initialStepsize = stepSize, stepsizeMin = stepsizeMin, stepsizeMax = stepsizeMax,
                                           GISTLinesearchCriterion = GISTLinesearchCriterion, GISTNonMonotoneNBack = GISTNonMonotoneNBack,
                                           maxIter_out = maxIter_out, maxIter_in = maxIter_in,
                                           break_outer = break_outer,
                                           scaleLambdaWithN = scaleLambdaWithN, sampleSize = sampleSize,
                                           approxFirst = approxFirst,
                                           numStart = numStart, approxOpt = approxOpt, approxMaxIt = approxMaxIt,
                                           extraTries = extraTries,
                                           differenceApprox = differenceApprox,
                                           verbose = verbose,
                                           progressBar = progressBar, parallelProgressBar = parallelProgressBar))
    }
    if(!any(class(regModel) == "try-error")){
      # save results
      fitAndParameters[parameterLabels,] <- regModel$thetas
      fitAndParameters["m2LL",] <- regModel$m2LL
      fitAndParameters["regM2LL",] <- regModel$regM2LL
    }

    # save fit indices
    if(returnFitIndices & !(any(class(regModel) == "try-error"))){
      fitIndicesTable <- regCtsem::exact_getFitIndices(mxObject = mxObject, parameterLabels = parameterLabels,
                                                       fitAndParameters = fitAndParameters, lambdas = lambdas, sampleSize = sampleSize)

      fitAndParameters[rownames(fitIndicesTable),] <- fitIndicesTable[rownames(fitIndicesTable),]
    }


    # cross-validation
    if(!is.null(cvSample)){

      fitCVTable <- regCtsem::exact_getCVFit(objective = objective, ctsemObject = ctsemObject, mxObject = mxObject, parameterLabels = parameterLabels,
                                             parameterValuesTable = as.matrix(fitAndParameters[parameterLabels,], ncol = length(lambdas)),
                                             lambdas = lambdas, cvSample = cvSample)

      fitAndParameters["cvM2LL",] <- fitCVTable["cvM2LL",]
    }

    returnList <- rlist::list.append(returnList, "fitAndParameters" = fitAndParameters)
    return(returnList)
  }

  # multi core

  fit <- try(do.call(regCtsem::exact_parallelRegCtsem, exactArgsIn))

  returnList <- c(returnList, fit)
  return(returnList)

}

#' exact_parallelRegCtsem
#'
#' parallel processing for regCtsem with optimization = "exact"
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param ctsemObject if objective = "ML": Fitted object of class ctsem. If you want to use objective = "Kalman", pass an object of type ctsemInit from ctModel
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param dataset only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param lambdasAutoLength if lambdas == "auto", lambdasAutoLength will determine the number of lambdas tested.
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
#' @param approxFirst Should approximate optimization be used first to obtain start values for exact optimization? 1 = only for first lambda, 2 = for all lambdas
#' @param numStart Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param approxOpt Used if approxFirst = 3. Should each of the generated starting values be optimized slightly? This can substantially improve the fit of the generated starting values. 1 = optimization with optim, 2 = optimization with Rsolnp
#' @param approxMaxIt Used if approxFirst = 3 and approxOpt > 1. How many outer iterations should be given to each starting values vector? More will improve the selected starting values but slow down the computation. If approxFirst =  4, or approxFirst = 5 this will control the number of outer iteration in optim or solnp .
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

exact_parallelRegCtsem <- function(# model
  ctsemObject = NULL,
  mxObject = NULL,
  dataset = NULL,
  # penalty settings
  regOn = "DRIFT",
  regIndicators,
  lambdas,
  lambdasAutoLength = 50,
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
  approxFirst = 3,
  numStart = 10,
  approxOpt = T,
  approxMaxIt = 5,
  extraTries = 3,
  # additional settings
  cores = 1,
  verbose = 0,
  silent = FALSE,
  progressBar = TRUE,
  parallelProgressBar = NULL){

  # catch function call
  fun_call <- as.list(environment())

  ## extract control arguments
  if(fun_call$optimizer == "GIST"){
    expectedControlNames <- names(controlGIST())
  }else if(fun_call$optimizer == "GLMNET"){
    expectedControlNames <- names(controlGLMNET())
  }

  control <- fun_call[expectedControlNames]
  fun_call <- fun_call[!names(fun_call) %in% unique(names(c(controlGIST(), controlGLMNET())))]

  parameterLabels <- names(OpenMx::omxGetParameters(mxObject))
  fitLabels <- c("m2LL", "regM2LL", "AIC", "BIC", "cvM2LL", "estimatedParameters")


  # set up multicore execution

  # returns a global cluster variable
  regCtsemClusterVariable <<- snow::makeSOCKcluster(cores, outfile="")
  doSNOW::registerDoSNOW(regCtsemClusterVariable)


  # Split the lambdas in equal lambdaBins:

  lambdaBins <- split(lambdas,
                        cut(x = 1:length(lambdas),
                            breaks = cores,
                            labels = paste0("core", 1:cores)))

  # for progess printing
  parallelTempFiles <- vector(mode = "list", length = cores)
  names(parallelTempFiles) <- paste0("core", 1:cores)

  for(core in 1:cores){
    parallelTempFiles[core] <- tempfile(pattern = paste0("core", core), fileext = ".txt")
    write.csv2(0,parallelTempFiles[[core]], row.names = FALSE)
  }

  maxItSum <- length(lambdas)

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


  fitAndParameters_combined <- foreach::foreach(iteration = 1:length(lambdaBins), .combine = cbind,
                                                .multicombine = TRUE,.packages = c("OpenMx","regCtsem"),
                                                .inorder = TRUE,
                                                .errorhandling = "remove",
                                                .verbose = FALSE
  ) %dopar% {
    # lambdas in this core
    lambdaBin <- lambdaBins[[iteration]]

    # expected results
    fitTable <- matrix(NA, nrow = length(parameterLabels)+length(fitLabels), ncol = length(lambdaBin))
    rownames(fitTable) <- c(parameterLabels, fitLabels)
    colnames(fitTable) <- lambdaBin

    initialModel <- mxObject

    iterationExactArgsIn <- fun_call
    iterationExactArgsIn$mxObject <- initialModel
    iterationExactArgsIn$lambdas <- lambdaBin
    iterationExactArgsIn$cores <- 1
    iterationExactArgsIn$progressBar <- FALSE
    iterationExactArgsIn$parallelProgressBar = list("printProgress" = printProgress,
                                                    "parallelTempFiles" = parallelTempFiles,
                                                    "maxItSum" = maxItSum,
                                                    "cores" = cores,
                                                    "iteration" = iteration)
    iterationExactArgsIn$control <- control

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



#' approx_regCtsem
#'
#' creates a regCtsem object for approximate optimization
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param ctsemObject Fitted object of class ctsemOMX
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param dataset only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param lambdasAutoLength if lambdas == "auto", lambdasAutoLength will determine the number of lambdas tested.
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
  lambdas,
  lambdasAutoLength = 50,
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

  # set lambdas if lambdas == "auto"
  if(any(lambdas == "auto")){
    regIndicatorsString <- mxObject[[regOn]]$labels[regIndicators==1]
    maxLambda <- getMaxLambda(mxObject = mxObject,
                                  regIndicators = regIndicatorsString,
                                  differenceApprox = "central",
                                  adaptiveLassoWeights = adaptiveLassoWeights)
    lambdas <- seq(0,maxLambda$maxLambda, length.out = lambdasAutoLength)
    # if scaleLambdaWithN = TRUE, the lambda values will be multiplied with N later on
    # we need to ensure that this will not change the values:
    if(scaleLambdaWithN){
      lambdas <- lambdas/sampleSize
    }
    returnList$setup$lambdas <- lambdas
    approxArgsIn$lambdas <- lambdas
  }

  if(cores == 1){

    # start fitting models
    fitAndParameters <- try(regCtsem::approx_iterateOverLambdas(ctsemObject = ctsemObject, mxObject = mxObject, sampleSize = sampleSize,
                                                                  # penalty settings
                                                                  regOn = regOn, regIndicators = regIndicators, lambdas = lambdas,
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



#' approx_parallelRegCtsem
#'
#' parallel processing for regCtsem with optimization = "approx"
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param ctsemObject Fitted object of class ctsemOMX
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @param lambdasAutoLength if lambdas == "auto", lambdasAutoLength will determine the number of lambdas tested.
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

approx_parallelRegCtsem <- function(
  # model
  ctsemObject = NULL,
  mxObject = NULL,
  dataset = NULL,
  # penalty settings
  regOn = "DRIFT",
  regIndicators,
  lambdas,
  lambdasAutoLength = 50,
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


  # Split the lambdas in equal lambdaBins:

  lambdaBins <- split(lambdas,
                        cut(x = 1:length(lambdas),
                            breaks = cores,
                            labels = paste0("core", 1:cores)))

  # for progress printing
  parallelTempFiles <- vector(mode = "list", length = cores)
  names(parallelTempFiles) <- paste0("core", 1:cores)

  for(core in 1:cores){
    parallelTempFiles[core] <- tempfile(pattern = paste0("core", core), fileext = ".txt")
    write.csv2(0,parallelTempFiles[[core]], row.names = FALSE)
  }

  maxItSum <- length(lambdas)

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


  fitAndParameters_combined <- foreach::foreach(iteration = 1:length(lambdaBins), .combine = cbind,
                                                .multicombine = TRUE,.packages = c("OpenMx","regCtsem"),
                                                .inorder = TRUE,
                                                .errorhandling = "remove",
                                                .verbose = FALSE
  ) %dopar% {
    # lambdas in this core
    lambdaBin <- lambdaBins[[iteration]]

    # set up model
    iterationParallelProgressBar = list("printProgress" = printProgress,
                                        "parallelTempFiles" = parallelTempFiles,
                                        "maxItSum" = maxItSum,
                                        "cores" = cores,
                                        "iteration" = iteration)

    iterationResult <- approx_regCtsem(ctsemObject = ctsemObject, mxObject = mxObject, dataset = dataset,
                                       # penalty settings
                                       regOn = regOn, regIndicators = regIndicators, lambdas = lambdaBin,
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


#### Cross - Validation ####

#' createCVFitTable
#'
#' sets up a table for CV fit values
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01)
#' @author Jannik Orzek
#' @export
createCVFitTable <- function(k, lambdas){
  CVFitTable <- matrix(NA, nrow = k+1, ncol = length(lambdas))
  rownames(CVFitTable) <- c(paste0("fold", 1:k), "mean CV fit")
  colnames(CVFitTable) <- lambdas
  return(CVFitTable)
}




#' createCVFoldsAndModels
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject.
#' @param dataset only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @author Jannik Orzek
#' @export
createCVFoldsAndModels <- function(mxObject, dataset = NULL, k){
  if(is.null(dataset)){
    fullData <- mxObject$data$observed}else{
      fullData <- dataset
    }

  cvFolds <- regCtsem::createFolds(mxObject = mxObject, dataset = dataset, k = k)

  testSets <- vector("list", length = k)
  names(testSets) <- names(cvFolds)

  trainSets <- vector("list", length = k)
  names(trainSets) <- names(cvFolds)

  cvModels <- vector("list", length = k)
  names(cvModels) <- names(cvFolds)

  return(list("fullData" = fullData, "cvFolds" = cvFolds, "testSets" = testSets, "trainSets" = trainSets, "cvModels" = cvModels))
}




#' createFolds
#'
#' created folds for automatic cross-validation
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param mxObject Object of class mxModel
#' @param dataset only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param k number of cross-validation folds (k-fold cross-validation)
#' @author Jannik Orzek
#' @export
createFolds = function(mxObject, dataset = NULL, k){
  if(is.null(dataset)){
    sampleSize <- mxObject$data$numObs
  }else{
    sampleSize <- nrow(dataset)
  }

  # shuffle
  Folds <- sample(1:sampleSize, sampleSize)

  Folds <- split(Folds, cut(x = 1:length(Folds), breaks = k, labels = paste0("fold", 1:k)))

  return(Folds)
}




#' iterateOverCVFolds
#'
#' computes results for automativ cross-validation
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param argsIn list of parameters passed to regCtsem
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param optimization type of optimization. Either exact or approx
iterateOverCVFolds <- function(argsIn, objective = "ML", optimization){
  # create folds
  cvFoldsAndModels <- createCVFoldsAndModels(mxObject = argsIn$mxObject, dataset = argsIn$dataset, k = argsIn$k)
  cvFolds <- cvFoldsAndModels$cvFolds
  fullData <- cvFoldsAndModels$fullData
  testSets <- cvFoldsAndModels$testSets
  trainSets <- cvFoldsAndModels$trainSets
  cvModels <-  cvFoldsAndModels$cvModels

  # get lambdas if lambdas == "auto"
  if(any(argsIn$lambdas == "auto")){
    maxLambda <- getMaxLambdaCV(ctsemObject = argsIn$ctsemObject,
                                    mxObject = argsIn$mxObject,
                                    fullData = fullData,
                                    trainSets = trainSets,
                                    KalmanStartValues = argsIn$KalmanStartValues,
                                    regOn = argsIn$regOn,
                                    regIndicators = argsIn$regIndicators,
                                    penalty = argsIn$penalty,
                                    adaptiveLassoWeights = argsIn$adaptiveLassoWeights,
                                    standardizeDrift = argsIn$standardizeDrift,
                                    k = argsIn$k,
                                    cvFolds = cvFolds,
                                    scaleLambdaWithN = argsIn$scaleLambdaWithN,
                                    objective = objective,
                                    optimization = optimization,
                                    differenceApprox = ifelse(is.null(argsIn$differenceApprox), "central", argsIn$differenceApprox)
    )
    if(optimization == "exact"){
      sparseParameterMatrix <- maxLambda$sparseParameterMatrix
    }
    maxLambda <- maxLambda$maxLambda
    argsIn$lambdas <- seq(0, maxLambda, length.out = argsIn$lambdasAutoLength)
  }

  # create fit table
  fit <- createCVFitTable(k = argsIn$k, lambdas = argsIn$lambdas)

  for(foldNumber in 1:argsIn$k){
    if(is.vector(fullData[cvFolds[[foldNumber]],])){
      # if a single row is selected, R extracts this row as vector, not as matrix
      testSets[[foldNumber]] <-  t(as.matrix(fullData[cvFolds[[foldNumber]],]))
    }else{
      testSets[[foldNumber]] <-  fullData[cvFolds[[foldNumber]],]
    }
    if(is.vector(fullData[-cvFolds[[foldNumber]],])){
      # if a single row is selected, R extracts this row as vector, not as matrix
      trainSets[[foldNumber]] <-  t(as.matrix(fullData[-cvFolds[[foldNumber]],]))
    }else{
      trainSets[[foldNumber]] <- fullData[-cvFolds[[foldNumber]],]
    }

    if(tolower(objective) == "ml"){
      currentModel <- argsIn$mxObject
      currentModel$data <- OpenMx::mxData(trainSets[[foldNumber]], type = "raw")
      currentModel <- mxRun(currentModel, useOptimizer = F, silent = T)

      # set input arguments
      currentModelArgsIn <- argsIn
      currentModelArgsIn$ctsemObject$mxobj <- currentModel
      currentModelArgsIn$mxObject <- currentModel
      currentModelArgsIn$autoCV <- FALSE
      currentModelArgsIn$cvSample <- OpenMx::mxData(testSets[[foldNumber]], type = "raw")
      currentModelArgsIn$returnFitIndices <- FALSE

      if(optimization == "exact"){
        if(exists("sparseParameterMatrix")){
          sparseParameterLabels <- rownames(sparseParameterMatrix)
          currentModelArgsIn$sparseParameters <- sparseParameterMatrix[,foldNumber]
          names(currentModelArgsIn$sparseParameters) <- sparseParameterLabels
        }
        currentModelFit <- try(do.call(regCtsem::exact_regCtsem, currentModelArgsIn))
      }else if(optimization == "approx"){
        currentModelFit <- try(do.call(regCtsem::approx_regCtsem, currentModelArgsIn))
      }
    }else if(tolower(objective) == "kalman"){
      # set input arguments
      currentModelArgsIn <- argsIn
      currentModelArgsIn$autoCV <- FALSE
      currentModelArgsIn$dataset <- trainSets[[foldNumber]]
      currentModelArgsIn$cvSample <- testSets[[foldNumber]]
      currentModelArgsIn$returnFitIndices <- FALSE
      currentModelArgsIn$mxObject <- createKalmanMultiSubjectModel(ctsemObject = argsIn$ctsemObject,
                                                                   dataset = trainSets[[foldNumber]],
                                                                   useOptimizer = FALSE,
                                                                   KalmanStartValues = argsIn$KalmanStartValues)

      if(optimization == "exact"){
        if(exists("sparseParameterMatrix")){
          sparseParameterLabels <- rownames(sparseParameterMatrix)
          currentModelArgsIn$sparseParameters <- sparseParameterMatrix[,foldNumber]
          names(currentModelArgsIn$sparseParameters) <- sparseParameterLabels
        }
        currentModelFit <- try(do.call(exact_regCtsem, currentModelArgsIn))
      }else if(optimization == "approx"){
        currentModelFit <- try(do.call(approx_regCtsem, currentModelArgsIn))
      }
    }

    cvModels[[foldNumber]] <- currentModelFit

    if(!any(class(currentModelFit) == "try-error")){
      if(any(!(colnames(currentModelFit$fitAndParameters) %in% colnames(fit) ))){
        warning("Error when binding cv results. Returing results as global variable currentModelFitError")
        currentModelFitError <<- currentModelFit
        next
      }
      if(! "cvM2LL" %in% rownames(currentModelFit$fitAndParameters)){
        warning("Error when binding cv results. Returing results as global variable currentModelFitError")
        currentModelFitError <<- currentModelFit
        next
      }
      fit[foldNumber,colnames(currentModelFit$fitAndParameters)] <- currentModelFit$fitAndParameters["cvM2LL",colnames(currentModelFit$fitAndParameters)]
      cat("\n")
      print(paste("Finished CV", foldNumber, "of", argsIn$k))
    }
  }
  return(list("fit" = fit, "folds" = list("models" = cvModels, "members" = cvFolds)))
}


#### Helper Functions ####

#' checkSetup
#'
#' internal checks
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param argsIn list with arguments
#' @author Jannik Orzek
#' @export
checkSetup <- function(argsIn){
  if(is.null(argsIn$ctsemObject) & is.null(argsIn$mxObject)){
    stop("Both ctsemObject and mxObect are missing. You need to provide at least one")
  }

  if(any(class(argsIn$ctsemObject)=="ctsemFit")){
    # check if ctsemObject is estimated with Kalman
    if(argsIn$ctsemObject$ctfitargs$objective == "Kalman"){
      stop("It seems like the provided ctsemObject was fitted with Kalman filter. To use the Kalman filter, provide the object of type ctsemInit from ctModel instead of the fitted model. Set the objective = 'Kalman' and provide a dataset")
    }
  }
  if(!is.null(argsIn$mxObject) & any(class(argsIn$mxObject)=="MxModel")){
    # check if ctsemObject is estimated with Kalman
    if(any(class(argsIn$mxObject$expectation) == "MxExpectationStateSpace") &  !any(class(argsIn$ctsemObject)=="ctsemInit")){
      stop("It seems like the provided mxObject was fitted with Kalman filter. To use the Kalman filter, provide the object of type ctsemInit from ctModel instead of the fitted model. Set the objective = 'Kalman' and provide a dataset")
    }
  }

  if(!(argsIn$optimization == "exact" || argsIn$optimization == "approx")){
    stop(paste("Optimization was set to", optimization, "however only exact or approx are supported."))
  }
  if(argsIn$optimization == "approx" && argsIn$extraTries == 1){
    warning(paste("Approximate optimization often requires multiple tries to find the optimal solution. Use extraTries to automatically try different statring values (e.g., extraTries = 5)."))
  }
  if(argsIn$cores > 1 && argsIn$verbose>0){
    stop("verbose > 0 only possible for single core execution. Set cores = 1 or verbose = 0.")
  }

  if(argsIn$optimization == "exact" && argsIn$optimizer == "GLMNET" && tolower(argsIn$lineSearch) == "armijo"){
    stop("lineSearch = 'armijo' is deprecated. Use lineSearch = 'Wolfe'.")
  }

  if(tolower(argsIn$objective) == "kalman"){
    if(is.null(argsIn$ctsemObject)){stop("Kalman filter requires a ctsemObject. Set up the model in ctsemOMX first")}
    if(!any(class(argsIn$ctsemObject) == "ctsemInit")){stop("ctsemObject has to be of class ctsemInit. You can get the ctsemInit object from ctModel.")}
    if(is.null(argsIn$dataset)){stop("No dataset was provided")}
  }

  if(tolower(argsIn$penalty) == "ridge"){
    if(argsIn$standardizeDrift){
      stop("Automatic drift standardization not supported for ridge regularization")
    }
    if(argsIn$optimization == "exact"){
      stop("Use optimization = approx for ridge regularization. This is confusing, but approx is the right optimization for ridge regularized SEM. The approx and exact differentiation is only relevant for lasso and adaptivelasso.")
    }
  }

  if(tolower(argsIn$penalty) == "adaptivelasso"){
    if(is.null(argsIn$adaptiveLassoWeights)){
      stop("Adaptive Lasso requires the adaptiveLassoWeights to be specified. auto uses the inverse of the absolute values of unregularized parameter estimates. This is only recommended with standardizeDrift = FALSE")
    }
    if(argsIn$standardizeDrift){
      stop("Combination of automatic standardization and adaptive lasso is not implemented. Standardization is a special variant of the adaptive lasso. Use penalty = 'lasso' or standardizeDrift = FALSE.")
    }
  }
  #if(any(argsIn$lambdas == "auto") & (tolower(argsIn$penalty) == "adaptivelasso" | argsIn$standardizeDrift)){
  #  stop("lambdas = 'auto' currently not supported for adative lasso or automatic standardization of drift parameters.")
  #}
}



#' createKalmanMultiSubjectModel
#'
#' Creates an mxModel Object for N >= 1 individuals using the Kalman filter
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param ctsemObject object of type ctsemInit from ctModel
#' @param dataset data set in wide format compatible to ctsem
#' @param useOptimizer Boolean: should the model be optimized
#' @export
createKalmanMultiSubjectModel <- function(ctsemObject, dataset, useOptimizer, silent = FALSE, KalmanStartValues = NULL){


  suppressMessages(invisible(capture.output(fit_kalmanModels <- ctsemOMX::ctFit(ctmodelobj = ctsemObject,
                                                                                dat = dataset,
                                                                                objective = "Kalman",
                                                                                fit = useOptimizer))))
  mxObject <- fit_kalmanModels$mxobj
  if(!is.null(KalmanStartValues)){
    parameterLabels <- names(OpenMx::omxGetParameters(mxObject))
    if(!all(names(KalmanStartValues)%in%parameterLabels)){
      stop("KalmanStartValues must have the same labels as the parameters in the model.")
    }
    mxObject <- OpenMx::omxSetParameters(mxObject, labels = parameterLabels, values = KalmanStartValues[parameterLabels])
  }
  return(mxObject)

  ######## Not used #######
  # create individual models
  # Note that we assume all persons to have the same parameter values
  individualModels <- vector("list", length = nrow(dataset))
  individualModelNames <- paste0("person", 1:nrow(dataset))

  for(person in 1:nrow(dataset)){
    if(!silent){
      cat('\r',paste0("Setting up the Kalman model: ", person, " of ", nrow(dataset)))
      flush.console()}
    suppressMessages(invisible(capture.output(individualModels[[person]] <- OpenMx::mxModel(name = individualModelNames[person],
                                                                                            ctsemOMX::ctFit(dat = t(as.matrix(dataset[person,])),
                                                                                                            ctmodelobj = ctsemObject,
                                                                                                            objective = 'Kalman',
                                                                                                            useOptimizer = FALSE)$mxobj))))
  }

  pointersToMatricesAndAlgebras <- vector("list", length = (length(individualModels[[1]]$algebras) + length(individualModels[[1]]$matrices)))
  namesOfMatricesAndAlgebras <- c(names(individualModels[[1]]$algebras), names(individualModels[[1]]$matrices))
  names(pointersToMatricesAndAlgebras) <- namesOfMatricesAndAlgebras
  for(i in 1:length(namesOfMatricesAndAlgebras)){
    pointersToMatricesAndAlgebras[[i]] <- mxAlgebraFromString(name = namesOfMatricesAndAlgebras[i], algString = paste0("person1.", namesOfMatricesAndAlgebras[i]))
  }

  mxIndividualModels <- OpenMx::mxModel(name = "MultiModel",
                                        submodels  = individualModels,
                                        # OpenMx::mxData(dataset, type = "raw"),
                                        OpenMx::mxFitFunctionMultigroup(individualModelNames),
                                        # for easy access to the matrices and algebras:
                                        pointersToMatricesAndAlgebras
  )

  mxIndividualModels <- OpenMx::omxAssignFirstParameters(mxIndividualModels)

  fit.mxIndividualModels <- mxRun(mxIndividualModels, silent = TRUE, useOptimizer = useOptimizer)
  cat("\n")
  return(fit.mxIndividualModels)
}




#' getAdaptiveLassoWeights
#'
#' Computes the adaptive lasso weights
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso for optimization = exact
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param standardizeDrift Boolean: Should Drift parameters be standardized automatically using T0VAR?
#' @author Jannik Orzek
#' @export
getAdaptiveLassoWeights <- function(mxObject, penalty, adaptiveLassoWeights, standardizeDrift){

  # if adaptiveLassoWeights were provided and no standardization is requested:
  if(tolower(penalty) == "adaptivelasso" & is.numeric(adaptiveLassoWeights) & (!standardizeDrift)){
    return(adaptiveLassoWeights)
  }

  # otherwise: set adaptiveLassoWeights
  ## if automatic drift standardization was requested
  if(standardizeDrift){
    # check if mxObject was optimized
    ## TODO

    # check if adaptiveLassoWeights were provided
    if(is.numeric(adaptiveLassoWeights)){ stop("standardizeDrift and provided adaptiveLassoWeights can not be combined automatically. Consider setting adaptiveLassoWeights = 'auto'") }

    # check if lasso or adaptiveLasso
    if(tolower(penalty) == "lasso"){
      thetaNames <- names(OpenMx::omxGetParameters(mxObject))
      adaptiveLassoWeights <- rep(1,length(thetaNames))
      names(adaptiveLassoWeights) <- thetaNames

      # compute standardizers
      T0VAR <- getT0VAR(mxObject)
      driftLabels <- mxObject$DRIFT$labels
      if(anyNA(driftLabels)){
        autoDriftLabels <- matrix(paste0("_autoDriftLabel_", rep(seq_len(nrow(driftLabels)), each = ncol(driftLabels)), "_", seq_len(ncol(driftLabels))),
                                  nrow = nrow(driftLabels), ncol = ncol(driftLabels), byrow = T)
        driftLabels[is.na(driftLabels)] <- autoDriftLabels[is.na(driftLabels)]
      }
      flatStandardizers <- regCtsem::getFlatStdizer(T0VAR = T0VAR, driftLabels = driftLabels)
      flatStandardizers <- flatStandardizers[rownames(flatStandardizers) %in% thetaNames,]

      for(thetaName in names(flatStandardizers)){
        adaptiveLassoWeights[thetaName] <- flatStandardizers[thetaName]*adaptiveLassoWeights[thetaName]
      }

      return(adaptiveLassoWeights)
    }

    if(tolower(penalty) == "adaptivelasso" && adaptiveLassoWeights == "auto"){
      stop("Automatic standardization and adaptive lasso weights can not be combined.")
      thetaNames <- names(OpenMx::omxGetParameters(mxObject))
      adaptiveLassoWeights <- abs(OpenMx::omxGetParameters(mxObject))^(-1)
      # compute standardizers
      T0VAR <- getT0VAR(mxObject)
      driftLabels <- mxObject$DRIFT$labels
      if(anyNA(driftLabels)){
        autoDriftLabels <- matrix(paste0("_autoDriftLabel_", rep(seq_len(nrow(driftLabels)), each = ncol(driftLabels)), "_", seq_len(ncol(driftLabels))),
                                  nrow = nrow(driftLabels), ncol = ncol(driftLabels), byrow = T)
        driftLabels[is.na(driftLabels)] <- autoDriftLabels[is.na(driftLabels)]
      }

      flatStandardizers <- getFlatStdizer(T0VAR = T0VAR, driftLabels = driftLabels)
      flatStandardizers <- flatStandardizers[rownames(flatStandardizers) %in% thetaNames,]

      for(thetaName in rownames(flatStandardizers)){
        adaptiveLassoWeights[thetaName] <- flatStandardizers[thetaName]*adaptiveLassoWeights[thetaName]
      }

      return(adaptiveLassoWeights)
    }
  }

  ## if lasso was requested, but no automatic standardization
  if(tolower(penalty) == "lasso"){
    thetaNames <- names(OpenMx::omxGetParameters(mxObject))
    adaptiveLassoWeights <- rep(1,length(thetaNames))
    names(adaptiveLassoWeights) <- thetaNames
    return(adaptiveLassoWeights)
  }

  if(tolower(penalty) == "ridge"){
    thetaNames <- names(OpenMx::omxGetParameters(mxObject))
    adaptiveLassoWeights <- rep(1,length(thetaNames))
    names(adaptiveLassoWeights) <- thetaNames
    return(adaptiveLassoWeights)
  }

  ## if adatpive lasso with automatic weigths was requested, but no automatic standardization
  if(tolower(penalty) == "adaptivelasso" && adaptiveLassoWeights == "auto"){
    adaptiveLassoWeights <- abs(OpenMx::omxGetParameters(mxObject))^(-1)
    return(adaptiveLassoWeights)
  }

  # else
  stop("Error while computing adaptive lasso weigths")

}




#' getFinalParameters
#'
#' returns the final parameters for a regularized model
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param regCtsemObject fitted regularized continuous time model
#' @param criterion select a criterion. Possible are AIC, BIC, cvM2LL
#' @param raw boolean: should the raw parameters be returned (raw = TRUE) or should the log-transformed variances be re-transformed (raw = FALSE)
#' @author Jannik Orzek
#' @import OpenMx
#' @export
getFinalParameters <- function(regCtsemObject, criterion = NULL, raw = FALSE){
  if(!regCtsemObject$setup$autoCV){
    minCriterionValue <- max(which(regCtsemObject$fit[criterion,] == min(regCtsemObject$fit[criterion,], na.rm = TRUE)))
    lambdas <- regCtsemObject$setup$lambdas
    bestLambda <- lambdas[minCriterionValue]
    if(raw){
      parameters <- regCtsemObject$parameterEstimatesRaw[,as.character(bestLambda)]
    }else{
      parameters <- regCtsemObject$parameters[,as.character(bestLambda)]
    }
    return(list("criterion" = criterion,
                "lambda" = bestLambda,
                "parameters" = parameters))
  }
  minCriterionValue <- max(which(regCtsemObject$fit["mean CV fit",] == min(regCtsemObject$fit["mean CV fit",], na.rm = TRUE)))
  lambdas <- regCtsemObject$setup$lambdas
  bestLambda <- lambdas[minCriterionValue]
  return(list("criterion" = "mean CV fit",
              "lambda" = bestLambda))
}

#' getFinalModel
#'
#' returns the final model for a regularized model. Note: Returns as mxObject!
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param regCtsemObject fitted regularized continuous time model
#' @param criterion select a criterion. Possible are AIC, BIC, cvM2LL
#' @author Jannik Orzek
#' @import OpenMx
#' @export
getFinalModel <- function(regCtsemObject, criterion = NULL){
  if(regCtsemObject$setup$autoCV){
    stop("getFinalModel not supported for automatic cross-validation. At the moment, you have to manually re-run the model with the best lambda value using the whole sample.")
  }
  bestPars <- getFinalParameters(regCtsemObject, criterion = criterion, raw = TRUE)
  message(paste0("Best fit for ", criterion, " was observed for lambda = ", bestPars$lambda, "."))
  finalModel <- OpenMx::omxSetParameters(regCtsemObject$setup$mxObject, values = bestPars$parameters, labels = names(bestPars$parameters))
  finalModel <- OpenMx::mxRun(finalModel, useOptimizer = FALSE, silent = TRUE)
  return(finalModel)
}


#' getFlatStdizer
#'
#' returns the standardizer for the standardized drift parameters if standardizeDrift = TRUE. Computes T0SD_predictor/T0SD_dependent
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param T0VAR matrix with T0VAR values
#' @param driftLabels vector with drift names
#' @export
getFlatStdizer <- function(T0VAR, driftLabels){
  stdizer <- matrix(1, nrow = nrow(T0VAR), ncol = ncol(T0VAR))
  stdizer <- stdizer%*%diag(sqrt(diag(T0VAR))) # times predictor sd
  stdizer <- t(t(stdizer)%*%diag(sqrt(diag(T0VAR))^(-1))) # divided by dependent sd
  flatStdizer <- OpenMx::cvectorize(stdizer) # flatten
  rownames(flatStdizer) <- as.vector(driftLabels)
  return(flatStdizer)
}




#' getMaxLambda
#'
#' computes an approximation of the lowest lambda which will set all regularized parameters to zero. This function is adapted from Murphy (2012) Machine learning: a probabilistic perspective. See p. 434 for more details.
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param regIndicators Labels of the regularized parameters (e.g. drift_eta1_eta2)
#' @param differenceApprox Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the unregularized parameter estimates.
#' @author Jannik Orzek
#' @import OpenMx
#' @export
getMaxLambda <- function(mxObject, regIndicators, differenceApprox, adaptiveLassoWeights){
  # This function is adapted from Murphy (2012) Machine learning: a probabilistic perspective. See p. 434 for more details.
  cat("Computing maximally required lambda ... ")

  converged <- FALSE
  # warning("automatically determining the maximal lambda with getLambdaMax is still experimental! It only produces an approximation of the required lambdaMax which can be too large or too small.")

  it <- 0

  while(!converged){
    if(it == length(regIndicators)){
      stop("Error while automatically setting the lambdas: The models did not converge. Try setting the lambdas manually.")
    }
    # extract parameter vector
    param <- OpenMx::omxGetParameters(mxObject)

    numberRegularized <- length(regIndicators) - it

    # Problem: setting all regularized parameters to zero might result in an impossible model
    # if this error occurs, regCtsem will iteratively try to set a subset of the parameters to zero
    # the size of the parameters determines the order in which they are set to zero
    # however, this is only a rough approximation and might result in an unsatisfactory maxLambda
    regIndicatorsCurrent <- names(sort(abs(param[regIndicators]))[1:numberRegularized])

    # step 1: set the regularized parameters to zero and estimate the model:

    param[regIndicatorsCurrent] <- 0
    freeParam <- rep(TRUE, length(param))
    names(freeParam) <- names(param)
    freeParam[regIndicatorsCurrent] <- FALSE

    sparseModel <- mxObject
    sparseModel <- OpenMx::omxSetParameters(model = sparseModel,
                                            labels = names(param),
                                            free = freeParam,
                                            values = param)
    sparseModel <- try(OpenMx::mxRun(sparseModel, silent = TRUE))
    if(any(class(sparseModel) == "try-error") | sparseModel$output$status$code == 10 | sparseModel$output$status$code == 5){
      cat("using TryHardctsem ...\n")
      sparseModel <- try(OpenMx::mxTryHardctsem(sparseModel, silent = TRUE))
    }
    if(any(class(sparseModel) == "try-error") | sparseModel$output$status$code == 10){
      if(it ==  0){
        warning("Error when determining the lambdas automatically: Setting all regularized parameters to zero resulted in an impossible model. regCtsem will try to at least set a subset of the regularized parameters to zero; however, this might result in a wrong maximum for the lambdas! Consider setting the lambdas manually.")
      }
      it <- it + 1
      next
    }

    nonZeroParam <- OpenMx::omxGetParameters(sparseModel)
    namesNonZeroParam <- names(nonZeroParam)

    # step 2: compute gradients with regularized parameters to zero and unregularized parameters set to nonZeroParam estimates
    param[namesNonZeroParam] <- nonZeroParam

    gradientModel <- OpenMx::mxModel(mxObject,
                                     OpenMx::mxComputeSequence(steps=list(OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
                                                                                                        hessian = FALSE))
                                     ))
    gradientModel <- OpenMx::omxSetParameters(model = gradientModel, labels = names(param), values = param)
    gradientModel <- try(OpenMx::mxRun(gradientModel, silent = TRUE))

    if(any(class(gradientModel) == "try-error")){
      if(it ==  0){
        warning("Error when determining the lambdas automatically: Setting all regularized parameters to zero resulted in an impossible model. regCtsem will try to at least set a subset of the regularized parameters to zero; however, this might result in a wrong maximum for the lambdas! Consider setting the lambdas manually.")
      }
      it <- it + 1
      next
    }

    # extract the gradient
    grad <- gradientModel$compute$steps[[1]]$output[["gradient"]]
    gradLabels <- rownames(grad)
    grad <- grad[,differenceApprox] # use specified gradient approximation
    names(grad) <- gradLabels

    converged <- TRUE
  }
  if(it > 0){
    warning(paste0("regCtsem did set ", numberRegularized, " of the ", length(regIndicators), " regularized parameters to zero when determining the maximal lambda."))
  }
  # define maxLambda as the maximal gradient of the regularized parameters
  maxLambda <- max(abs(grad[regIndicators]) * adaptiveLassoWeights[regIndicators]^(-1))

  cat("DONE \n")

  return(list("maxLambda" = maxLambda, "sparseParameters" = param))
}




#' getMaxLambdaCV
#'
#'
#' computes an approximation of the lowest lambda which will set all regularized parameters to zero. This function is adapted from Murphy (2012) Machine learning: a probabilistic perspective. See p. 434 for more details.
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param fullData Dataset for all samples combined
#' @param KalmanStartValues Optional starting values for the parameters when using Kalman filter
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso and adaptiveLasso for optimization = exact
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the unregularized parameter estimates.
#' @param standardizeDrift Boolean: Should Drift parameters be standardized automatically using T0VAR?
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @param cvFolds list of numeric vectors indicating which data row in fullData belongs to which cv sample
#' @param trainSets empty list of the same length as k. The data sets for training the models will be saved here
#' @param regIndicators Labels of the regularized parameters (e.g. drift_eta1_eta2)
#' @param differenceApprox Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the unregularized parameter estimates.
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param optimization which optimization procedure should be used. Possible are  "exact" or "approx".
#' @param differenceApprox Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#'
#' @author Jannik Orzek
#' @import OpenMx
#' @export

getMaxLambdaCV <- function(ctsemObject, mxObject, fullData, KalmanStartValues, regOn, regIndicators, penalty, adaptiveLassoWeights, standardizeDrift, k, cvFolds, trainSets, scaleLambdaWithN, objective, optimization, differenceApprox = "central"){

  maxLambda <- 0

  for(foldNumber in 1:k){
    cat(paste0("Sample ", foldNumber, ": "))
    if(is.vector(fullData[-cvFolds[[foldNumber]],])){
      # if a single row is selected, R extracts this row as vector, not as matrix
      trainSets[[foldNumber]] <-  t(as.matrix(fullData[-cvFolds[[foldNumber]],]))
    }else{
      trainSets[[foldNumber]] <- fullData[-cvFolds[[foldNumber]],]
    }
    sampleSize <- length(unlist((cvFolds))) - length(cvFolds[[foldNumber]])

    if(tolower(objective) == "ml"){
      currentModel <- mxObject
      currentModel$data <- OpenMx::mxData(trainSets[[foldNumber]], type = "raw")
      currentModel <- mxRun(currentModel, silent = TRUE)

      currentAdaptiveLassoWeights <- getAdaptiveLassoWeights(mxObject = currentModel, penalty = penalty, adaptiveLassoWeights = adaptiveLassoWeights, standardizeDrift = standardizeDrift)

    }else if(tolower(objective) == "kalman"){
      # set input arguments
      currentModel <- createKalmanMultiSubjectModel(ctsemObject = ctsemObject,
                                                    dataset = trainSets[[foldNumber]],
                                                    useOptimizer = TRUE,
                                                    KalmanStartValues = KalmanStartValues)

      currentAdaptiveLassoWeights <- getAdaptiveLassoWeights(mxObject = currentModel,
                                                             penalty = penalty,
                                                             adaptiveLassoWeights = adaptiveLassoWeights,
                                                             standardizeDrift = standardizeDrift)
    }
    if(tolower(optimization) == "approx"){
      regIndicatorsString <- mxObject[[regOn]]$labels[regIndicators == 1]
      currentMaxLambda <- getMaxLambda(mxObject = currentModel,
                                           regIndicators = regIndicatorsString,
                                           differenceApprox = differenceApprox,
                                           adaptiveLassoWeights = currentAdaptiveLassoWeights)$maxLambda
      sparseParameterMatrix <- NULL
    }else{
      parameterLabels <- names(OpenMx::omxGetParameters(currentModel))
      currentMaxLambda <- getMaxLambda(mxObject = currentModel,
                                           regIndicators = regIndicators,
                                           differenceApprox = differenceApprox,
                                           adaptiveLassoWeights = currentAdaptiveLassoWeights)
      currentSparseParameters <- currentMaxLambda$sparseParameters
      currentMaxLambda <- currentMaxLambda$maxLambda
      if(foldNumber == 1){
        sparseParameterMatrix <- matrix(NA, nrow = length(parameterLabels), ncol = k)
        rownames(sparseParameterMatrix) <- parameterLabels
        sparseParameterMatrix[parameterLabels,foldNumber] <- currentSparseParameters[parameterLabels]
      }else{
        sparseParameterMatrix[parameterLabels,foldNumber] <- currentSparseParameters[parameterLabels]
      }

    }

    if(scaleLambdaWithN){
      currentMaxLambda <- currentMaxLambda/sampleSize
    }

    maxLambda <- max(maxLambda,
                       currentMaxLambda)
  }

  return(list("maxLambda" = maxLambda, "sparseParameterMatrix" = sparseParameterMatrix))
}




#' getT0VAR
#'
#' extracts the T0VAR from an mxObject
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param mxObject mxObject
#' @export
getT0VAR <- function(mxObject){

  # case 1: T0VAR not stationary
  if(!is.null(mxObject$T0VAR$result)){
    return(mxObject$T0VAR$result)
  }
  # case 2: T0VAR stationary
  stop("stop from getT0VAR: Automatic standardization not yet implemented for stationarity T0VAR")

}




#' getVarianceFromVarianceBase2
#'
#' returns variance from varianceBase
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param varianceBaseValues values of varianceBase
#' @import OpenMx
#' @author Jannik Orzek
#' @export
getVarianceFromVarianceBase2 <- function(varianceBaseValues){
  varianceCholValues <- OpenMx::vec2diag(exp(OpenMx::diag2vec(varianceBaseValues))) +
    varianceBaseValues - OpenMx::vec2diag(OpenMx::diag2vec(varianceBaseValues))
  varianceValues <- varianceCholValues %*% t(varianceCholValues)
  return(varianceValues)
}




#' getVariancesInParameterEstimates
#'
#' computes the variances from the parameterEstimates matrix
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param mxObject mxObject
#' @param parameterEstimates parameter estimates from regularized ctsem
#' @export
getVariancesInParameterEstimates <- function(mxObject, parameterEstimates){
  if(!is.matrix(parameterEstimates)){
    stop("parameterEstimates has to be of class matrix")
  }

  latentNames <- diag(mxObject$DRIFT$labels)
  latentNames <- sub(x = latentNames, pattern = "drift_", replacement = "")
  manifestNames <- rownames(mxObject$LAMBDA$values)

  for(lambda in 1:ncol(parameterEstimates)){
    tempMxObject <- mxObject
    tempMxObject <- OpenMx::omxSetParameters(model = tempMxObject,
                                             labels = rownames(parameterEstimates),
                                             values = parameterEstimates[,lambda]
    )

    # T0VAR
    if(any(mxObject$T0VARbase$free)){
      T0VAR <- regCtsem::getVarianceFromVarianceBase2(varianceBaseValues = tempMxObject$T0VARbase$values)

      if(lambda == 1){
        T0VARBaseLabels <- tempMxObject$T0VARbase$labels[!is.na(tempMxObject$T0VARbase$labels)]
        T0VARLabels <- paste0("T0VAR_",rep(latentNames, each = length(latentNames)),
                              "_",
                              rep(latentNames, length(latentNames)))
        T0VARLabels <- matrix(T0VARLabels,
                              nrow = length(latentNames),
                              ncol = length(latentNames),
                              byrow = TRUE)

        T0VARs <- OpenMx::cvectorize(T0VAR)
      }else{
        T0VARs <- cbind(T0VARs, OpenMx::cvectorize(T0VAR))
      }

    }

    # DIFFUSION

    if(any(mxObject$DIFFUSIONbase$free)){
      DIFFUSION <- regCtsem::getVarianceFromVarianceBase2(varianceBaseValues = tempMxObject$DIFFUSIONbase$values)

      if(lambda == 1){
        DIFFUSIONBaseLabels <- tempMxObject$DIFFUSIONbase$labels[!is.na(tempMxObject$DIFFUSIONbase$labels)]
        DIFFUSIONLabels <- paste0("DIFFUSION_",rep(latentNames, each = length(latentNames)),
                                  "_",
                                  rep(latentNames, length(latentNames)))
        DIFFUSIONLabels <- matrix(DIFFUSIONLabels,
                                  nrow = length(latentNames),
                                  ncol = length(latentNames),
                                  byrow = TRUE)

        DIFFUSIONs <- OpenMx::cvectorize(DIFFUSION)
      }else{
        DIFFUSIONs <- cbind(DIFFUSIONs, OpenMx::cvectorize(DIFFUSION))
      }

    }

    # MANIFESTVAR

    if(any(mxObject$MANIFESTVARbase$free)){
      MANIFESTVAR <- regCtsem::getVarianceFromVarianceBase2(varianceBaseValues = tempMxObject$MANIFESTVARbase$values)

      if(lambda == 1){
        MANIFESTVARBaseLabels <- tempMxObject$MANIFESTVARbase$labels[!is.na(tempMxObject$MANIFESTVARbase$labels)]
        MANIFESTVARLabels <- paste0("MANIFESTVAR_",rep(manifestNames, each = length(manifestNames)),
                                    "_",
                                    rep(manifestNames, length(manifestNames)))
        MANIFESTVARLabels <- matrix(MANIFESTVARLabels,
                                    nrow = length(manifestNames),
                                    ncol = length(manifestNames),
                                    byrow = TRUE)

        MANIFESTVARs <- OpenMx::cvectorize(MANIFESTVAR)
      }else{
        MANIFESTVARs <- cbind(MANIFESTVARs, OpenMx::cvectorize(MANIFESTVAR))
      }

    }

    if(any(grepl("asymDIFFUSIONalg", mxObject$T0VAR$labels))){
      # in this case: T0VAR is set to stationarity
      DRIFTHATCH <- tempMxObject$DRIFT$values %x% tempMxObject$II$values + tempMxObject$II$values %x% tempMxObject$DRIFT$values
      asymDIFFUSIONalg <- -solve(DRIFTHATCH) %*% cvectorize(DIFFUSION)
      T0VAR <- matrix(asymDIFFUSIONalg, nrow = length(latentNames), byrow = F)

      if(lambda == 1){
        T0VARBaseLabels <- c()
        T0VARLabels <- paste0("T0VAR_",rep(latentNames, each = length(latentNames)),
                              "_",
                              rep(latentNames, length(latentNames)))
        T0VARLabels <- matrix(T0VARLabels,
                              nrow = length(latentNames),
                              ncol = length(latentNames),
                              byrow = TRUE)

        T0VARs <- OpenMx::cvectorize(T0VAR)
      }else{
        T0VARs <- cbind(T0VARs, OpenMx::cvectorize(T0VAR))
      }

    }

  }
  # replace values in parameterEstimates
  # T0VAR
  if(any(mxObject$T0VARbase$free)){
    rownames(T0VARs) <- OpenMx::cvectorize(T0VARLabels)
    parameterEstimates <- parameterEstimates[!(rownames(parameterEstimates) %in% T0VARBaseLabels),]
    if(!is.matrix(parameterEstimates)){
      parLabels <- names(parameterEstimates)
      parameterEstimates <- matrix(parameterEstimates, nrow = length(parLabels))
      rownames(parameterEstimates) <- parLabels
    }
    parameterEstimates <- rbind(parameterEstimates, T0VARs)
  }
  if(any(grepl("asymDIFFUSIONalg", mxObject$T0VAR$labels))){
    rownames(T0VARs) <- OpenMx::cvectorize(T0VARLabels)
    if(!is.matrix(parameterEstimates)){
      parLabels <- names(parameterEstimates)
      parameterEstimates <- matrix(parameterEstimates, nrow = length(parLabels))
      rownames(parameterEstimates) <- parLabels
    }
    parameterEstimates <- rbind(parameterEstimates, T0VARs)
  }

  # DIFFUSION
  if(any(mxObject$DIFFUSIONbase$free)){
    rownames(DIFFUSIONs) <- OpenMx::cvectorize(DIFFUSIONLabels)
    parameterEstimates <- parameterEstimates[!(rownames(parameterEstimates) %in% DIFFUSIONBaseLabels),]
    if(!is.matrix(parameterEstimates)){
      parLabels <- names(parameterEstimates)
      parameterEstimates <- matrix(parameterEstimates, nrow = length(parLabels))
      rownames(parameterEstimates) <- parLabels
    }
    parameterEstimates <- rbind(parameterEstimates, DIFFUSIONs)
  }

  # MANIFESTVAR
  if(any(mxObject$MANIFESTVARbase$free)){
    rownames(MANIFESTVARs) <- OpenMx::cvectorize(MANIFESTVARLabels)
    parameterEstimates <- parameterEstimates[!(rownames(parameterEstimates) %in% MANIFESTVARBaseLabels),]
    if(!is.matrix(parameterEstimates)){
      parLabels <- names(parameterEstimates)
      parameterEstimates <- matrix(parameterEstimates, nrow = length(parLabels))
      rownames(parameterEstimates) <- parLabels
    }
    parameterEstimates <- rbind(parameterEstimates, MANIFESTVARs)
  }

  return(parameterEstimates)

}




#' regIndicatorsFromNameToMatrix
#'
#' transforms string indicated regIndicators to matrix indicated regIndicators
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param regOn string specifying which matrix should be regularized. Currently only supports DRIFT
#' @param regIndicators matrix with ones and zeros specifying which parameters in regOn should be regularized. Must be of same size as the regularized matrix. 1 = regularized, 0 = not regularized. Alternatively, labels for the regularized parameters can be used (e.g. drift_eta1_eta2)
regIndicatorsFromNameToMatrix <- function(mxObject, regOn, regIndicators){
  regIndicatorsMatrix <- matrix(as.numeric(c(mxObject[[regOn]]$labels %in% regIndicators)),
                                ncol = ncol(mxObject[[regOn]]$labels),
                                nrow = nrow(mxObject[[regOn]]$labels))
  return(regIndicatorsMatrix)
}




#' separateFitAndParameters
#'
#' separates fit results from parameter estimates
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param regCtsemObject Object from exact_regCtsem or approx_regCtsem
#' @import OpenMx
#' @author Jannik Orzek
#' @export
separateFitAndParameters <- function(regCtsemObject){
  # get parameter labels:

  parameterLabels <- names(OpenMx::omxGetParameters(regCtsemObject$setup$mxObject))
  fitAndParametersLabels <- rownames(regCtsemObject$fitAndParameters)
  fitLabels <- fitAndParametersLabels[!(fitAndParametersLabels %in% parameterLabels)]

  fit <- regCtsemObject$fitAndParameters[fitLabels,]
  parameterEstimates <- regCtsemObject$fitAndParameters[parameterLabels,]

  if(!is.matrix(parameterEstimates)){
    parameterEstimates <- matrix(parameterEstimates, nrow = length(parameterLabels))
    rownames(parameterEstimates) <- parameterLabels
  }
  if(!is.matrix(fit)){
    fit <- matrix(fit, nrow = length(fitAndParametersLabels[fitLabels]))
    rownames(fit) <- fitLabels
  }
  return(list("fit" = fit, "parameterEstimates" = parameterEstimates))

}




#' generateDRIFTPlot
#'
#' generates plot of drift values
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param model model from regCtsem
#' @param ylab label for y-axis. auto will set to lambda
#' @param xlab label for x-axis. auto will set to drift
#' @param skiptYlabComp boolean: is ylab a string
#' @param skiptXlabComp boolean: is xlab a string
#' @export
generateDRIFTPlot <- function(model, ylab = "auto", xlab = "auto", skiptYlabComp, skiptXlabComp){
  drifts <- model$parameters[grepl("drift", rownames(model$parameters)),]
  regIndicators <- model$setup$regIndicators
  drifts_regularized <- drifts[regIndicators,]
  lambdas <- model$setup$lambdas

  colsWithNA <- apply(model$parameters,2,function(x) any(is.na(x)))

  color <- ifelse(rownames(drifts) %in% regIndicators, yes = "white", "black")

  par(mar=c(4, 3, 5, 2), xpd=TRUE)

  matplot(x = lambdas[!colsWithNA], t(drifts[!colsWithNA]), lty = 2,
          lwd = 2, type = "l",
          col = color, xlab = ifelse(skiptXlabComp, xlab,
                                     ifelse(xlab == "auto",expression(lambda),ylab)),
          ylab = ifelse(skiptYlabComp, ylab,
                        ifelse(ylab == "auto","drift",ylab)))
  matplot(x = lambdas[!colsWithNA], t(drifts_regularized[!colsWithNA]), lty = 1, lwd = 2, type = "l", col = "#2166AC", add = TRUE)
  tickat <- seq(1, length(lambdas[!colsWithNA]), length.out = 10)
  axis(3, at = lambdas[tickat],
       labels=apply(drifts == 0,2,sum)[tickat],
       outer= F,
       line=1,col="black",col.ticks="black",col.axis="black")
  mtext("# zeroed parameters",3,line=3,at=mean(lambdas),col="black", cex = 1)

  par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
}

#' %dopar%
#' define operator dopar
#' @export
`%dopar%` <- foreach::`%dopar%`



