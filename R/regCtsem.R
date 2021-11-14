#' regCtsem
#'
#' main function: performs regularized ctsem
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param ctsemObject Fitted object of class ctsemFit
#' @param dataset Data set in wide format compatible with ctsemOMX
#' @param regIndicators Labels of the regularized parameters (e.g. drift_eta1_eta2).
#' @param targetVector named vector with values towards which the parameters are regularized (Standard is regularization towards zero)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01). Alternatively, lambdas can be set to "auto". regCtsem will then compute an upper limit for lambda and test lambdasAutoLength increasing lambda values
#' @param lambdasAutoLength if lambdas == "auto", lambdasAutoLength will determine the number of lambdas tested.
#' @param penalty Currently supported are lasso, ridge and adaptiveLasso
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the inverse of unregularized parameter estimates.
#' @param adaptiveLassoPower power for the adaptive lasso weights. The weights will be set to parameterValues^{adaptiveLassoPower}, where parameterValues refers to the unregularized maximum likelihood estimates
#' @param cvSample cross-validation sample. Has to be in wide format and compatible with ctsemOMX
#' @param autoCV Should automatic cross-validation be used? Possible are "No", "kFold" or "Blocked". kFold splits the dataset in k groups by selecting independent units from the rows. Blocked is a within-unit split, where for each person blocks of observations are deleted. See Bulteel, K., Mestdagh, M., Tuerlinckx, F., & Ceulemans, E. (2018). VAR(1) based models do not always outpredict AR(1) models in typical psychological applications. Psychological Methods, 23(4), 740–756. https://doi.org/10.1037/met0000178 for a more detailed explanation
#' @param k number of cross-validation folds if autoCV = "kFold" or autoCV = "Blocked"
#' @param sparseParameters labeled vector with parameter estimates of the most sparse model. If regValues = "auto" the sparse parameters will be computed automatically.
#' @param subjectSpecificParameters EXPERIMENTAL! A vector of parameter labels for parameters which should be estimated person-specific. If these parameter labels are also passed to regIndicators, all person-specific parameters will be regularized towards a group-parameter. This is a 2-step-procedure: In step 1 all parameters are constrained to equality between individuals to estimate the group parameters. In step 2 the parameters are estimated person-specific, but regularized towards the group parameter from step 1.
#' @param standardizeDrift Should Drift parameters be standardized automatically? Set to 'No' for no standardization, 'T0VAR' for standardization using the T0VAR or 'asymptoticDiffusion' for standardization using the asymptotic diffusion
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended as the likelihood is also sample size dependent
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param BICWithNAndT Boolean: TRUE = Use N and T in the formula for the BIC (-2log L + log(N+T)*k, where k is the number of parameters in the model). FALSE = Use N in the formula for the BIC (-2log L + log(N)). Defaults to FALSE
#' @param optimization which optimization procedure should be used. Possible are  "exact" or "approx". exact is recommended for sparsity inducing penalty functions (lasso and adaptive lasso)
#' @param optimizer for exact optimization: Either GIST or GLMNET. When using optimization = "approx", Rsolnp or any of the optimizers in optimx can be used. See ?optimx
#' @param control List with control arguments for the optimizer. See ?controlGIST, ?controlGLMNET and ?controlApprox for the respective parameters
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress.
#' @param trainingWheels If set to FALSE all bells and whistles used to keep regCtsem on track are turned off (no multiple starting values, no initial optimization with solnp or optimx). The focus is speed instead of accuracy. This might work in simulated data, but is NOT recommended with real data. The optimizer is quite likely to get stuck in local minima.
#' @return returns an object of class regCtsem. Without cross-validation, this object will have the fields setup (all arguments passed to the function), fitAndParameters (used internally to store the fit and the raw (i.e., untransformed) parameters), fit (fit indices, ect.), parameterEstimatesRaw (raw, i.e. untransformed parameters; used internally), and parameters (transformed parameters)
#' @examples
#' \donttest{
#'
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
#' regIndicators <- c("drift_eta2_eta1", "drift_eta1_eta2")
#' # Note: If you are unsure what the parameters are called in
#' # your model, check: fit_myModel$ctmodelobj$DRIFT for the drift or
#' # omxGetParameters(fit_myModel$ctmodelobj) for all parameters
#'
#' # Optimize model using GIST with lasso penalty
#' regModel <- regCtsem::regCtsem(ctsemObject = fit_myModel,
#'                                dataset = dat,
#'                                regIndicators = regIndicators,
#'                                lambdas = "auto",
#'                                lambdasAutoLength = 20)
#' summary(regModel)
#' plot(regModel)
#' plot(regModel, what = "fit", criterion = c("AIC", "BIC", "m2LL"))
#'
#' # The best parameter estimates and the final model as mxObject can be extracted with:
#' # getFinalParameters(regCtsemObject = regModel, criterion = "BIC")
#' # bestModel <- getFinalModel(regCtsemObject = regModel, criterion = "BIC")
#' # WARNING: The returned model is of type cpptsem. You can access it's elements with the
#' # $ operator. For example: bestModel$DRIFTValues
#'
#' # WARNING: If you load an existing regCtsem object, the underlying C++ model will no longer
#' # exist. You can restore this model with restore(). Example:
#' # save(regModel, file = "regModel.RData")
#' # load("regModel.RData")
#' # regModel <- restore(regModel)
#' # regModel$setup$cpptsemObject$DRIFTValues
#'
#' # Optimize model using GLMNET with lasso penalty
#' regModel <- regCtsem::regCtsem(ctsemObject = fit_myModel,
#'                                dataset = dat,
#'                                regIndicators = regIndicators,
#'                                lambdas = "auto",
#'                                lambdasAutoLength = 20,
#'                                optimizer = "GLMNET")
#'
#' summary(regModel, criterion = "BIC")
#' plot(regModel, what = "drift")
#' plot(regModel, what = "fit", criterion = c("AIC", "BIC", "m2LL"))
#' plot(regModel, what = "drift_eta1_eta2")
#'
#' # The same regularization can be performed with the approximate optimization:
#' regModelApprox <- regCtsem::regCtsem(ctsemObject = fit_myModel,
#'                                      dataset = dat,
#'                                      regIndicators = regIndicators,
#'                                      lambdas = "auto",
#'                                      lambdasAutoLength = 20,
#'                                      optimization = "approx",
#'                                      control = list(
#'                                        epsilon = .001, # epsilon is used to transform the non-differentiable
#'                                        #lasso penalty to a differentiable one if optimization = approx
#'                                        zeroThresh = .04 # threshold below which parameters will be evaluated as == 0
#'                                      ))
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
#' fit_myModel <- ctFit(dat = traindata, ctmodelobj = myModel, objective = "Kalman")
#'
#' # select DRIFT values:
#' regIndicators <- c("drift_eta2_eta1", "drift_eta1_eta2")
#' # Note: If you are unsure what the parameters are called in
#' # your model, check: fit_myModel$ctmodelobj$DRIFT for the drift or
#' # omxGetParameters(fit_myModel$ctmodelobj) for all parameters
#'
#' ## Optimization with GIST:
#' regModel <- regCtsem::regCtsem(ctsemObject = fit_myModel,
#'                                dataset = traindata,
#'                                regIndicators = regIndicators,
#'                                lambdas = "auto",
#'                                lambdasAutoLength = 20,
#'                                cvSample = testdata # data set for cross-validation
#' )
#'
#' summary(regModel, criterion = "cvM2LL")
#' plot(regModel, what = "fit", criterion = "cvM2LL")
#'
#' #### EXPERIMENTAL FEATURES: USE WITH CAUTION! ####
#'
#' library(regCtsem)
#'
#' ## Example 4: Kalman Filter with person specific parameter values
#' ## WARNING: THIS WILL TAKE A WHILE TO RUN
#' set.seed(175446)
#'
#' ## define the population model:
#'
#' # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
#' ct_drift <- matrix(c(-.3,0,.2,-.2),2,2,TRUE)
#' dataset <- c()
#' indpars <- c()
#'
#' # We will simulate data for 10 individuals with person-specific parameters
#' # These person-specific parameters will then be regularized towards a
#' # common group parameter
#' for(i in 1:10){
#'   while(TRUE){
#'     DRIFT <- ct_drift + matrix(c(0,rnorm(1,0,.5),0,0),2,2,TRUE)
#'     if(!any(Re(eigen(DRIFT)$values) > 0)){break}
#'   }
#'   indpars <- c(indpars, DRIFT[1,2])
#'   generatingModel<-ctsem::ctModel(Tpoints=50,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                                   MANIFESTVAR=diag(0,2),
#'                                   LAMBDA=diag(1,2),
#'                                   DRIFT=DRIFT,
#'                                   DIFFUSION=matrix(c(.5,0,0,.5),2),
#'                                   CINT=matrix(0,nrow = 2, ncol = 1),
#'                                   T0MEANS=matrix(0,ncol=1,nrow=2),
#'                                   T0VAR=diag(1,2), type = "omx")
#'   dataset <- rbind(dataset, ctsem::ctGenerate(generatingModel,n.subjects = 1, wide = TRUE))
#'
#' }
#'
#' ## Build the analysis model.
#' myModel <- ctsem::ctModel(Tpoints=50,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
#'                           LAMBDA=diag(1,2),
#'                           MANIFESTVAR=diag(0,2),
#'                           DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
#'                           T0MEANS=matrix(0,ncol=1,nrow=2),
#'                           T0VAR="auto", type = "omx")
#' myModel <- ctFit(myModel, dat = dataset, objective = "Kalman",
#'                  useOptimizer = TRUE)
#' regIndicators <- c("drift_eta2_eta1", "drift_eta1_eta2")
#' # the following parameters will be estimated person-specific and (as we specified this above)
#' # regularized. The regularization will be towards a group parameter
#' subjectSpecificParameters <- c("drift_eta2_eta1", "drift_eta1_eta2")
#' regModel <- regCtsem(ctsemObject = myModel,
#'                      dataset = dataset,
#'                      regIndicators = regIndicators,
#'                      lambdasAutoLength = 5, # 5 will not be enough, but this takes some time to execute
#'                      subjectSpecificParameters = subjectSpecificParameters
#' )
#' summary(regModel, criterion = "BIC")
#' plot(regModel, what = "drift")
#' }
#'
#' @author Jannik Orzek
#' @import ctsemOMX rlist optimx Rsolnp
#' @export
regCtsem <- function(
  # model
  ctsemObject,  dataset,
  # penalty settings
  regIndicators,
  targetVector = NULL,
  lambdas = "auto",
  lambdasAutoLength = 50,
  penalty = "lasso",
  adaptiveLassoWeights = NULL,
  adaptiveLassoPower = -1,
  cvSample = NULL,
  autoCV = "No",
  k = 5,
  sparseParameters = NULL,
  subjectSpecificParameters = NULL,

  standardizeDrift = "No",
  scaleLambdaWithN = TRUE,
  returnFitIndices = TRUE,
  BICWithNAndT = FALSE,

  optimization = "exact",
  optimizer = "GIST",
  control = list(),
  # additional settings
  verbose = 0,
  trainingWheels = TRUE
){
  if(is.null(dataset)){stop("Data set in wide format required!")}

  # save input in list
  argsIn <- as.list(environment())

  # set objective
  argsIn$objective <- ifelse(ctsemObject$ctfitargs$objective == "Kalman", "Kalman", "ML")
  cat(paste0("Object with objective = ", ifelse(ctsemObject$ctfitargs$objective == "Kalman", "Kalman", "ML"), " detected.\n"))

  # save number of time points
  argsIn$Tpoints <- ctsemObject$ctmodelobj$Tpoints

  if(!is.null(subjectSpecificParameters)){
    if(argsIn$objective != "Kalman"){
      stop("Subject-specific parameters are only supported for models fitted with Kalman filter.")
    }
    if(!is.null(cvSample) || (autoCV == "kFold")){
      stop("kFold Cross-validation not supported with subject-specific parameters. Use 'No' or 'Blocked' instead.")
    }
    if(tolower(argsIn$optimizer) != "gist"){
      stop("Only GIST optimizer allows for person-specific parameters")
    }
    warning("Subject-specific parameters are a very experimental feature. Usage not recommended!")
  }

  ## Defaults for optimizer
  if(optimization == "exact" && penalty == "ridge"){
    optimization <- "approx"
    argsIn$optimization <- optimization
  }
  if(penalty == "ridge"){
    returnFitIndices <- FALSE
    argsIn$returnFitIndices <- returnFitIndices
    if(any(lambdas == "auto")){
      stop("Ridge regularization requires explicit specification of lambdas. lambdas = 'auto' is not allowed.")
    }
  }

  if(!trainingWheels){
    if(optimization == "approx"){
      controlTemp <- controlApprox(forceCpptsem = TRUE, nMultistart = 0)
      if(tolower(optimizer) == "gist" || tolower(optimizer) == "glmnet"){
        if(controlTemp$controlApproxOptimizer$package == "optimx"){
          optimizer <- controlTemp$controlApproxOptimizer$method
        }else{
          optimizer <- "solnp"
        }
      }
    } else if(optimization == "exact" && optimizer == "GIST"){
      controlTemp <- controlGIST(forceCpptsem = TRUE, approxFirst = FALSE, numStart = 0, nMultistart = 0)
    }else if(optimization == "exact" && optimizer == "GLMNET"){
      controlTemp <- controlGLMNET(tryCpptsem = TRUE, approxFirst = FALSE, numStart = 0, nMultistart = 0)
    }else{
      stop("Unknown optimization or optimizer specified. Possible are exact or approx for optimization and GIST or GLMNET for optimizer")
    }
  }else{
    if(optimization == "approx"){
      controlTemp <- controlApprox()
      if(tolower(optimizer) == "gist" || tolower(optimizer) == "glmnet"){
        if(controlTemp$controlApproxOptimizer$package == "optimx"){
          optimizer <- controlTemp$controlApproxOptimizer$method
        }else{
          optimizer <- "solnp"
        }
      }
    } else if(optimization == "exact" && optimizer == "GIST"){
      controlTemp <- controlGIST()
    }else if(optimization == "exact" && optimizer == "GLMNET"){
      controlTemp <- controlGLMNET()
    }else{
      stop("Unknown optimization or optimizer specified. Possible are exact or approx for optimization and GIST or GLMNET for optimizer")
    }
  }
  controlExpectedNames <- names(controlTemp)
  controlReceivedNames <- names(control)
  if(!all(controlReceivedNames %in% controlExpectedNames)){
    stop(paste0("Unknown parameters passed to control. Expected any of the following arguments: ", paste0(controlExpectedNames, collapse = ", ")))
  }
  controlTemp[controlReceivedNames] <- control[controlReceivedNames]

  argsIn <- c(argsIn, controlTemp)
  if(optimization == "exact" && optimizer == "GLMNET" && !is.matrix(argsIn$initialHessianApproximation)){
    if(argsIn$initialHessianApproximation == "OpenMx"){
      if(!is.null(ctsemObject$mxobj$output$hessian)){
        argsIn$initialHessianApproximation <- ctsemObject$mxobj$output$hessian
      }else{
        warning("Extracting Hessian from OpenMx failed. Using identity")
        argsIn$initialHessianApproximation <- "identity"
      }
    }

    if(!is.matrix(argsIn$initialHessianApproximation) && tolower(argsIn$initialHessianApproximation) == "identity"){
        argsIn$initialHessianApproximation <- diag(length(omxGetParameters(ctsemObject$mxobj)))
    }
  }

  if(!is.null(targetVector) && any(targetVector!=0)){
    if(optimizer == "GLMNET"){
      warning("Optimzer = GLMNET does not support a target vector with non-zero target values. Switching to GIST")
      optimizer <- "GIST"
      argsIn$optimizer <- optimizer
    }
  }
  if(is.null(targetVector)){
    targetVector <- rep(0, length(regIndicators))
    names(targetVector) <- regIndicators
    argsIn$targetVector <- targetVector
  }
  # check setup
  #checkSetup(argsIn)
  if(!is.null(cvSample) && (autoCV != "No")){stop("Either provide a cvSample or use autoCV. Combinations of both are currently not supported.")}
  if(autoCV == "kFold"){
    message("k-Fold cross-validation requested. Please note that between-person standardizing of data (mean centering or scaling) prior to the analysis can invalidate the results of cross-validation as it will result in dependencies between the training and validation set.")
  }
  if(autoCV == "Blocked"){
    message("Blocked cross-validation requested. Please note that within-person standardizing of data (mean centering or scaling) prior to the analysis can invalidate the results of cross-validation as it will result in dependencies between the training and validation set.")
  }

  # check if all targets are in the regIndicators
  if(!is.null(targetVector)){
    if((!all(names(targetVector) %in% regIndicators)) ||
       (!all(regIndicators %in% names(targetVector)))
    )
      stop("Names of targetVector and regIndicators do not match.")
  }

  # translate model to C++
  if(tolower(argsIn$objective)  == "ml"){

    cpptsemObject <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = dataset))
    argsIn$cpptsemObject <- cpptsemObject

    # if there is a cvSample: generate model for cvSample as well
    if(!is.null(cvSample)){
      cvSampleCpptsemObject <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = cvSample, silent = TRUE))
      argsIn$cvSampleCpptsemObject <- cvSampleCpptsemObject
      if(any(class(cvSampleCpptsemObject) == "try-error")){
        stop("Setting up the CV sample model failed.")
      }
    }else{
      argsIn$cvSampleCpptsemObject <- NULL
    }
    if(any(class(cpptsemObject) == "try-error")){
      stop("Setting up the cpptsem model failed. Try cpptsemFromCtsem function on your ctsemObject directly.")
    }else{
      # check fit
      cpptsemObject$computeRAM()
      cpptsemObject$fitRAM()
      m2LLcpp <- cpptsemObject$m2LL
      testM2LL <- round(ctsemObject$mxobj$fitfunction$result[[1]] - m2LLcpp,3) == 0
      if (!testM2LL & !argsIn$forceCpptsem){
        stop(paste0("Differences in fit between ctsem and cpptsem object: ", ctsemObject$mxobj$fitfunction$result[[1]], " vs ", m2LLcpp, ". Did you pass the same data to the ctsemObject and as dataset? This can also be a result of different implementations of the matrix exponential in the Eigen mathematical library used by OpenMx and the Armadillo library used by regCtsem. See vignette(topic = 'MatrixExponential', package = 'regCtsem') for more details. Use control = list('forceCpptsem' = TRUE) to ignore this warning."))
      }
    }
  }else if (tolower(argsIn$objective)  == "kalman"){
    if(!is.null(subjectSpecificParameters)){
      cpptsemObject <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = dataset, group = 1:nrow(dataset), groupSpecificParameters = subjectSpecificParameters))
    }else{
      cpptsemObject <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = dataset))
      argsIn$cpptsemObject <- cpptsemObject
    }
    # if there is a cvSample: generate model for cvSample as well
    if(!is.null(cvSample)){
      cvSampleCpptsemObject <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = cvSample, silent = TRUE))
      argsIn$cvSampleCpptsemObject <- cvSampleCpptsemObject
      if(any(class(cvSampleCpptsemObject) == "try-error")){
        stop("Setting up the CV sample model failed.")
      }
    }else{
      argsIn$cvSampleCpptsemObject <- NULL
    }
    if (any(class(cpptsemObject) == "try-error")){
      stop("Setting up the cpptsem model failed.")
    }else{
      # test fit
      cpptsemObject$computeAndFitKalman()
      m2LLcpp <- cpptsemObject$m2LL
      testM2LL <- round(ctsemObject$mxobj$fitfunction$result[[1]] - m2LLcpp,3) == 0
      if (!testM2LL & !argsIn$forceCpptsem){
        stop(paste0("Differences in fit between ctsem and cpptsem object: ", ctsemObject$mxobj$fitfunction$result[[1]], " vs ", m2LLcpp, ". Did you pass the same data to the ctsemObject and as dataset? This can also be a result of different implementations of the matrix exponential in the Eigen mathematical library used by OpenMx and the Armadillo library used by regCtsem. See vignette(topic = 'MatrixExponential', package = 'regCtsem') for more details. Use control = list('forceCpptsem' = TRUE) to ignore this warning."))
      }
    }
  }else{
    stop("Unknown objective")
  }

  if(any(c("DIFFUSIONbase", "T0VARbase", "MANIFESTVARbase") %in% unique(cpptsemObject$parameterTable$matrix[cpptsemObject$parameterTable$label %in% regIndicators]))){
    warning("You seem to be regularizing covariance parameters. Please not that regCtsem uses the log-Cholesky implementation to estimate covariance matrices. See Pinheiro, J. C., & Bates, D. M. (1996). Unconstrained parametrizations for variance-covariance matrices. Statistics and Computing, 6(3), 289–296. https://doi.org/10.1007/BF00140873. Check that the regularization actually does what it should!")
  }

  ## Additional settings for person specific parameter estimates
  if(!is.null(subjectSpecificParameters)){
    # optimize subject specific model
    message("Fitting model with person-specific parameter estimates.")
    startingValues <- cpptsemObject$getParameterValues()

    # optimize
    CpptsemFit <- try(optimizeCpptsem(cpptsemObject = cpptsemObject, nMultistart = argsIn$nMultistart),
                      silent = TRUE)

    if(any(class(CpptsemFit) == "try-error")){stop("Rsolnp for the model with person-specific parameter estimates resulted in errors.")}
    if(CpptsemFit$convergence > 0){warning(paste0("Rsolnp reports convcode  > 0 for the model with person-specific parameter estimates: ", CpptsemFit$convcode, ". See ?optimx for more details."))}

    # compute unregularized fit
    cpptsemObject$setParameterValues(CpptsemFit$pars, names(CpptsemFit$pars))
    cpptsemObject$computeAndFitKalman()

    ## add person-specific parameter estimates to the regularized parameters if requested
    subjectSpecificParameterLabels <- c()
    for(grPar in subjectSpecificParameters){
      subjectSpecificParameterLabels <- paste0(grPar, "_G", 1:nrow(dataset))
      if(grPar %in% regIndicators){
        regIndicators <- regIndicators[-which(regIndicators == grPar)]
        regIndicators <- c(regIndicators, subjectSpecificParameterLabels)
        targetVector <- targetVector[-which(names(targetVector) == grPar)]
        # add group parameter value to targetValues
        targetVector <- c(targetVector, startingValues[subjectSpecificParameterLabels])
      }
    }
    argsIn$cpptsemObject <- cpptsemObject
    argsIn$regIndicators <- regIndicators
    argsIn$targetVector <- targetVector

  }

  # set adaptiveLassoWeights
  if(autoCV == "No"){
    argsIn$adaptiveLassoWeights <- getAdaptiveLassoWeights(cpptsemObject = cpptsemObject,
                                                           penalty = argsIn$penalty,
                                                           adaptiveLassoWeights = argsIn$adaptiveLassoWeights,
                                                           adaptiveLassoPower =  argsIn$adaptiveLassoPower,
                                                           standardizeDrift = argsIn$standardizeDrift)
    # if autoCV the adaptive lasso weights have to be set for each sub-sample (see below)
  }

  #### Exact Optimization ####
  #### without automatic cross-validation ####
  if(autoCV == "No" && argsIn$optimization == "exact"){
    regCtsemObject <- regCtsem::exact_regCtsem(cpptsemObject = cpptsemObject,
                                               dataset = argsIn$dataset,
                                               regIndicators = argsIn$regIndicators,
                                               targetVector = argsIn$targetVector,
                                               lambdas = argsIn$lambdas,
                                               lambdasAutoLength = argsIn$lambdasAutoLength,
                                               penalty = argsIn$penalty,
                                               adaptiveLassoWeights = argsIn$adaptiveLassoWeights,
                                               returnFitIndices = argsIn$returnFitIndices,
                                               BICWithNAndT = argsIn$BICWithNAndT,
                                               Tpoints = argsIn$Tpoints,
                                               cvSampleCpptsemObject = argsIn$cvSampleCpptsemObject,
                                               optimizer = argsIn$optimizer,
                                               objective = argsIn$objective,
                                               sparseParameters = argsIn$sparseParameters,
                                               stepSize = argsIn$stepSize,
                                               lineSearch = argsIn$lineSearch,
                                               c1 = argsIn$c1,
                                               c2 = argsIn$c2,
                                               sig = argsIn$sig,
                                               gam = argsIn$gam,
                                               initialHessianApproximation = argsIn$initialHessianApproximation,
                                               maxIter_out = argsIn$maxIter_out,
                                               maxIter_in = argsIn$maxIter_in,
                                               maxIter_line = argsIn$maxIter_line,
                                               eps_out = argsIn$eps_out,
                                               eps_in = argsIn$eps_in,
                                               eps_WW = argsIn$eps_WW,
                                               eta = argsIn$eta,
                                               stepsizeMin = argsIn$stepsizeMin,
                                               stepsizeMax = argsIn$stepsizeMax,
                                               GISTLinesearchCriterion = argsIn$GISTLinesearchCriterion,
                                               GISTNonMonotoneNBack = argsIn$GISTNonMonotoneNBack,
                                               break_outer = argsIn$break_outer,
                                               scaleLambdaWithN = argsIn$scaleLambdaWithN,
                                               approxFirst = argsIn$approxFirst,
                                               numStart = argsIn$numStart,
                                               nMultistart = argsIn$nMultistart,
                                               controlApproxOptimizer = argsIn$controlApproxOptimizer,
                                               verbose = argsIn$verbose)

    fitAndParametersSeparated <- try(separateFitAndParameters(regCtsemObject))
    if(!any(class(fitAndParametersSeparated) == "try-error")){
      regCtsemObject$fit <- fitAndParametersSeparated$fit
      regCtsemObject$parameterEstimatesRaw <- fitAndParametersSeparated$parameterEstimates
      regCtsemObject$parameters <- try(getParameterEstimates(regCtsemObject = regCtsemObject,
                                                             parameterEstimatesRaw = fitAndParametersSeparated$parameterEstimates),
                                       silent = TRUE)
      if(any(class(regCtsemObject$parameters) == "try-error")){
        warning("Could not compute the variances and covariances from the DIFFUSIONbase and T0VARbase. You will only find the raw paramter values in the output.")
      }
    }

    cat("\n")
    regCtsemObject$setup <- append(regCtsemObject$setup, argsIn[!(names(argsIn) %in% names(regCtsemObject$setup))])
    if(!is.null(subjectSpecificParameters)){
      class(regCtsemObject) <- "regCtsemMultiSubject"
    }else{
      class(regCtsemObject) <- "regCtsem"
    }
    return(regCtsemObject)
  }

  #### with automatic cross-validation ####
  # generate splits
  if((autoCV == "kFold" || autoCV == "Blocked") && argsIn$optimization == "exact"){
    cvFoldsAndModels <- regCtsem::createCVFoldsAndModels(dataset = argsIn$dataset,
                                                         Tpoints = argsIn$Tpoints,
                                                         manifestNames = ctsemObject$ctmodelobj$manifestNames,
                                                         k = argsIn$k, autoCV,
                                                         initialPars = ("T0MEANS" %in%  argsIn$cpptsemObject$parameterTable$matrix) || ("T0VARbase" %in%  argsIn$cpptsemObject$parameterTable$matrix))
    if(!is.numeric(lambdas) && lambdas == "auto"){
      maxLambdas <- matrix(NA, nrow = 1, ncol = argsIn$k)
      sparseParameters <- matrix(NA, nrow = length(cpptsemObject$getParameterValues()), ncol = argsIn$k)
      rownames(sparseParameters) <- names(cpptsemObject$getParameterValues())
    }

    # prepare adaptive lasso weights
    if(is.null(argsIn$adaptiveLassoWeights)){
      argsIn$adaptiveLassoWeights <- matrix(NA, nrow = k, ncol = length(cpptsemObject$getParameterValues()))
      rownames(argsIn$adaptiveLassoWeights) <- paste("CV", 1:k)
      colnames(argsIn$adaptiveLassoWeights) <- names(cpptsemObject$getParameterValues())
    }else{
      coln <- names(argsIn$adaptiveLassoWeights)
      argsIn$adaptiveLassoWeights <- matrix(rep(argsIn$adaptiveLassoWeights, k), nrow = k, ncol = length(argsIn$adaptiveLassoWeights), byrow = TRUE)
      colnames(argsIn$adaptiveLassoWeights) <- coln
    }

    for(i in 1:argsIn$k){

      if(!is.null(subjectSpecificParameters)){
        cvFoldsAndModels$trainModels[[i]] <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = cvFoldsAndModels$trainSets[[i]], group = 1:nrow(dataset), groupSpecificParameters = subjectSpecificParameters, silent = TRUE))
        cvFoldsAndModels$trainModels[[i]]$setParameterValues(cpptsemObject$getParameterValues(), names(cpptsemObject$getParameterValues()))
      }else{
        cvFoldsAndModels$trainModels[[i]] <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = cvFoldsAndModels$trainSets[[i]], silent = TRUE))
      }

      # optimize
      invisible(capture.output(fitTrain <- try(optimizeCpptsem(cpptsemObject = cvFoldsAndModels$trainModels[[i]], nMultistart = argsIn$nMultistart),
                                               silent = TRUE), type = c("output", "message")))

      if(!any(class(fitTrain) == "try-error")){
        cvFoldsAndModels$trainModels[[i]]$setParameterValues(fitTrain$par, names(fitTrain$par))
      }else{
        stop("Error while optimizing training models. This is often a sign that the sample size of the cross-validation training sets is to small. Consider using a higher k or rerun the analysis with a differnt seed to get a new split of the samples.")
      }

      if(!is.null(subjectSpecificParameters)){
        cvFoldsAndModels$testModels[[i]] <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = cvFoldsAndModels$testSets[[i]], group = 1:nrow(dataset), groupSpecificParameters = subjectSpecificParameters, silent = TRUE))
      }else{
        cvFoldsAndModels$testModels[[i]] <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = cvFoldsAndModels$testSets[[i]], silent = TRUE))
      }


      # set adaptive lasso weights based on training-sample
      if(all(is.na(argsIn$adaptiveLassoWeights[i,]))){
        argsIn$adaptiveLassoWeights[i,] <- getAdaptiveLassoWeights(cpptsemObject = cvFoldsAndModels$trainModels[[i]],
                                                                   penalty = argsIn$penalty,
                                                                   adaptiveLassoWeights = NULL,
                                                                   adaptiveLassoPower =  argsIn$adaptiveLassoPower,
                                                                   standardizeDrift = argsIn$standardizeDrift)[colnames(argsIn$adaptiveLassoWeights)]
      }

      # compute lambda_max
      if(!is.numeric(lambdas) && lambdas == "auto"){
        cat(paste0("Fold ", i, ": "))
        maxLambda <- regCtsem::getMaxLambda(cpptsemObject = cvFoldsAndModels$trainModels[[i]],
                                            objective = argsIn$objective,
                                            regIndicators = argsIn$regIndicators,
                                            targetVector = argsIn$targetVector,
                                            adaptiveLassoWeights = argsIn$adaptiveLassoWeights[i,],
                                            nMultistart = argsIn$nMultistart)
        maxLambdas[i] <- maxLambda$maxLambda
        sparseParameters[names(maxLambda$sparseParameters),i] <- maxLambda$sparseParameters
        argsIn$sparseParameters <- sparseParameters
      }
    }

    if(!is.numeric(lambdas) && lambdas == "auto"){
      if(scaleLambdaWithN){
        sampleSize <- nrow(dataset)
        maxLambdas <- maxLambdas/sampleSize
      }
      maxLambda <- maxLambdas + maxLambdas/25 # adding some wiggle room as there will always be some deviations
      argsIn$lambdas <- seq(0, max(maxLambdas, na.rm = TRUE), length.out = argsIn$lambdasAutoLength)
    }

    ## fit models
    cvFit <- matrix(NA, nrow = argsIn$k+1, length(argsIn$lambdas))
    rownames(cvFit) <- c(paste0("CV", 1:argsIn$k), "mean")
    subModels <- vector("list", argsIn$k)
    names(subModels) <- paste0("CVModel", 1:argsIn$k)

    for(i in 1:argsIn$k){
      message(paste0("\n Fitting CV Model ", i, " of ", argsIn$k, "."))
      regCtsemObject <- regCtsem::exact_regCtsem(cpptsemObject = cvFoldsAndModels$trainModels[[i]],
                                                 dataset = cvFoldsAndModels$trainSets[[i]],
                                                 regIndicators = argsIn$regIndicators,
                                                 targetVector = argsIn$targetVector,
                                                 lambdas = argsIn$lambdas,
                                                 lambdasAutoLength = argsIn$lambdasAutoLength,
                                                 penalty = argsIn$penalty,
                                                 adaptiveLassoWeights = argsIn$adaptiveLassoWeights[i,],
                                                 returnFitIndices = argsIn$returnFitIndices,
                                                 BICWithNAndT = argsIn$BICWithNAndT,
                                                 Tpoints = argsIn$Tpoints,
                                                 cvSampleCpptsemObject = cvFoldsAndModels$testModels[[i]],
                                                 optimizer = argsIn$optimizer,
                                                 objective = argsIn$objective,
                                                 sparseParameters = argsIn$sparseParameters[,i],
                                                 stepSize = argsIn$stepSize,
                                                 lineSearch = argsIn$lineSearch,
                                                 c1 = argsIn$c1,
                                                 c2 = argsIn$c2,
                                                 sig = argsIn$sig,
                                                 gam = argsIn$gam,
                                                 initialHessianApproximation = argsIn$initialHessianApproximation,
                                                 maxIter_out = argsIn$maxIter_out,
                                                 maxIter_in = argsIn$maxIter_in,
                                                 maxIter_line = argsIn$maxIter_line,
                                                 eps_out = argsIn$eps_out,
                                                 eps_in = argsIn$eps_in,
                                                 eps_WW = argsIn$eps_WW,
                                                 eta = argsIn$eta,
                                                 stepsizeMin = argsIn$stepsizeMin,
                                                 stepsizeMax = argsIn$stepsizeMax,
                                                 GISTLinesearchCriterion = argsIn$GISTLinesearchCriterion,
                                                 GISTNonMonotoneNBack = argsIn$GISTNonMonotoneNBack,
                                                 break_outer = argsIn$break_outer,
                                                 scaleLambdaWithN = argsIn$scaleLambdaWithN,
                                                 approxFirst = argsIn$approxFirst,
                                                 numStart = argsIn$numStart,
                                                 nMultistart = argsIn$nMultistart,
                                                 controlApproxOptimizer = argsIn$controlApproxOptimizer,
                                                 verbose = argsIn$verbose)

      fitAndParametersSeparated <- try(separateFitAndParameters(regCtsemObject))
      if(!any(class(fitAndParametersSeparated) == "try-error")){
        cvFit[i,] <- fitAndParametersSeparated$fit["cvM2LL",]
        regCtsemObject$fit <- fitAndParametersSeparated$fit
        regCtsemObject$parameterEstimatesRaw <- fitAndParametersSeparated$parameterEstimates
        regCtsemObject$parameters <- try(getParameterEstimates(regCtsemObject = regCtsemObject,
                                                               parameterEstimatesRaw = fitAndParametersSeparated$parameterEstimates),
                                         silent = TRUE)
        if(any(class(regCtsemObject$parameters) == "try-error")){
          warning("Could not compute the variances and covariances from the DIFFUSIONbase and T0VARbase. You will only find the raw paramter values in the output.")
        }
      }
      subModels[[i]] <- regCtsemObject
    }
    cvFit["mean",] <- apply(cvFit[1:argsIn$k,], 2, mean, na.rm = TRUE)
    regCtsemCVObject <- list("fit" = cvFit, "subModels" = subModels, "cvFoldsAndModels" = cvFoldsAndModels, "setup" = argsIn)
    class(regCtsemCVObject) <- "regCtsemCV"
    cat("\n")
    return(regCtsemCVObject)
  }


  #### Approximate Optimization ####

  #### without cross-validation ####
  if(autoCV == "No" && tolower(optimization) == "approx"){
    regCtsemObject <- regCtsem::approx_regCtsem(cpptsemObject = cpptsemObject,
                                                dataset = argsIn$dataset,
                                                # penalty settings
                                                regIndicators = argsIn$regIndicators,
                                                lambdas = argsIn$lambdas,
                                                lambdasAutoLength = argsIn$lambdasAutoLength,
                                                targetVector = argsIn$targetVector,
                                                penalty = argsIn$penalty,
                                                adaptiveLassoWeights = argsIn$adaptiveLassoWeights,
                                                # fit settings
                                                returnFitIndices = argsIn$returnFitIndices,
                                                BICWithNAndT = argsIn$BICWithNAndT,
                                                Tpoints = argsIn$Tpoints,
                                                cvSampleCpptsemObject = argsIn$cvSampleCpptsemObject,
                                                # optimization settings
                                                objective = argsIn$objective,
                                                epsilon = argsIn$epsilon,
                                                zeroThresh = argsIn$zeroThresh,
                                                controlApproxOptimizer = argsIn$controlApproxOptimizer,
                                                nMultistart = argsIn$nMultistart,
                                                # additional settings
                                                scaleLambdaWithN = argsIn$scaleLambdaWithN,
                                                verbose = argsIn$verbose)

    fitAndParametersSeparated <- try(separateFitAndParameters(regCtsemObject))
    if(!any(class(fitAndParametersSeparated) == "try-error")){
      regCtsemObject$fit <- fitAndParametersSeparated$fit
      regCtsemObject$parameterEstimatesRaw <- fitAndParametersSeparated$parameterEstimates
      regCtsemObject$parameters <- try(getParameterEstimates(regCtsemObject = regCtsemObject,
                                                             parameterEstimatesRaw = fitAndParametersSeparated$parameterEstimates),
                                       silent = TRUE)
      if(any(class(regCtsemObject$parameters) == "try-error")){
        warning("Could not compute the variances and covariances from the DIFFUSIONbase and T0VARbase. You will only find the raw paramter values in the output.")
      }
    }
    regCtsemObject$setup <- append(regCtsemObject$setup, argsIn[!(names(argsIn) %in% names(regCtsemObject$setup))])
    if(!is.null(subjectSpecificParameters)){
      class(regCtsemObject) <- "regCtsemMultiSubject"
    }else{
      class(regCtsemObject) <- "regCtsem"
    }
    cat("\n")
    return(regCtsemObject)
  }

  #### with automatic cross-validation ####
  # generate splits
  if((autoCV == "kFold" || autoCV == "Blocked") && tolower(optimization) == "approx"){
    cvFoldsAndModels <- regCtsem::createCVFoldsAndModels(dataset = argsIn$dataset,
                                                         Tpoints = argsIn$Tpoints,
                                                         manifestNames = ctsemObject$ctmodelobj$manifestNames,
                                                         k = argsIn$k,
                                                         autoCV = autoCV,
                                                         initialPars = ("T0MEANS" %in%  argsIn$cpptsemObject$parameterTable$matrix) || ("T0VARbase" %in%  argsIn$cpptsemObject$parameterTable$matrix))
    if(!is.numeric(lambdas) && lambdas == "auto"){
      maxLambdas <- matrix(NA, nrow = 1, ncol = argsIn$k)
      sparseParameters <- matrix(NA, nrow = length(cpptsemObject$getParameterValues()), ncol = argsIn$k)
      rownames(sparseParameters) <- names(cpptsemObject$getParameterValues())
    }

    # prepare adaptive lasso weights
    if(is.null(argsIn$adaptiveLassoWeights)){
      argsIn$adaptiveLassoWeights <- matrix(NA, nrow = k, ncol = length(cpptsemObject$getParameterValues()))
      rownames(argsIn$adaptiveLassoWeights) <- paste("CV", 1:k)
      colnames(argsIn$adaptiveLassoWeights) <- names(cpptsemObject$getParameterValues())
    }else{
      coln <- names(argsIn$adaptiveLassoWeights)
      argsIn$adaptiveLassoWeights <- matrix(rep(argsIn$adaptiveLassoWeights, k), nrow = k, ncol = length(argsIn$adaptiveLassoWeights), byrow = TRUE)
      colnames(argsIn$adaptiveLassoWeights) <- coln
    }

    for(i in 1:argsIn$k){
      if(!is.null(subjectSpecificParameters)){
        cvFoldsAndModels$trainModels[[i]] <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = cvFoldsAndModels$trainSets[[i]], group = 1:nrow(dataset), groupSpecificParameters = subjectSpecificParameters, silent = TRUE))
        cvFoldsAndModels$trainModels[[i]]$setParameterValues(cpptsemObject$getParameterValues(), names(cpptsemObject$getParameterValues()))
      }else{
        cvFoldsAndModels$trainModels[[i]] <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = cvFoldsAndModels$trainSets[[i]], silent = TRUE))
      }

      # optimize
      invisible(capture.output(fitTrain <- try(optimizeCpptsem(cpptsemObject = cvFoldsAndModels$trainModels[[i]], nMultistart = argsIn$nMultistart),
                                               silent = TRUE), type = c("output", "message")))

      if(!any(class(fitTrain) == "try-error")){
        cvFoldsAndModels$trainModels[[i]]$setParameterValues(fitTrain$par, names(fitTrain$par))
      }else{
        stop("Error while optimizing training models. This is often a sign that the sample size of the cross-validation training sets is too small. Consider using a higher k or rerun the analysis with a different seed to get a new split of the samples.")
      }

      if(!is.null(subjectSpecificParameters)){
        cvFoldsAndModels$testModels[[i]] <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = cvFoldsAndModels$testSets[[i]], group = 1:nrow(dataset), groupSpecificParameters = subjectSpecificParameters, silent = TRUE))
      }else{
        cvFoldsAndModels$testModels[[i]] <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = cvFoldsAndModels$testSets[[i]], silent = TRUE))
      }

      # set adaptive lasso weights based on training-sample
      if(all(is.na(argsIn$adaptiveLassoWeights[i,]))){
        argsIn$adaptiveLassoWeights[i,] <- getAdaptiveLassoWeights(cpptsemObject = cvFoldsAndModels$trainModels[[i]],
                                                                   penalty = argsIn$penalty,
                                                                   adaptiveLassoWeights = NULL,
                                                                   adaptiveLassoPower =  argsIn$adaptiveLassoPower,
                                                                   standardizeDrift = argsIn$standardizeDrift)[colnames(argsIn$adaptiveLassoWeights)]
      }

      # compute lambda_max
      if(!is.numeric(lambdas) && lambdas == "auto"){
        maxLambda <- regCtsem::getMaxLambda(cpptsemObject = cvFoldsAndModels$trainModels[[i]],
                                            objective = argsIn$objective,
                                            regIndicators = argsIn$regIndicators,
                                            targetVector = argsIn$targetVector,
                                            adaptiveLassoWeights = argsIn$adaptiveLassoWeights[i,],
                                            nMultistart = argsIn$nMultistart)
        maxLambdas[i] <- maxLambda$maxLambda
        sparseParameters[names(maxLambda$sparseParameters),i] <- maxLambda$sparseParameters
        argsIn$sparseParameters <- sparseParameters
      }
    }

    if(!is.numeric(lambdas) && lambdas == "auto"){
      if(scaleLambdaWithN){
        sampleSize <- nrow(dataset)
        maxLambdas <- maxLambdas/sampleSize
      }
      maxLambda <- maxLambdas + maxLambdas/25 # adding some wiggle room as there will always be some deviations
      argsIn$lambdas <- seq(0, max(maxLambdas, na.rm = TRUE), length.out = argsIn$lambdasAutoLength)
    }

    ## fit models
    cvFit <- matrix(NA, nrow = argsIn$k+1, length(argsIn$lambdas))
    rownames(cvFit) <- c(paste0("CV", 1:argsIn$k), "mean")
    subModels <- vector("list", argsIn$k)
    names(subModels) <- paste0("CVModel", 1:argsIn$k)

    for(i in 1:argsIn$k){
      message(paste0("\n Fitting CV Model ", i, " of ", argsIn$k, "."))
      regCtsemObject <- regCtsem::approx_regCtsem(cpptsemObject = cvFoldsAndModels$trainModels[[i]],
                                                  dataset = cvFoldsAndModels$trainSets[[i]],
                                                  # penalty settings
                                                  regIndicators = argsIn$regIndicators,
                                                  lambdas = argsIn$lambdas,
                                                  lambdasAutoLength = argsIn$lambdasAutoLength,
                                                  targetVector = argsIn$targetVector,
                                                  penalty = argsIn$penalty,
                                                  adaptiveLassoWeights = argsIn$adaptiveLassoWeights[i,],
                                                  # fit settings
                                                  returnFitIndices = argsIn$returnFitIndices,
                                                  BICWithNAndT = argsIn$BICWithNAndT,
                                                  Tpoints = argsIn$Tpoints,
                                                  cvSampleCpptsemObject = cvFoldsAndModels$testModels[[i]],
                                                  # optimization settings
                                                  objective = argsIn$objective,
                                                  epsilon = argsIn$epsilon,
                                                  zeroThresh = argsIn$zeroThresh,
                                                  nMultistart = argsIn$nMultistart,
                                                  controlApproxOptimizer = argsIn$controlApproxOptimizer,
                                                  # additional settings
                                                  scaleLambdaWithN = argsIn$scaleLambdaWithN,
                                                  verbose = argsIn$verbose)

      fitAndParametersSeparated <- try(separateFitAndParameters(regCtsemObject))
      if(!any(class(fitAndParametersSeparated) == "try-error")){
        cvFit[i,] <- fitAndParametersSeparated$fit["cvM2LL",]
        regCtsemObject$fit <- fitAndParametersSeparated$fit
        regCtsemObject$parameterEstimatesRaw <- fitAndParametersSeparated$parameterEstimates
      }
      subModels[[i]] <- regCtsemObject
    }
    cvFit["mean",] <- apply(cvFit[1:argsIn$k,], 2, mean, na.rm = TRUE)
    regCtsemCVObject <- list("fit" = cvFit, "subModels" = subModels, "cvFoldsAndModels" = cvFoldsAndModels, "setup" = argsIn)
    class(regCtsemCVObject) <- "regCtsemCV"
    cat("\n")
    return(regCtsemCVObject)
  }
  stop("Something went wrong. Check your model specification")
}


#' exact_regCtsem
#'
#' creates a regCtsem object for exact optimization
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param cpptsemObject Object of type cpptsem
#' @param dataset Please provide a data set in wide format compatible to ctsemOMX
#' @param regIndicators Labels for the regularized parameters (e.g. drift_eta1_eta2)
#' @param targetVector named vector with values towards which the parameters are regularized
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01). Alternatively, lambdas can be set to "auto". regCtsem will then compute an upper limit for lambda and test lambdasAutoLength increasing lambda values
#' @param lambdasAutoLength if lambdas == "auto", lambdasAutoLength will determine the number of lambdas tested.
#' @param penalty Currently supported are ridge, lasso, and adaptiveLasso
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the inverse of unregularized parameter estimates.
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param BICWithNAndT Boolean: TRUE = Use N and T in the formula for the BIC (-2log L + log(N+T)*k, where k is the number of parameters in the model). FALSE = Use both N in the formula for the BIC (-2log L + log(N))
#' @param Tpoints Number of time points (used for BICWithNAndT)
#' @param cvSampleCpptsemObject cppstem for cross-validation
#' @param optimizer Either "GIST" or "GLMNET"
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param sparseParameters labeled vector with parameter estimates of the most sparse model. Required for approxFirst = 3. If regValues = "auto" the sparse parameters will be computed automatically.
#' @param stepSize GLMNET & GIST: initial step size of the outer iteration
#' @param lineSearch GLMNET: String indicating which linesearch should be used. Defaults to the one described in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Alternatively (not recommended) Wolfe conditions (lineSearch = "Wolfe") can be used in the outer iteration. Setting to "none" is also not recommended!
#' @param c1 GLMNET: c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
#' @param c2 GLMNET: c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
#' @param sig GLMNET & GIST: GLMNET: only relevant when lineSearch = 'GLMNET' | GIST: sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
#' @param gam GLMNET when lineSearch = 'GLMNET'. Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param initialHessianApproximation GLMNET: Which initial hessian approximation should be used? Possible are: 'identity' for an identity matrix and 'OpenMx' (here the hessian approxmiation from the mxObject is used). If the Hessian from 'OpenMx' is not positive definite, the negative Eigenvalues will be 'flipped' to positive Eigenvalues. "estimate" will estimate the Hessian using optimHess. All of these approaches work most of the time, but not always. Alternatively, a matrix can be provided which will be used as initial Hessian
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
#' @param break_outer criterion for breaking outer iterations of GIST. See ?controlGIST for more information
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended as the likelihood is also sample size dependent
#' @param approxFirst Should approximate optimization be used first to obtain start values for exact optimization?
#' @param numStart Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param controlApproxOptimizer settings passed to optimx or Rsolnp
#' @param nMultistart controls how many different starting values are tried when estimating lambda_max
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#'
#' @author Jannik Orzek
#' @import ctsemOMX rlist
#' @export
exact_regCtsem <- function(  # model
  cpptsemObject,# = NULL,
  dataset,# = NULL,
  # penalty settings
  regIndicators,
  targetVector,
  lambdas,# = "auto",
  lambdasAutoLength,# = 50,
  penalty,# = "lasso",
  adaptiveLassoWeights,# = NULL,
  # fit settings
  returnFitIndices,# = TRUE,
  BICWithNAndT,
  Tpoints,# = NULL,
  cvSampleCpptsemObject,# = NULL,
  # optimization settings
  optimizer,# = "GIST",
  objective,
  sparseParameters,# = NULL,
  # settings for optimization
  stepSize,# = 1,
  lineSearch,# = "GLMNET", # only used in GLMNET
  c1,# = .0001, # c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
  c2,# = .9, # c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
  sig,# = 10^(-5),
  gam,# = 0,
  initialHessianApproximation,# = NULL,
  maxIter_out,# = 100,
  maxIter_in,# = 1000,
  maxIter_line,# = 500,
  eps_out,# = .0000000001,
  eps_in,# = .0000000001,
  eps_WW,# = .0001,
  # settings for GIST
  eta,# = 2,
  stepsizeMin,# = 1/(10^30),
  stepsizeMax,# = 10^30,
  GISTLinesearchCriterion,# = "monotone",
  GISTNonMonotoneNBack,# = 5,
  break_outer,# = c("parameterChange" = 10^(-5)),
  # general
  scaleLambdaWithN,# = TRUE,
  approxFirst,# = FALSE,
  numStart,# = 0,
  nMultistart,
  controlApproxOptimizer,
  # additional settings
  verbose # = 0
){

  exactArgsIn <- as.list(environment())

  returnList <- list("setup" = exactArgsIn)

  parameterLabels <- names(cpptsemObject$getParameterValues())

  returnList$setup$parameterLabels <- parameterLabels

  sampleSize <- nrow(dataset)

  # set lambdas if lambdas == "auto"
  if(any(lambdas == "auto")){
    if(is.null(targetVector)){
      targetVector <- rep(0, length(regIndicators))
      names(targetVector) <- regIndicators
    }
    maxLambda <- getMaxLambda(cpptsemObject = cpptsemObject,
                              objective = objective,
                              regIndicators = regIndicators,
                              targetVector = targetVector,
                              adaptiveLassoWeights = adaptiveLassoWeights,
                              nMultistart = nMultistart)
    sparseParameters <- maxLambda$sparseParameters
    maxLambda <- maxLambda$maxLambda + maxLambda$maxLambda/25 # adding some wiggle room as there will always be some deviations
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

  fits <- c("regM2LL", "m2LL")
  if(returnFitIndices){
    fits <- c(fits, "AIC", "BIC", "estimatedParameters")
  }
  if(!is.null(cvSampleCpptsemObject)){
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
    regModel <- try(regCtsem::exact_bfgsGLMNET(cpptsemObject = cpptsemObject, dataset = dataset, objective = objective,
                                               regIndicators = regIndicators, lambdas = lambdas,
                                               adaptiveLassoWeights = adaptiveLassoWeights,
                                               # additional settings
                                               sparseParameters = sparseParameters,
                                               stepSize = stepSize, lineSearch = lineSearch, c1 = c1, c2 = c2,
                                               sig = sig, gam = gam,
                                               initialHessianApproximation = initialHessianApproximation,
                                               maxIter_out = maxIter_out,
                                               maxIter_in = maxIter_in, maxIter_line = maxIter_line,
                                               eps_out = eps_out, eps_in = eps_in, eps_WW = eps_WW,
                                               scaleLambdaWithN = scaleLambdaWithN, sampleSize = sampleSize,
                                               approxFirst = approxFirst,
                                               numStart = numStart, controlApproxOptimizer = controlApproxOptimizer,
                                               verbose = verbose))
  }
  if(tolower(optimizer) == "gist"){
    regModel <- try(regCtsem::exact_GIST(cpptsemObject = cpptsemObject, dataset = dataset, objective = objective,
                                         regIndicators = regIndicators, targetVector = targetVector,
                                         lambdas = lambdas, adaptiveLassoWeights = adaptiveLassoWeights,
                                         # additional settings
                                         sparseParameters = sparseParameters,
                                         eta = eta, sig = sig,
                                         initialStepsize = stepSize, stepsizeMin = stepsizeMin, stepsizeMax = stepsizeMax,
                                         GISTLinesearchCriterion = GISTLinesearchCriterion, GISTNonMonotoneNBack = GISTNonMonotoneNBack,
                                         maxIter_out = maxIter_out, maxIter_in = maxIter_in,
                                         break_outer = break_outer,
                                         scaleLambdaWithN = scaleLambdaWithN, sampleSize = sampleSize,
                                         approxFirst = approxFirst,
                                         numStart = numStart, controlApproxOptimizer = controlApproxOptimizer,
                                         verbose = verbose
    )
    )
  }
  if(!any(class(regModel) == "try-error")){
    # save results
    fitAndParameters[parameterLabels,] <- regModel$thetas
    fitAndParameters["m2LL",] <- regModel$m2LL
    fitAndParameters["regM2LL",] <- regModel$regM2LL
  }

  # save fit indices
  if(returnFitIndices & !(any(class(regModel) == "try-error"))){
    if(is.null(targetVector) || all(targetVector == 0)){
      fitIndicesTable <- regCtsem::exact_getFitIndices(parameterLabels = parameterLabels,
                                                       fitAndParameters = fitAndParameters, lambdas = lambdas,
                                                       sampleSize = ifelse(BICWithNAndT,Tpoints*sampleSize, sampleSize))
    }else{
      fitIndicesTable <- regCtsem::exact_getFitIndicesWithTarget(parameterLabels = parameterLabels, regIndicators = regIndicators,
                                                                 fitAndParameters = fitAndParameters, targetVector = targetVector,
                                                                 lambdas = lambdas,
                                                                 sampleSize = ifelse(BICWithNAndT,Tpoints*sampleSize, sampleSize))

    }
    fitAndParameters[rownames(fitIndicesTable),] <- fitIndicesTable[rownames(fitIndicesTable),]
  }


  # cross-validation
  if(!is.null(cvSampleCpptsemObject)){

    fitCVTable <- regCtsem::exact_getCVFit(objective = objective, cvSampleCpptsemObject = cvSampleCpptsemObject, parameterLabels = parameterLabels,
                                           parameterValuesTable = as.matrix(fitAndParameters[parameterLabels,], ncol = length(lambdas)),
                                           lambdas = lambdas)

    fitAndParameters["cvM2LL",] <- fitCVTable["cvM2LL",]
  }

  returnList <- rlist::list.append(returnList, "fitAndParameters" = fitAndParameters)
  return(returnList)
}


#' approx_regCtsem
#'
#' creates a regCtsem object for approximate optimization
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param cpptsemObject Fitted object of class cpptsem
#' @param dataset Please provide a data set in wide format compatible to ctsemOMX
#' @param regIndicators Labels for the regularized parameters (e.g. drift_eta1_eta2)
#' @param targetVector named vector with values towards which the parameters are regularized
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01). Alternatively, lambdas can be set to "auto". regCtsem will then compute an upper limit for lambda and test lambdasAutoLength increasing lambda values
#' @param lambdasAutoLength if lambdas == "auto", lambdasAutoLength will determine the number of lambdas tested.
#' @param penalty Currently supported are ridge, lasso, and adaptiveLasso
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the inverse of unregularized parameter estimates.
#' @param returnFitIndices Boolean: should fit indices be returned?
#' @param BICWithNAndT Boolean: TRUE = Use N and T in the formula for the BIC (-2log L + log(N+T)*k, where k is the number of parameters in the model). FALSE = Use both N in the formula for the BIC (-2log L + log(N))
#' @param Tpoints Number of time points (used for BICWithNAndT)
#' @param cvSampleCpptsemObject cppstem for cross-validation
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
#' @param zeroThresh threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
#' @param nMultistart controls how many different starting values are tried when estimating lambda_max
#' @param controlApproxOptimizer settings passed to optimx or Rsolnp
#' @param scaleLambdaWithN Boolean: Should the penalty value be scaled with the sample size? True is recommended, as the likelihood is also sample size dependent
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#'
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
approx_regCtsem <- function(  # model
  cpptsemObject,# = NULL,
  dataset,# = NULL,
  # penalty settings
  regIndicators,
  lambdas,
  lambdasAutoLength,# = 50,
  targetVector,
  penalty,# = "lasso",
  adaptiveLassoWeights,# = NULL,
  # fit settings
  returnFitIndices,# = TRUE,
  BICWithNAndT,
  Tpoints,# = NULL,
  cvSampleCpptsemObject,# = NULL,
  # optimization settings
  objective,# = "ML",
  epsilon,# = .001,
  zeroThresh,# = .001,
  controlApproxOptimizer,
  nMultistart,
  # additional settings
  scaleLambdaWithN,# = TRUE,
  verbose# = 0
){

  approxArgsIn <- as.list(environment())

  returnList <- list("setup" = approxArgsIn)

  parameterLabels <- names(cpptsemObject$getParameterValues())

  returnList$setup$parameterLabels <- parameterLabels

  sampleSize <- nrow(dataset)

  # set lambdas if lambdas == "auto"
  if(any(lambdas == "auto")){
    if(is.null(targetVector)){
      targetVector <- rep(0, length(regIndicators))
      names(targetVector) <- regIndicators
    }
    maxLambda <- getMaxLambda(cpptsemObject = cpptsemObject,
                              objective = objective,
                              regIndicators = regIndicators,
                              targetVector = targetVector,
                              adaptiveLassoWeights = adaptiveLassoWeights,
                              nMultistart = nMultistart)
    sparseParameters <- maxLambda$sparseParameters
    maxLambda <- maxLambda$maxLambda + maxLambda$maxLambda/25 # adding some wiggle room as there will always be some deviations
    lambdas <- seq(0,maxLambda, length.out = lambdasAutoLength)
    # if scaleLambdaWithN = TRUE, the lambda values will be multiplied with N later on
    # we need to ensure that this will not change the values:
    if(scaleLambdaWithN){
      lambdas <- lambdas/sampleSize
    }
    returnList$setup$lambdas <- lambdas
    approxArgsIn$lambdas <- lambdas
    returnList$setup$sparseParameters <- sparseParameters
    approxArgsIn$sparseParameters <- sparseParameters
  }


  # start fitting models
  fitAndParameters <- try(regCtsem::approx_iterateOverLambdas(cpptsemObject = cpptsemObject, dataset = dataset, sampleSize = sampleSize,
                                                              # penalty settings
                                                              regIndicators = regIndicators, lambdas = lambdas,
                                                              penalty = penalty, adaptiveLassoWeights = adaptiveLassoWeights,
                                                              targetVector = targetVector,
                                                              # fit settings
                                                              returnFitIndices = returnFitIndices, BICWithNAndT = BICWithNAndT, Tpoints = Tpoints,
                                                              # optimization settings
                                                              objective = objective, epsilon = epsilon, zeroThresh = zeroThresh, controlApproxOptimizer = controlApproxOptimizer,
                                                              # additional settings
                                                              scaleLambdaWithN = scaleLambdaWithN, verbose = verbose))

  # cross-validation
  if(!is.null(cvSampleCpptsemObject)){

    fitCVTable <- regCtsem::exact_getCVFit(objective = objective, cvSampleCpptsemObject = cvSampleCpptsemObject, parameterLabels = parameterLabels,
                                           parameterValuesTable = as.matrix(fitAndParameters[parameterLabels,], ncol = length(lambdas)),
                                           lambdas = lambdas)

    fitAndParameters <- rbind(fitCVTable, fitAndParameters)
  }

  return(rlist::list.append(returnList, "fitAndParameters" = fitAndParameters))

}


#### Cross - Validation ####

#' createCVFoldsAndModels
#'
#' creates cv folds and models
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param dataset Please provide a data set in wide format compatible to ctsemOMX
#' @param Tpoints Number of time points
#' @param manifestNames names of manifest variables
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @param autoCV Form of cross-validation: "kFold" or "Blocked"
#' @param initialPars are any of the initial parameters (T0MEANS or T0VAR) estimated? If so, 10 percent of the initial observations will not be used in cross-validation when autoCV = "Blocked". Otherwise the initial parameters will be very poorly estimated
#' @author Jannik Orzek
#' @export
createCVFoldsAndModels <- function(dataset, Tpoints, manifestNames, k, autoCV, initialPars){
  fullData <- dataset

  if(autoCV == "kFold"){
    cvFolds <- regCtsem::createFolds(nrow(fullData), k = k)

    testSets <- vector("list", length = k)
    names(testSets) <- names(cvFolds)
    for(foldNumber in 1:k){
      if(is.vector(fullData[cvFolds[[foldNumber]],])){
        # if a single row is selected, R extracts this row as vector, not as matrix
        testSets[[foldNumber]] <-  t(as.matrix(fullData[cvFolds[[foldNumber]],]))
      }else{
        testSets[[foldNumber]] <-  fullData[cvFolds[[foldNumber]],]
      }
    }

    trainSets <- vector("list", length = k)
    names(trainSets) <- names(cvFolds)
    for(foldNumber in 1:k){
      if(is.vector(fullData[-cvFolds[[foldNumber]],])){
        # if a single row is selected, R extracts this row as vector, not as matrix
        trainSets[[foldNumber]] <-  t(as.matrix(fullData[-cvFolds[[foldNumber]],]))
      }else{
        trainSets[[foldNumber]] <- fullData[-cvFolds[[foldNumber]],]
      }
    }

    trainModels <- vector("list", length = k)
    names(trainModels) <- names(cvFolds)

    testModels <- vector("list", length = k)
    names(testModels) <- names(cvFolds)

    return(list("fullData" = fullData, "cvFolds" = cvFolds, "testSets" = testSets, "trainSets" = trainSets, "trainModels" = trainModels, "testModels" = testModels))
  }

  if(autoCV == "Blocked"){
    warning("Blocked CV naively builds blocks of the data of each individual and does not account for dependencies which still exist in the training and test set (see Bergmeir, C., & Benítez, J. M. (2012). On the use of cross-validation for time series predictor evaluation. Information Sciences, 191, 192–213. https://doi.org/10.1016/j.ins.2011.12.028 and Bulteel, K., Mestdagh, M., Tuerlinckx, F., & Ceulemans, E. (2018). VAR(1) based models do not always outpredict AR(1) models in typical psychological applications. Psychological Methods, 23(4), 740–756. https://doi.org/10.1037/met0000178). Also, different time intervals are not accounted for. More sophisticated forms of cross-validation currently have to be implemented manually.")
    if(initialPars){
      # If initial parameters are estimated we will force the first few observations to always be used in the training
      # set. This is to ensure that the estimates for the initial parameters are somewhat more reliable.
      Folds <- (max(1, Tpoints %/% 10)):(Tpoints-1)
      #ignore <- paste0(paste0(manifestNames, "_T"), rep(0:(max(1, Tpoints %/% 10)-1), each = length(manifestNames)))
      warning(paste0("Bergmeir & Benítez, (2012) do not recommend the use of cross-validation with non-stationary time series! However, your model seems to not assume stationarity. Keep this in mind when interpreting the results. max(1, Tpoints %/% 10) = ", max(1, Tpoints %/% 10), " of the initial observations will not be used in blocked CV. Otherwise the initial parameters will be very poorly estimated."))
    }else{
      #ignore  <- c()
      # If stationarity is assumed, all time points are used in the cross-validation. No special treatment
      # for the first few observations.
      Folds <- 0:(Tpoints-1)
    }


    cvFolds <- split(Folds, cut(x = 1:length(Folds), breaks = k, labels = paste0("fold", 1:k)))

    testSets <- vector("list", length = k)
    names(testSets) <- names(cvFolds)
    for(foldNumber in 1:k){
      # select all data except the test set
      testColumns <- paste0(paste0(manifestNames, "_T"), rep(cvFolds[[foldNumber]], each = length(manifestNames)))
      dTs <- colnames(fullData)[grepl("dT", colnames(fullData))]
      testSet <- fullData
      testSet[,!(colnames(testSet) %in% c(testColumns, dTs))] <- NA # remove all but the test set
      testSets[[foldNumber]] <-  testSet
    }

    trainSets <- vector("list", length = k)
    names(trainSets) <- names(cvFolds)
    for(foldNumber in 1:k){
      testColumns <- paste0(paste0(manifestNames, "_T"), rep(cvFolds[[foldNumber]], each = length(manifestNames)))
      trainSet <- fullData
      trainSet[,testColumns] <- NA # remove test set
      trainSets[[foldNumber]] <- trainSet
    }

    trainModels <- vector("list", length = k)
    names(trainModels) <- names(cvFolds)

    testModels <- vector("list", length = k)
    names(testModels) <- names(cvFolds)

    return(list("fullData" = fullData, "cvFolds" = cvFolds, "testSets" = testSets, "trainSets" = trainSets, "trainModels" = trainModels, "testModels" = testModels))
  }
  stop("Error in createCVFoldsAndModels: autoCV must be set to 'kFold' or 'Blocked'")
}


#' createFolds
#'
#' created folds for automatic cross-validation
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param sampleSize sample size
#' @param k number of cross-validation folds (k-fold cross-validation)
#' @author Jannik Orzek
#' @export
createFolds <- function(sampleSize, k){
  # shuffle
  Folds <- sample(1:sampleSize, sampleSize)

  Folds <- split(Folds, cut(x = 1:length(Folds), breaks = k, labels = paste0("fold", 1:k)))

  return(Folds)
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
  warning("NOT YET ADJUSTED FOR NEW IMPLEMENTATION")
  if(is.null(argsIn$ctsemObject) & is.null(argsIn$mxObject)){
    stop("Both ctsemObject and mxObect are missing. You need to provide at least one")
  }
  if(!any(class(argsIn$ctsemObject) == "ctsemFit")){
    stop("ctsemObject has to be of class ctsemFit")
  }

  #if(!is.null(argsIn$mxObject) & any(class(argsIn$mxObject)=="MxModel")){
  #  # check if ctsemObject is estimated with Kalman
  #  if(any(class(argsIn$mxObject$expectation) == "MxExpectationStateSpace") &  !any(class(argsIn$ctsemObject)=="ctsemInit")){
  #    stop("It seems like the provided mxObject was fitted with Kalman filter. To use the Kalman filter, provide the object of type ctsemInit from ctModel instead of the fitted model. Set the objective = 'Kalman' and provide a dataset")
  #  }
  #}

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
    #if(!any(class(argsIn$ctsemObject) == "ctsemInit")){stop("ctsemObject has to be of class ctsemInit. You can get the ctsemInit object from ctModel.")}
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


#' getAdaptiveLassoWeights
#'
#' Computes the adaptive lasso weights
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param cpptsemObject Fitted object of class cpptsem
#' @param penalty type
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param adaptiveLassoPower power of the adaptive lasso weights
#' @param standardizeDrift Should Drift parameters be standardized automatically? Set to 'No' for no standardization, 'T0VAR' for standardization using the T0VAR or 'asymptoticDiffusion' for standardization using the asymptotic diffusion
#' @author Jannik Orzek
#' @export
getAdaptiveLassoWeights <- function(cpptsemObject, penalty, adaptiveLassoWeights, adaptiveLassoPower, standardizeDrift){
  if(!(standardizeDrift == "No" || standardizeDrift == "T0VAR" || standardizeDrift == "asymptoticDiffusion" )){
    stop("standardizeDrift has to be set to 'No', 'T0VAR' or 'asymptoticDiffusion'")
  }

  if(standardizeDrift == "No"){
    warning("Not using any automatic standardization.")
  }

  # if adaptiveLassoWeights were provided and no standardization is requested:
  if(tolower(penalty) == "adaptivelasso" && is.numeric(adaptiveLassoWeights) && (standardizeDrift == "No")){
    return(adaptiveLassoWeights)
  }

  # otherwise: set adaptiveLassoWeights
  ## if automatic drift standardization was requested
  if(standardizeDrift == "T0VAR" || standardizeDrift == "asymptoticDiffusion"){
    cat(paste0("Standardizing drift parameters with ", standardizeDrift, "\n"))

    # check if adaptiveLassoWeights were provided
    if(is.numeric(adaptiveLassoWeights)){ stop("standardizeDrift and provided adaptiveLassoWeights can not be combined automatically. Consider setting adaptiveLassoWeights = 'auto'") }

    thetaNames <- names(cpptsemObject$getParameterValues())

    # check if lasso or adaptiveLasso
    if(tolower(penalty) == "lasso" || tolower(penalty) == "ridge"){
      adaptiveLassoWeights <- rep(1,length(thetaNames))
      names(adaptiveLassoWeights) <- thetaNames
    }else if(tolower(penalty) == "adaptivelasso"){
      adaptiveLassoWeights <- abs(cpptsemObject$getParameterValues())^(-1)
    }else{
      stop(paste0("Unexpected penalty: ", penalty, ". Possible are lasso, ridge and adaptiveLasso."))
    }

    # compute standardizers
    if(standardizeDrift == "T0VAR"){
      VARIs <- cpptsemObject$T0VARValues[]
    }
    if(standardizeDrift == "asymptoticDiffusion"){
      VARIs <- cpptsemObject$asymptoticDIFFUSION[]
    }
    DRIFTS <- cpptsemObject$parameterTable[cpptsemObject$parameterTable$matrix == "DRIFT",]
    driftLabels <- matrix("", nrow = nrow(cpptsemObject$DRIFTValues), ncol = ncol(cpptsemObject$DRIFTValues))
    for(i in 1:nrow(DRIFTS)){
      driftLabels[DRIFTS$row[i]+1, DRIFTS$col[i]+1] <- DRIFTS$label[i]
    }

    if(anyNA(driftLabels)){
      autoDriftLabels <- matrix(paste0("_autoDriftLabel_", rep(seq_len(nrow(driftLabels)), each = ncol(driftLabels)), "_", seq_len(ncol(driftLabels))),
                                nrow = nrow(driftLabels), ncol = ncol(driftLabels), byrow = TRUE)
      driftLabels[is.na(driftLabels)] <- autoDriftLabels[is.na(driftLabels)]
    }

    flatStandardizers <- regCtsem::getFlatStdizer(VARIs = VARIs, driftLabels = driftLabels)
    flatStandardizers <- flatStandardizers[rownames(flatStandardizers) %in% thetaNames,]

    for(thetaName in names(flatStandardizers)){
      adaptiveLassoWeights[thetaName] <- flatStandardizers[thetaName]*adaptiveLassoWeights[thetaName]
    }
    if(tolower(penalty) == "adaptivelasso"){
      # use exponent
      adaptiveLassoWeights <- adaptiveLassoWeights^(abs(adaptiveLassoPower))
    }

    return(adaptiveLassoWeights)
  }

  ## if lasso was requested, but no automatic standardization
  if(tolower(penalty) == "lasso"){
    thetaNames <- names(cpptsemObject$getParameterValues())
    adaptiveLassoWeights <- rep(1,length(thetaNames))
    names(adaptiveLassoWeights) <- thetaNames
    return(adaptiveLassoWeights)
  }

  ## if ridge was requested, but no automatic standardization
  if(tolower(penalty) == "ridge"){
    thetaNames <- names(cpptsemObject$getParameterValues())
    adaptiveLassoWeights <- rep(1,length(thetaNames))
    names(adaptiveLassoWeights) <- thetaNames
    return(adaptiveLassoWeights)
  }

  ## if adaptivelasso lasso was requested, but no automatic standardization
  if(tolower(penalty) == "adaptivelasso" && is.null(adaptiveLassoWeights)){
    adaptiveLassoWeights <- abs(cpptsemObject$getParameterValues())^(-abs(adaptiveLassoPower))
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
getFinalParameters <- function(regCtsemObject, criterion, raw = TRUE){
  if(regCtsemObject$setup$autoCV == "No"){
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
  minCriterionValue <- max(which(regCtsemObject$fit["mean",] == min(regCtsemObject$fit["mean",], na.rm = TRUE)))
  lambdas <- regCtsemObject$setup$lambdas
  bestLambda <- lambdas[minCriterionValue]
  return(list("criterion" = "mean",
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
getFinalModel <- function(regCtsemObject, criterion){
  if(regCtsemObject$setup$autoCV != "No"){
    stop("getFinalModel not supported for automatic cross-validation. At the moment, you have to manually re-run the model with the best lambda value using the whole sample.")
  }

  bestPars <- getFinalParameters(regCtsemObject, criterion = criterion, raw = TRUE)
  message(paste0("Best fit for ", criterion, " was observed for lambda = ", bestPars$lambda, "."))
  regCtsemObject$setup$cpptsemObject$setParameterValues(bestPars$parameters, names(bestPars$parameters))
  if(tolower(regCtsemObject$setup$objective) == "ml"){
    regCtsemObject$setup$cpptsemObject$computeRAM()
    regCtsemObject$setup$cpptsemObject$fitRAM()
  }else{
    regCtsemObject$setup$cpptsemObject$computeAndFitKalman()
  }

  return(regCtsemObject$setup$cpptsemObject)
}


#' getFlatStdizer
#'
#' returns the standardizer for the standardized drift parameters. Computes SD_predictor/SD_dependent
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param VARIs matrix with variances for standardization values
#' @param driftLabels vector with drift names
#' @export
getFlatStdizer <- function(VARIs, driftLabels){
  stdizer <- matrix(1, nrow = nrow(VARIs), ncol = ncol(VARIs))
  stdizer <- stdizer%*%diag(sqrt(diag(VARIs))) # times predictor sd
  stdizer <- t(t(stdizer)%*%diag(sqrt(diag(VARIs))^(-1))) # divided by dependent sd
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
#' @param cpptsemObject Fitted object of class cpptsem
#' @param objective ML or Kalman
#' @param regIndicators Labels of the regularized parameters (e.g. drift_eta1_eta2)
#' @param targetVector vector with target values
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the unregularized parameter estimates.
#' @param nMultistart number of multi-start iterations
#' @author Jannik Orzek
#' @import OpenMx
#' @export
getMaxLambda <- function(cpptsemObject, objective, regIndicators, targetVector, adaptiveLassoWeights, nMultistart){
  # This function is adapted from Murphy (2012) Machine learning: a probabilistic perspective. See p. 434 for more details.
  cat("Computing lambda_max ... ")
  converged <- FALSE

  it <- 0
  # save parameters
  parameters <- cpptsemObject$getParameterValues()
  while(!converged){
    if(it == length(regIndicators)){
      stop("Error while automatically setting the lambdas: The models did not converge. Try setting the lambdas manually.")
    }

    # change parameters to target values:
    param <- parameters

    numberRegularized <- length(regIndicators) - it

    # Problem: setting all regularized parameters to zero might result in an impossible model
    # if this error occurs, regCtsem will iteratively try to set a subset of the parameters to zero
    # the size of the parameters determines the order in which they are set to zero
    # however, this is only a rough approximation and might result in an unsatisfactory maxLambda
    regIndicatorsCurrent <- names(sort(abs(param[regIndicators]))[1:numberRegularized])

    # step 1: set the regularized parameters to their respective target and estimate the model:

    param[regIndicatorsCurrent] <- targetVector[regIndicatorsCurrent]
    freeParam <- rep(TRUE, length(param))
    names(freeParam) <- names(param)
    freeParam[regIndicatorsCurrent] <- FALSE

    cpptsemObject$setParameterValues(param, names(param))

    # check if start values are feasible
    tryFit <- fitCpptsem(parameterValues = param[freeParam],
                         cpptsemObject = cpptsemObject,
                         objective = objective,
                         free = freeParam,
                         failureReturns = NA)
    if(is.na(tryFit)){
      # set to defaults of ctsemOMX
      paramterTable <- cpptsemObject$parameterTable
      for(p in 1:length(paramterTable$label)){
        if(!freeParam[paramterTable$label[p]]){next}
        if(paramterTable$matrix[p] == "DRIFT"){
          if(paramterTable$row[p] == paramterTable$col[p]){
            param[paramterTable$label[p]] <- -.45
          }else{
            param[paramterTable$label[p]] <- -.05
          }
        }
        if(paramterTable$matrix[p] == "DIFFUSIONbase"){
          if(paramterTable$row[p] == paramterTable$col[p]){
            param[paramterTable$label[p]] <- log(10)
          }else{
            param[paramterTable$label[p]] <- 0
          }
        }
        if(paramterTable$matrix[p] == "T0VARbase"){
          if(paramterTable$row[p] == paramterTable$col[p]){
            param[paramterTable$label[p]] <- log(1)
          }else{
            param[paramterTable$label[p]] <- 0
          }
        }
        if(paramterTable$matrix[p] == "MANIFESTMEANS"){
          param[paramterTable$label[p]] <- 0
        }
      }
      tryFit <- fitCpptsem(parameterValues = param[freeParam],
                           cpptsemObject = cpptsemObject,
                           objective = objective,
                           free = freeParam,
                           failureReturns = NA)
      if(is.na(tryFit)){
        if(it ==  0){
          warning("Error when determining the lambdas automatically: Setting all regularized parameters to their target values resulted in an impossible model. regCtsem will try to at least set a subset of the regularized parameters to their target; however, this might result in a wrong maximum for the lambdas! Consider setting the lambdas manually.")
        }
        it <- it + 1
        next
      }
    }

    # optimize
    sparseModel <- try(optimizeCpptsem(cpptsemObject = cpptsemObject, free = freeParam, nMultistart = nMultistart),
                       silent = TRUE)

    if(any(class(sparseModel) == "try-error")){
      if(it ==  0){
        warning("Error when determining the lambdas automatically: Setting all regularized parameters to their target values resulted in an impossible model. regCtsem will try to at least set a subset of the regularized parameters to their target; however, this might result in a wrong maximum for the lambdas! Consider setting the lambdas manually.")
      }
      it <- it + 1
      next
    }

    nonZeroParam <- sparseModel$pars
    namesNonZeroParam <- names(nonZeroParam)

    # step 2: compute gradients with regularized parameters set to target and unregularized parameters set to nonZeroParam estimates
    param[namesNonZeroParam] <- nonZeroParam
    cpptsemObject$setParameterValues(param, names(param))
    grad <- regCtsem::exact_getCppGradients(cpptsemObject, objective = objective)

    if(any(class(grad) == "try-error") || anyNA(grad)){
      if(it ==  0){
        warning("Error when determining the lambdas automatically: Setting all regularized parameters to zero resulted in an impossible model. regCtsem will try to at least set a subset of the regularized parameters to zero; however, this might result in a wrong maximum for the lambdas! Consider setting the lambdas manually.")
      }
      it <- it + 1
      next
    }

    # extract the gradient
    gradLabels <- rownames(grad)
    converged <- TRUE
  }
  if(it > 0){
    warning(paste0("regCtsem did set ", numberRegularized, " of the ", length(regIndicators), " regularized parameters to zero when determining the maximal lambda."))
  }
  # define maxLambda as the maximal gradient of the regularized parameters
  maxLambda <- max(abs(grad[regIndicators]) * adaptiveLassoWeights[regIndicators]^(-1))

  cat("DONE \n")

  # reset parameters
  cpptsemObject$setParameterValues(parameters, names(parameters))
  f <- fitCpptsem(cpptsemObject, objective = objective, parameterValues = parameters, failureReturns = NA)
  return(list("maxLambda" = maxLambda, "sparseParameters" = param))
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


#' getParameterEstimates
#'
#' computes the parameters given the raw parameter estimates
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param regCtsemObject regCtsemObject
#' @param parameterEstimatesRaw raw parameter estimates from regularized ctsem
#' @export
getParameterEstimates <- function(regCtsemObject, parameterEstimatesRaw){
  parameterTable <- regCtsemObject$setup$cpptsemObject$parameterTable
  parameterTable[] <- parameterTable[]
  nlatent <- nrow(regCtsemObject$setup$cpptsemObject$DRIFTValues)
  nmanifest <- nrow(regCtsemObject$setup$cpptsemObject$MANIFESTVARValues)
  if(!is.null(parameterTable$groupID) && length(unique(parameterTable$groupID))>1){
    parameterEstimatesList <- vector("list", length(unique(parameterTable$groupID)))
    names(parameterEstimatesList) <- paste0("group", unique(parameterTable$groupID))
    for(group in unique(parameterTable$groupID)){
      groupParameterTable <- parameterTable[parameterTable$groupID == group,]
      #groupParameterTable$label <- unlist(strsplit(groupParameterTable$label, paste0("_G", group)))
      parameterEstimates <- c()
      for(parMat in unique(groupParameterTable$matrix)){
        if(grepl("base", parMat)){
          matName <- strsplit(parMat, "base")[[1]]
          if(matName == "MANIFESTVAR"){
            parameterEstimates <- rbind(parameterEstimates,
                                        getVariances(parameterEstimatesRaw = parameterEstimatesRaw,
                                                     matName = matName,
                                                     baseMatName = parMat,
                                                     parameterTable = groupParameterTable,
                                                     nVariables = nmanifest,
                                                     variableNames = "Y")
            )
          }else{
            parameterEstimates <- rbind(parameterEstimates,
                                        getVariances(parameterEstimatesRaw = parameterEstimatesRaw,
                                                     matName = matName,
                                                     baseMatName = parMat,
                                                     parameterTable = groupParameterTable,
                                                     nVariables = nlatent,
                                                     variableNames = "eta")
            )
          }
        }else{
          parameterEstimates <- rbind(parameterEstimates,
                                      t(t(parameterEstimatesRaw[groupParameterTable$label[groupParameterTable$matrix == parMat],]))
          )
        }
      }
      colnames(parameterEstimates) <- colnames(parameterEstimatesRaw)
      parameterEstimates <- parameterEstimates[sort(rownames(parameterEstimates)), ]

      parameterEstimatesList[[paste0("group", group)]] <- parameterEstimates
    }
    return(parameterEstimatesList)
  }

  parameterEstimates <- c()
  for(parMat in unique(parameterTable$matrix)){
    if(grepl("base", parMat)){
      matName <- strsplit(parMat, "base")[[1]]
      if(matName == "MANIFESTVAR"){
        parameterEstimates <- rbind(parameterEstimates,
                                    getVariances(parameterEstimatesRaw = parameterEstimatesRaw,
                                                 matName = matName,
                                                 baseMatName = parMat,
                                                 parameterTable = parameterTable,
                                                 nVariables = nmanifest,
                                                 variableNames = "Y")
        )
      }else{
        parameterEstimates <- rbind(parameterEstimates,
                                    getVariances(parameterEstimatesRaw = parameterEstimatesRaw,
                                                 matName = matName,
                                                 baseMatName = parMat,
                                                 parameterTable = parameterTable,
                                                 nVariables = nlatent,
                                                 variableNames = "eta")
        )
      }
    }else{
      parameterEstimates <- rbind(parameterEstimates,
                                  t(t(parameterEstimatesRaw[parameterTable$label[parameterTable$matrix == parMat],]))
      )
    }
  }
  colnames(parameterEstimates) <- colnames(parameterEstimatesRaw)
  parameterEstimates <- parameterEstimates[sort(rownames(parameterEstimates)), ]
  return(parameterEstimates)
}

#' getVariances
#'
#' computes the variances for varianceBase matrices
#'
#' NOTE: Function located in file regCtsem.R
#'
#' @param parameterEstimatesRaw raw parameter estimates from regularized ctsem
#' @param matName name of the new matrix
#' @param baseMatName name of the base matrix
#' @param parameterTable parameterTable
#' @param nVariables number of variables in the matrix
#' @param variableNames names of the variables
#' @export
getVariances <- function(parameterEstimatesRaw, matName, baseMatName, parameterTable, nVariables, variableNames){
  baseMat <- diag(-999, nrow = nVariables, ncol = nVariables)
  VARLabels <- matrix(paste0(matName,"_", variableNames,seq_len(nVariables),
                             rep(paste0("_", variableNames,seq_len(nVariables)), each = nVariables)
  ), nrow = nVariables, ncol = nVariables, byrow = FALSE)
  varPars <- matrix(NA, nrow = sum(lower.tri(VARLabels, diag = TRUE)), ncol = ncol(parameterEstimatesRaw))
  rownames(varPars) <- VARLabels[lower.tri(VARLabels, diag = TRUE)]
  for(lambda in 1:ncol(parameterEstimatesRaw)){
    for(parLab in parameterTable$label[parameterTable$matrix == baseMatName]){
      baseMat[parameterTable$row[parameterTable$label == parLab]+1, parameterTable$col[parameterTable$label == parLab]+1] <- parameterEstimatesRaw[parLab,1]
    }
    VAR <- getVarianceFromVarianceBase2(baseMat)
    varPars[,lambda] <- VAR[lower.tri(VARLabels, diag = TRUE)]
  }
  return(varPars)
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

  parameterLabels <- names(regCtsemObject$setup$cpptsemObject$getParameterValues())
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
       outer = FALSE,
       line=1,col="black",col.ticks="black",col.axis="black")
  mtext("# zeroed parameters",3,line=3,at=mean(lambdas),col="black", cex = 1)

  par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
}



