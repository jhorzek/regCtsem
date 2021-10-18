#' restore
#'
#' restore the cpptsem object in a regCtsem object.
#'
#' After saving and loading a regCtsem object, the underlying C++ model will be lost. This function resores that model
#' @param regCtsemObject object of type regCtsem
#' @export
restore <- function(regCtsemObject){
  if(regCtsemObject$setup$autoCV != "No"){
    stop("Restoring currently not supported for models with automatic cross-validation")
  }
  if(is.null(regCtsemObject$setup$subjectSpecificParameters)){
    regCtsemObject$setup$cpptsemObject <- cpptsemFromCtsem(ctsemModel = regCtsemObject$setup$ctsemObject, wideData = regCtsemObject$setup$dataset, removeD = TRUE)
  }else{
    regCtsemObject$setup$cpptsemObject <- cpptsemFromCtsem(ctsemModel = regCtsemObject$setup$ctsemObject, wideData = regCtsemObject$setup$dataset, removeD = TRUE, group = seq(1,nrow(regCtsemObject$setup$dataset)), groupSpecificParameters = regCtsemObject$setup$subjectSpecificParameters)
  }
  if(!is.null(regCtsemObject$parameterEstimatesRaw)){
    regCtsemObject$setup$cpptsemObject$setParameterValues(regCtsemObject$parameterEstimatesRaw[,1], rownames(regCtsemObject$parameterEstimatesRaw))
  }

  if(regCtsemObject$setup$objective == "Kalman"){
    regCtsemObject$setup$cpptsemObject$computeAndFitKalman()
  }

  if(regCtsemObject$setup$objective == "ML"){
    regCtsemObject$setup$cpptsemObject$computeRAM()
    regCtsemObject$setup$cpptsemObject$fitRAM()
  }

  return(regCtsemObject)
}


#' profileLikelihood
#'
#' computes the profile likelihood for a selected parameter
#'
#' @param cpptsemObject model of type cpptsem
#' @param parameter label of the parameter for which the profile likelihood should be computed
#' @param values values for the parameter at which the profile likelihood will be computed
#' @param lambda penalty value
#' @param startingValues starting values
#' @export
profileLikelihood <- function(cpptsemObject,
                              parameter,
                              values,
                              lambda,
                              startingValues = cpptsemObject$getParameterValues(),
                              adaptiveLassoWeights,
                              N,
                              regIndicators,
                              targetVector,
                              epsilon = 1e-4,
                              controlRsolnp = regCtsem::controlRsolnp()
){
  if(lambda > 0){
    warning("Currently only using approximate regularization and lasso type regularization!")
  }
  # set objective function
  if(any(class(cpptsemObject) == "Rcpp_cpptsemRAMmodel")){
    objective <- "ML"
  }else if (any(class(cpptsemObject) ==  "Rcpp_cpptsemKalmanModel")){
    objective <- "Kalman"
  }else{
    stop("Object has to be of type Rcpp_cpptsemRAMmodel or Rcpp_cpptsemKalmanModel")
  }

  # set optimizer
  if("failureReturns" %in% names(controlRsolnp)){
    failureReturns <- controlRsolnp$failureReturns
  }else{
    warning("No failureReturns for optimx. Using .Machine$double.xmax/2")
    failureReturns <- .Machine$double.xmax/2
  }
  if("eqfun" %in% names(controlRsolnp)){
    eqfun <- controlRsolnp$eqfun
  }else{
    eqfun <- NULL
  }
  if("eqB" %in% names(controlRsolnp)){
    eqB <- controlRsolnp$eqB
  }else{
    eqB <- NULL
  }
  if("ineqfun" %in% names(controlRsolnp)){
    ineqfun <- controlRsolnp$ineqfun
  }else{
    ineqfun <- NULL
  }
  if("ineqLB" %in% names(controlRsolnp)){
    ineqLB <- controlRsolnp$ineqLB
  }else{
    ineqLB <- NULL
  }
  if("ineqUB" %in% names(controlRsolnp)){
    ineqUB <- controlRsolnp$ineqUB
  }else{
    ineqUB <- NULL
  }
  if("LB" %in% names(controlRsolnp)){
    LB <- controlRsolnp$LB
  }else{
    LB <- NULL
  }
  if("UB" %in% names(controlRsolnp)){
    UB <- controlRsolnp$UB
  }else{
    UB <- NULL
  }
  if("control" %in% names(controlRsolnp)){
    control <- controlRsolnp$control
  }else{
    "control" = list("trace" = 0)
  }

  # set free parameters
  free <- rep(TRUE, length(startingValues))
  names(free) <- names(startingValues)
  free[parameter] <- FALSE

  fitFunction <- function(parameters, cpptsemObject, free, adaptiveLassoWeights, N, lambda_, regIndicators, targetVector, epsilon, objective, failureReturns){
    currentValues <- cpptsemObject$getParameterValues()
    # change free parameter only
    currentValues[free] <- parameters[free][names(currentValues[free])]

    # call actual fit function
    regM2LL <- ifelse(tolower(objective) == "ml",
                      regCtsem::approx_RAMRegM2LLCpptsem(parameters = currentValues,
                                                         cpptsemmodel = cpptsemObject,
                                                         adaptiveLassoWeights = adaptiveLassoWeights,
                                                         N = N,
                                                         lambda_ = lambda_,
                                                         regIndicators = regIndicators,
                                                         targetVector = targetVector,
                                                         epsilon = epsilon,
                                                         objective = objective,
                                                         failureReturns = failureReturns),
                      regCtsem::approx_KalmanRegM2LLCpptsem(parameters = currentValues,
                                                            cpptsemmodel = cpptsemObject,
                                                            adaptiveLassoWeights = adaptiveLassoWeights,
                                                            N = N,
                                                            lambda_ = lambda_,
                                                            regIndicators = regIndicators,
                                                            targetVector = targetVector,
                                                            epsilon = epsilon,
                                                            objective = objective,
                                                            failureReturns = failureReturns))
    return(regM2LL)
  }

  lik <- rep(NA, length(values))
  pb <- txtProgressBar(min = 1, max = length(values))

  for(value in values){
    setTxtProgressBar(pb, which(values == value))
    # set parameter
    startingValues[parameter] <- value
    cpptsemObject$setParameterValues(startingValues, names(startingValues))

    # optimize model
    invisible(capture.output(currentFit <- try(Rsolnp::solnp(par = startingValues,
                                                             fun = fitFunction,
                                                             #gr = gradCpptsem,
                                                             eqfun = eqfun, eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB,
                                                             ineqUB = ineqUB, LB = LB, UB = UB, control = control,
                                                             cpptsemObject = cpptsemObject, free = free,
                                                             adaptiveLassoWeights = adaptiveLassoWeights,
                                                             N = N, lambda_ = lambda, regIndicators = regIndicators, targetVector = targetVector,
                                                             epsilon = epsilon, objective = objective, failureReturns = failureReturns),
                                               silent = TRUE), type = c("output", "message")))
    if(any(class(currentFit) == "try-error")){
      warning(paste0("Optimization for ", parameter, " = ", value, " failed."))
      next}
    if(currentFit$convergence > 0 && !silent){
      warning(paste0("Rsolnp reports convcode = ", currentFit$convergence, " for ", parameter, " = ", value, ". See ?Rsolnp::solnp for more details."))}
    lik[which(value == values)] <- currentFit$values[length(currentFit$values)]
  }

  return(data.frame("parameterValue" = values, "Likelihood" = lik))

}

#' profileLikelihoodFromRegCtsem
#'
#' computes the profile likelihood for a selected parameter based on a regularized model
#'
#' @param regCtsemObject model of type regCtsem
#' @param parameter label of the parameter for which the profile likelihood should be computed
#' @param values values for the parameter at which the profile likelihood will be computed
#' @param lambda penalty value
#' @examples
#' \dontrun{
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
#'                                lambdasAutoLength = 5)
#'
#' # compute profile likelihood for drift_eta2_eta1 and lambda = 0
#' pl1 <- profileLikelihoodFromRegCtsem(regCtsemObject = regModel, parameter = "drift_eta2_eta1",
#'                                      values = seq(-1,1,.1),
#'                                      lambda = 0)
#'
#' # compute profile likelihood for drift_eta2_eta1 and lambda = 10
#' pl2 <- profileLikelihoodFromRegCtsem(regCtsemObject = regModel, parameter = "drift_eta2_eta1",
#'                                      values = seq(-1,1,.1),
#'                                      lambda = 10)
#'  }
#' @export
profileLikelihoodFromRegCtsem <- function(regCtsemObject, parameter, values, lambda){
  if(!(regCtsemObject$setup$penalty == "lasso" || regCtsemObject$setup$penalty == "adaptiveLasso")){
    stop("Only implemented for lasso or adaptive lasso")
  }
  lambdas <- regCtsemObject$setup$lambdas
  closestLambda <- which(abs(lambdas - lambda) == min(abs(lambdas - lambda)))[1]
  startingValues <- regCtsemObject$parameterEstimatesRaw[,closestLambda]
  pl <- profileLikelihood(cpptsemObject = regCtsemObject$setup$cpptsemObject,
                          parameter = parameter,
                          values = values,
                          lambda = lambda,
                          startingValues = startingValues,
                          adaptiveLassoWeights = regCtsemObject$setup$adaptiveLassoWeights,
                          N = ifelse(regCtsemObject$setup$scaleLambdaWithN,
                                     nrow(regCtsemObject$setup$dataset),
                                     1),
                          regIndicators = regCtsemObject$setup$regIndicators,
                          targetVector = regCtsemObject$setup$targetVector
  )
  plot(x = pl$parameterValue, y = pl$Likelihood,
       main = "Profile Likelihood",
       xlab = parameter,
       ylab = ifelse(lambda > 0, "regularized likelihood value", "likelihood value"),
       type = "l")
  return(pl)
}

#' profileLikelihoodFromCpptsem
#'
#' computes the profile likelihood for a selected parameter based on a cpptsem
#'
#' @param regCtsemObject model of type regCtsem
#' @param parameter label of the parameter for which the profile likelihood should be computed
#' @param values values for the parameter at which the profile likelihood will be computed
#' @param lambda penalty value
#' @examples
#' \dontrun{
#' library(regCtsem)
#'
#' CINT = matrix(c(0,0), nrow = 2, ncol = 1)
#'
#' stationary <- c('T0TRAITEFFECT','T0TIPREDEFFECT')
#'
#' ## ctsem model without trait
#' AnomAuthmodel1 <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
#'                           Tpoints = 5, n.latent = 2, n.manifest = 2,
#'                           MANIFESTVAR=diag(0, 2),
#'                           TRAITVAR = NULL,
#'                           CINT = CINT)
#' AnomAuthfit1 <- ctFit(AnomAuth, AnomAuthmodel1, useOptimizer = TRUE, stationary = stationary)
#'
#' ## with cpptsem
#' cpptsemModel <- cpptsemFromCtsem(ctsemModel = AnomAuthfit1, wideData = AnomAuth)
#' pl <- profileLikelihoodFromCpptsem(cpptsemObject = cpptsemModel, parameter = "drift_eta1_eta2", values = seq(-1,1,.2))
#' }
#' @export
profileLikelihoodFromCpptsem <- function(cpptsemObject, parameter, values){
  if(any(class(cpptsemObject) == "Rcpp_cpptsemRAMmodel")){
    objective <- "ML"
  }else if (any(class(cpptsemObject) ==  "Rcpp_cpptsemKalmanModel")){
    objective <- "Kalman"
  }else{
    stop("Object has to be of type Rcpp_cpptsemRAMmodel or Rcpp_cpptsemKalmanModel")
  }

  lambda <- 0
  startingValues <- cpptsemObject$getParameterValues()
  resetTo <- startingValues[]

  adaptiveLassoWeights <- startingValues
  adaptiveLassoWeights[] <- 0
  regIndicators <- c()
  targetVector <- startingValues
  targetVector[] <- 0
  pl <- profileLikelihood(cpptsemObject = cpptsemObject,
                          parameter = parameter,
                          values = values,
                          lambda = lambda,
                          startingValues = startingValues,
                          adaptiveLassoWeights = adaptiveLassoWeights,
                          N = 1,
                          regIndicators = regIndicators,
                          targetVector = targetVector
  )
  plot(x = pl$parameterValue, y = pl$Likelihood,
       main = "Profile Likelihood",
       xlab = parameter,
       ylab = ifelse(lambda > 0, "regularized likelihood value", "likelihood value"),
       type = "l")

  cpptsemObject$setParameterValues(resetTo, names(resetTo))
  if(any(class(cpptsemObject) == "Rcpp_cpptsemRAMmodel")){
    cpptsemObject$computeRAM()
    cpptsemObject$fitRAM()
  }else if (any(class(cpptsemObject) ==  "Rcpp_cpptsemKalmanModel")){
    cpptsemObject$computeAndFitKalman()
  }else{
    stop("Object has to be of type Rcpp_cpptsemRAMmodel or Rcpp_cpptsemKalmanModel")
  }
  points(resetTo[parameter], cpptsemObject$m2LL, col = "blue")

  return(pl)
}


#' checkIdentification
#'
#' Checks the local identification of a model around the provided parameters
#'
#' The function is based on the procedure described in Huang, P.-H. (2020). lslx: Semi-Confirmatory Structural Equation Modeling via Penalized Likelihood. Journal of Statistical Software, 93(7). https://doi.org/10.18637/jss.v093.i07
#' @param model model of type ctsemFit, Rcpp_cpptsemRAMmodel, or MxRAMModel. WARNING: Currently only supported for RAM based models
#' @param pars named vector with parameter values
#' @return TRUE if identified, FALSE otherwise
#' @export
checkIdentification <- function(model, pars){
  if(any(class(model) == "ctsemFit")){
    model <- model$mxobj
    if(OpenMx::imxHasDefinitionVariable(model)){stop("Model has definition variables")}
    model <- OpenMx::omxSetParameters(model, labels = names(pars), values = pars)
    getMoments <- function(model, pars){
      model <- omxSetParameters(model, labels = names(pars), values = pars)
      model <- mxRun(model, useOptimizer = FALSE, silent = TRUE)
      moments <- c(OpenMx::mxGetExpected(model, "means"), OpenMx::mxGetExpected(model, "covariance"))
      return(moments)
    }
  }else if(any(class(model) == "Rcpp_cpptsemRAMmodel")){
    model$setParameterValues(pars, names(pars))
    getMoments <- function(model, pars){
      model$setParameterValues(pars, names(pars))
      invisible(capture.output(o <- try(model$computeRAM(), silent = T), type = "message"))
      if(any(class(o) == "try-error")){return(NA)}
      invisible(capture.output(o <- try(model$fitRAM(), silent = T), type = "message"))
      if(any(class(o) == "try-error")){return(NA)}
      moments <- c(model$expectedMeans, model$expectedCovariance)
      return(moments)
    }
  }else if(any(class(model) == "Rcpp_cpptsemKalmanModel")){
    stop("Currently only implemented for RAM models")
  }else if(any(class(model) == "MxModel")){
    if(OpenMx::imxHasDefinitionVariable(model)){stop("Model has definition variables")}
    model <- OpenMx::omxSetParameters(model, labels = names(pars), values = pars)
    getMoments <- function(model, pars){
      model <- omxSetParameters(model, labels = names(pars), values = pars)
      model <- mxRun(model, useOptimizer = FALSE, silent = TRUE)
      moments <- c(OpenMx::mxGetExpected(model, "means"), OpenMx::mxGetExpected(model, "covariance"))
      return(moments)
    }
  }else{stop("model of unknown class")}

  parsNonZero <- pars != 0

  jac <- numDeriv::jacobian(func = getMoments, x = pars, model = model)

  jacSubset <- jac[,parsNonZero]

  return(min(svd(jacSubset)$d) > 0)
}


#' checkNonConvexity
#'
#' checks the non-convexity of a provided model within the provided parameter bounds numerically.
#'
#' This is a very simple implementation of the concept described in Tamura, K., & Gallagher, M. (2019). Quantitative measure of nonconvexity for black-box continuous functions. Information Sciences, 476, 64–82. https://doi.org/10.1016/j.ins.2018.10.009
#' The basic idea is as follows: The function generates nSample points within the provided parameter bounds. For each combination of these sample points the following condition is evaluated:
#' f(x1) + f(x2) >= 2*f((x1 + x2)/2)  (midpoint convexity, see Equation 2 in Tamura, K., & Gallagher, M., 2019)
#' This condition should evaluate to TRUE for all possible combinations of x1 and x2. If it evaluates to FALSE, this implies that the model fitting function is not convex.
#' WARNING: Even if all evaluations indicate that midpoint convexity is given, this is no proof that the function is convex! However, if a counter-example is found (i.e., at least one evaluation results in FALSE), this implies
#' that the function might be nonconvex. Keep in mind that this is only a numerical check and susceptible to numerical errors
#' @param model model of type ctsemFit, Rcpp_cpptsemRAMmodel, Rcpp_cpptsemKalmanModel, or MxRAMModel
#' @param lowerBound vector of the same length as the number of parameters in the model specifying the lower bound for each parameter
#' @param upperBound vector of the same length as the number of parameters in the model specifying the upper bound for each parameter
#' @param nSample number of sample points between the lower and upper bound
#' @param sampleWith possible are runif (random points between lower and upper bound) and seq (equidistant points between lower and upper bound)
#' @return returns the generated samplePoints, fitAtSamplePoints which gives the -2log-Lokelihood at the sample points, a logical vector indicating if the midpoint non-convexity was given for the comparison of these sample points (isNonConvex) and the difference f(x1) + f(x2) - 2*f((x1 + x2)/2) for each combination of sample points. This difference should be positive for all comparisons!
#' @return
checkNonConvexity <- function(model, lowerBound = NULL, upperBound = NULL, nSample = 50, sampleWith = "runif"){
  if(any(class(model) == "ctsemFit")){
    model <- model$mxobj
    pars <- omxGetParameters(model)
    getFit <- function(mxObject, pars){
      mxObject <- omxSetParameters(mxObject, labels = names(pars), values = pars)
      mxObject <- mxRun(mxObject, useOptimizer = FALSE, silent = TRUE)
      return(mxObject$fitfunction$result[[1]])
    }
  }else if(any(class(model) == "Rcpp_cpptsemRAMmodel")){
    pars <- model$getParameterValues()
    getFit <- function(Rcpp_cpptsemRAMmodel, pars){
      Rcpp_cpptsemRAMmodel$setParameterValues(pars, names(pars))
      invisible(capture.output(o <- try(Rcpp_cpptsemRAMmodel$computeRAM(), silent = T), type = "message"))
      if(any(class(o) == "try-error")){return(NA)}
      invisible(capture.output(o <- try(Rcpp_cpptsemRAMmodel$fitRAM(), silent = T), type = "message"))
      if(any(class(o) == "try-error")){return(NA)}
      return(Rcpp_cpptsemRAMmodel$m2LL)
    }
  }else if(any(class(model) == "Rcpp_cpptsemKalmanModel")){
    pars <- model$getParameterValues()
    getFit <- function(Rcpp_cpptsemKalmanModel, pars){
      Rcpp_cpptsemKalmanModel$setParameterValues(pars, names(pars))
      invisible(capture.output(o <- try(Rcpp_cpptsemKalmanModel$computeAndFitKalman(), silent = T), type = "message"))
      if(any(class(o) == "try-error")){return(NA)}
      return(Rcpp_cpptsemKalmanModel$m2LL)
    }
  }else if(any(class(model) == "MxRAMModel")){
    pars <- omxGetParameters(model)
    getFit <- function(mxObject, pars){
      mxObject <- omxSetParameters(mxObject, labels = names(pars), values = pars)
      mxObject <- mxRun(mxObject, useOptimizer = FALSE, silent = TRUE)
      return(mxObject$fitfunction$result[[1]])
    }
  }else{stop("model of unknown class")}

  if(is.null(lowerBound)){lowerBound <- pars - 3*abs(pars)}
  if(is.null(upperBound)){upperBound <- pars + 3*abs(pars)}

  # sample points between lower and upper bounds
  bounds <- cbind(lowerBound, upperBound)
  if(sampleWith == "runif"){
    samplePoints <- apply(bounds, 1,
                          function(x) runif(n = nSample,
                                            min = x[1],
                                            max = x[2]))
  }else if(sampleWith == "seq"){
    samplePoints <- apply(bounds, 1,
                          function(x) seq(from = x[1],
                                          to = x[2],
                                          length.out = nSample))
  }


  message("Computing fit at sample points: \n")
  fitAtSamplePoints <- c()
  pb <- txtProgressBar(min = 0, max = nSample, initial = 0, char = "=",
                       width = NA, title, label, style = 3, file = "")
  it <- 1
  for(samplePoint in 1:nSample){
    it <- it+1
    setTxtProgressBar(pb, it)
    fitAtSamplePoint <- try(getFit(model, pars = samplePoints[samplePoint,]), silent = T)
    if(any(class(fitAtSamplePoint) == "try-error") || is.na(fitAtSamplePoint)){
      observedNA <- TRUE
      fitAtSamplePoints[samplePoint] <- NA
      next}
    fitAtSamplePoints[samplePoint] <- fitAtSamplePoint
  }
  cat("\n")

  message("Evaluating midpoint convexity for all pairs of sample points: \n")
  isNonConvex <- c()
  elemNames <- c()
  midpointDifference <- c()
  observedNA <- FALSE
  maxChecks <- (sum(!is.na(fitAtSamplePoints))-1)*(sum(!is.na(fitAtSamplePoints)))/2
  pb <- txtProgressBar(min = 0, max = maxChecks, initial = 0, char = "=",
                       width = NA, title, label, style = 3, file = "")
  it <- 1
  for(samplePoint1 in 1:nSample){
    f_1 <- fitAtSamplePoints[samplePoint1]
    if(samplePoint1 == nSample){break}
    if(is.na(f_1)){next}

    for(samplePoint2 in (samplePoint1+1):nSample){
      f_2 <- fitAtSamplePoints[samplePoint2]

      if(is.na(f_2)){next}
      it <- it+1
      setTxtProgressBar(pb, it)

      f_3 <- try(getFit(model, pars = (samplePoints[samplePoint1,] + samplePoints[samplePoint2,])/2), silent = T)
      if(any(class(f_3) == "try-error")){
        observedNA <- TRUE
        next}
      if(anyNA(c(f_1,f_2,f_3))){
        observedNA <- TRUE
        next}
      # check midpoint convexity
      isNonConvex <- c(isNonConvex, !(f_1 + f_2 >= 2*f_3))
      elemNames <- c(elemNames, paste0(samplePoint1, " vs ", samplePoint2))
      midpointDifference <- c(midpointDifference, f_1 + f_2 - 2*f_3)
    }
  }
  cat("\n")
  if(observedNA){warning("Some sample points resulted in errors / NA. Only the results for non-NA evaluations are reported.")}
  names(isNonConvex) <- elemNames
  names(midpointDifference) <- elemNames
  rownames(samplePoints) <- 1:nSample

  message(paste0("Of in total ", length(isNonConvex), " non-NA checks ", sum(isNonConvex), " resulted in non-convex evaluations at the mid point."))
  return(list("samplePoints" = samplePoints,"fitAtSamplePoints" = fitAtSamplePoints ,"isNonConvex" = isNonConvex, "midpointDifference" = midpointDifference))
}

#' checkNonConvexity3D
#'
#' provides a 3D plot of the likelihood surface for two selected variables
#'
#' The function generates nSamples points within the provided parameter bounds. For each combination of sample point
#' @param model model of type ctsemFit, Rcpp_cpptsemRAMmodel, Rcpp_cpptsemKalmanModel, or MxRAMModel
#' @param parnames vector of length 2 with the names of the parameters for which the likelihood surface should be plotted
#' @param lowerBound1 double: lower bound for the first parameter
#' @param lowerBound2 double: lower bound for the second parameter
#' @param upperBound1 double: upper bound for the first parameter
#' @param upperBound2 double: upper bound for the second parameter
#' @param nSamples number of sample points between the lower and upper bound
#' @return returns the a list with arguments to pass to plot3D::persp (e.g. do.call(args, plot3D::persp))
#' @return
checkNonConvexity3D <- function(model, parnames, lowerBound1, upperBound1, lowerBound2, upperBound2, nSamples){
  warning("Very experimental function. Don't use...")
  if(any(class(model) == "ctsemFit")){
    model <- model$mxobj
    pars <- omxGetParameters(model)
    getFit <- function(mxObject, pars){
      mxObject <- omxSetParameters(mxObject, labels = names(pars), values = pars)
      mxObject <- mxRun(mxObject, useOptimizer = FALSE, silent = TRUE)
      return(mxObject$fitfunction$result[[1]])
    }
  }else if(any(class(model) == "Rcpp_cpptsemRAMmodel")){
    pars <- model$getParameterValues()
    getFit <- function(Rcpp_cpptsemRAMmodel, pars){
      Rcpp_cpptsemRAMmodel$setParameterValues(pars, names(pars))
      invisible(capture.output(o <- try(Rcpp_cpptsemRAMmodel$computeRAM(), silent = T), type = "message"))
      if(any(class(o) == "try-error")){return(NA)}
      invisible(capture.output(o <- try(Rcpp_cpptsemRAMmodel$fitRAM(), silent = T), type = "message"))
      if(any(class(o) == "try-error")){return(NA)}
      return(Rcpp_cpptsemRAMmodel$m2LL)
    }
  }else if(any(class(model) == "Rcpp_cpptsemKalmanModel")){
    pars <- model$getParameterValues()
    getFit <- function(Rcpp_cpptsemKalmanModel, pars){
      Rcpp_cpptsemKalmanModel$setParameterValues(pars, names(pars))
      invisible(capture.output(o <- try(Rcpp_cpptsemKalmanModel$computeAndFitKalman(), silent = T), type = "message"))
      if(any(class(o) == "try-error")){return(NA)}
      return(Rcpp_cpptsemKalmanModel$m2LL)
    }
  }else if(any(class(model) == "MxRAMModel")){
    pars <- omxGetParameters(model)
    getFit <- function(mxObject, pars){
      mxObject <- omxSetParameters(mxObject, labels = names(pars), values = pars)
      mxObject <- mxRun(mxObject, useOptimizer = FALSE, silent = TRUE)
      return(mxObject$fitfunction$result[[1]])
    }
  }else{stop("model of unknown class")}

  if(!any(parnames %in% names(pars))){stop("parnames not found in model")}

  pars1 <- seq(lowerBound1, upperBound1, length.out = nSamples)
  pars2 <- seq(lowerBound2, upperBound2, length.out = nSamples)

  fitMatrix <- matrix(NA, nrow = length(pars1), ncol = length(pars2))
  rownames(fitMatrix) <- paste0(parnames[1], " = ", pars1)
  colnames(fitMatrix) <- paste0(parnames[2], " = ", pars2)

  pb <- txtProgressBar(min = 0, max = nSamples^2, initial = 0, char = "=",
                       width = NA, title, label, style = 3, file = "")
  it <- 1
  for(par1 in 1:length(pars1)){
    for(par2 in 1:length(pars2)){
      it <- it+1
      setTxtProgressBar(pb, it)
      pars[parnames[1]] <- pars1[par1]
      pars[parnames[2]] <- pars2[par2]
      fitMatrix[par1, par2] <- getFit(model, pars)
    }
  }
  fitMatrix[!is.finite(fitMatrix)] <- NA
  plotArgs <- list(x = pars1,
                   y = pars2,
                   z = fitMatrix, xlab = parnames[1], ylab = parnames[2],
                   theta = 30,
                   col = "springgreen", shade = 0.5)
  try(do.call(persp, plot3D::plotArgs))

  return(plotArgs)
}


#' shinyfy
#'
#' creates a shiny application from a regCtsemObject
#' @param regCtsemObject regCtsemObject
#' @param runShiny boolean: if TRUE, the app is started directly. If false, the code to generate the app is returned
#' @param figHeight height of the figure
#' @export
shinyfy <- function(regCtsemObject, runShiny = TRUE, figHeight = 600){
  if(regCtsemObject$setup$autoCV){stop("Shinyfy currently not supported for models with automatic cross-validation.")}
  filename <- tempfile(pattern = "modelClone", tmpdir = tempdir(), fileext = ".RData")
  save(regCtsemObject, file = filename)
  lambdas <- regCtsemObject$setup$lambdas
  criteria <- c("m2LL","AIC", "BIC", "cvM2LL")
  criteria <- criteria[criteria %in% rownames(regCtsemObject$fit)]

  selCriterium <- paste0("list(", paste0(paste0("'", criteria , "' = '" , criteria, "'"), collapse = ", "), ")")
  selLambda <- paste0("list(", paste0(paste0("'", lambdas , "' = " , seq_len(length(lambdas))), collapse = ", "), ")")

  shinyCode <- paste0(
    'library(shiny)
library(ggplot2)
library(gridExtra)
library(regCtsem)
library(qgraph)
load("', filename, '")

latentNames <- regCtsemObject$setup$ctsemObject$ctmodelobj$latentNames
lambdas <- regCtsemObject$setup$lambdas
fit <- regCtsemObject$fit
parameters <- regCtsemObject$parameterEstimatesRaw
parameterNames <- rownames(parameters)
DRIFTlabels <- regCtsemObject$setup$mxObject$DRIFT$labels
DRIFTvalues <- regCtsemObject$setup$mxObject$DRIFT$values
DIFFUSIONbaseLabels <- regCtsemObject$setup$mxObject$DIFFUSIONbase$labels
DIFFUSIONbaseValues <- regCtsemObject$setup$mxObject$DIFFUSIONbase$values

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("regCtsem model inspection"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput("criterion", h3("criterion"), ', selCriterium,'
                        ),
            checkboxInput("discrete", "discrete", value = FALSE),
            numericInput("dT",
                         h3("time interval"),
                         value = 1),
            checkboxInput("manualLambda", "manually select lambda", value = FALSE),
            selectInput("lambdaVal", h3("lambda"),', selLambda,', selected = 1)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("qplot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$qplot <- renderPlot({
        makeEdges <- function(mat){
            from <- 1:ncol(mat)
            to <- 1:nrow(mat)
            weights <- as.vector(mat)
            edges <- matrix(nrow = length(weights), ncol = 3)
            colnames(edges) <- c("from", "to", "weight")
            edges[,"weight"] <- weights
            edges[,"from"] <- rep(from, each = length(to))
            edges[,"to"] <- rep(to, length(from))
            return(edges)
        }

        manualLambda <- input$manualLambda
        dT <- as.numeric(input$dT)
        criterion <- input$criterion
        lambdaVal <- as.integer(input$lambdaVal)
        discrete <- input$discrete

        ## check if manualLambda was set
        if(manualLambda){
            lambda <- lambdas[lambdaVal]
        }else{
            # else use criterion
            finalPars <- getFinalParameters(regCtsemObject, criterion = criterion, raw = TRUE)
            lambda <- finalPars$lambda
        }

        # set values
        for(DRIFTlabel in unique(c(DRIFTlabels[!is.na(DRIFTlabels)]))){
            DRIFTvalues[DRIFTlabels == DRIFTlabel & !is.na(DRIFTlabels)] <- parameters[DRIFTlabel,lambdas == lambda]
        }
        DRIFTHASH <- regCtsem::computeDRIFTHASH(DRIFTvalues)

        for(DIFFUSIONbaseLabel in unique(c(DIFFUSIONbaseLabels[!is.na(DIFFUSIONbaseLabels)]))){
            DIFFUSIONbaseValues[DIFFUSIONbaseLabels == DIFFUSIONbaseLabel & !is.na(DIFFUSIONbaseLabels)] <- parameters[DIFFUSIONbaseLabel,lambdas == lambda]
        }
        DIFFUSION <- regCtsem::getVarianceFromVarianceBase2(varianceBaseValues = DIFFUSIONbaseValues)

        if(discrete){

        discreteDRIFT <- expm(DRIFTvalues*dT)
        discreteDiffusion <- matrix(solve(DRIFTHASH) %*% (expm(DRIFTHASH*dT) - diag(nrow(DRIFTHASH)))%*%c(DIFFUSION),
                                    nrow = nrow(DIFFUSION), ncol = ncol(DIFFUSION))

        layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), heights = c(2,1))
        edges <- makeEdges(mat = discreteDRIFT)
        edge.color <- ifelse(edges[,"weight"]>0, "#008080", "#800000")
        qgraph(title = "discrete autoregressive and cross-lagged parameters",
               edges,
               maximum = .6,
               fade=FALSE,
               layout="circular",
               labels = latentNames,
               lty=ifelse(edges[,"weight"]>0,1,5),
               edge.labels=T, font = 2,
               edge.color=edge.color)

        edges <- makeEdges(mat = discreteDiffusion)
        edge.color <- ifelse(edges[,"weight"]>0, "#008080", "#800000")
        qgraph(title = "discrete diffusion",
               edges,
               maximum = .6,
               fade=FALSE,
               bidirectional = TRUE,
               layout="circular",
               labels = latentNames,
               lty=ifelse(edges[,"weight"]>0,1,5),
               edge.labels=T, font = 2,
               edge.color=edge.color)
        plot(regCtsemObject = regCtsemObject, what = "fit", criterion = criterion)
        }else{
        edges <- makeEdges(mat = DRIFTvalues)

        edge.color <- ifelse(edges[,"weight"]>0, "#008080", "#800000")
        layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), heights = c(2,1))
        qgraph(title = "DRIFT",
               edges,
               maximum = .6,
               fade=FALSE,
               layout="circular",
               labels = latentNames,
               lty=ifelse(edges[,"weight"]>0,1,5),
               edge.labels=T, font = 2,
               edge.color=edge.color)

        edges <- makeEdges(mat = DIFFUSION)
        edge.color <- ifelse(edges[,"weight"]>0, "#008080", "#800000")
        qgraph(title = "DIFFUSION",
               edges,
               maximum = .6,
               fade=FALSE,
               bidirectional = TRUE,
               layout="circular",
               labels = latentNames,
               lty=ifelse(edges[,"weight"]>0,1,5),
               edge.labels=T, font = 2,
               edge.color=edge.color)
        plot(regCtsemObject = regCtsemObject, what = "fit", criterion = criterion)

        }
    }, height = ', figHeight, ')
  }

  # Run the application
shinyApp(ui = ui, server = server)
  ')

if(runShiny){
  filename <- tempfile(pattern = "regCtsem_Shiny", tmpdir = tempdir(), fileext = ".R")
  fileConn <- file(filename)
  writeLines(shinyCode, fileConn)
  close(fileConn)

  shiny::runApp(filename)

}else{
  return(shinyCode)
}
}
