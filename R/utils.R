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
    regCtsemObject$setup$cpptsemObject$computeAndFitKalman(0)
  }

  if(regCtsemObject$setup$objective == "ML"){
    regCtsemObject$setup$cpptsemObject$computeRAM()
    regCtsemObject$setup$cpptsemObject$fitRAM()
  }

  return(regCtsemObject)
}

#' showParameters
#'
#' shows the parameters of a model fitted with ctFit.
#' Importantly, the untransformed paraemters are returned.
#' For instance, the variance-covariance matrices are implemented in log-Cholesky form
#' (see Pinheiro, J. C., & Bates, D. M. (1996). Unconstrained parametrizations for variance-covariance matrices. Statistics and Computing, 6(3), 289–296. https://doi.org/10.1007/BF00140873).
#'
#' @param ctsemObject Fitted object of class ctsemFit
#' @examples
#' library(regCtsem)
#'
#' # The following example is taken from ?ctFit:
#' data(AnomAuth)
#' AnomAuthmodel <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
#'                          Tpoints = 5, n.latent = 2, n.manifest = 2, MANIFESTVAR=diag(0, 2), TRAITVAR = NULL)
#' AnomAuthfit <- ctFit(AnomAuth, AnomAuthmodel)
#'
#' showParameters(AnomAuthfit)
#' @export
showParameters <- function(ctsemObject){
  return(regCtsem::extractParameterTableFromMx(ctsemObject$mxobj))
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
  startingValues2 <- startingValues[] # startingValues2 will be used as an alternative set of starting values based on the previous iteration
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

    if(!is.null(startingValues2)){
      startingValues2[parameter] <- value
      cpptsemObject$setParameterValues(startingValues2, names(startingValues2))

      # optimize model
      invisible(capture.output(currentFit2 <- try(Rsolnp::solnp(par = startingValues2,
                                                                fun = fitFunction,
                                                                #gr = gradCpptsem,
                                                                eqfun = eqfun, eqB = eqB, ineqfun = ineqfun, ineqLB = ineqLB,
                                                                ineqUB = ineqUB, LB = LB, UB = UB, control = control,
                                                                cpptsemObject = cpptsemObject, free = free,
                                                                adaptiveLassoWeights = adaptiveLassoWeights,
                                                                N = N, lambda_ = lambda, regIndicators = regIndicators, targetVector = targetVector,
                                                                epsilon = epsilon, objective = objective, failureReturns = failureReturns),
                                                  silent = TRUE), type = c("output", "message")))

    }
    if(any(class(currentFit) == "try-error") && any(class(currentFit2) == "try-error")){
      warning(paste0("Optimization for ", parameter, " = ", value, " failed."))
      next}
    if((any(class(currentFit) == "try-error" || currentFit$convergence > 0) && (any(class(currentFit2) == "try-error") || currentFit2$convergence > 0))){
      warning(paste0("Rsolnp resulted in errors for ", parameter, " = ", value, ". Skipping this value."))
      next
    }
    if(any(class(currentFit) == "try-error") || currentFit$convergence > 0 || currentFit$values[length(currentFit$values)] > currentFit2$values[length(currentFit2$values)]){
      startingValues2 <- currentFit2$pars
      lik[which(value == values)] <- currentFit2$values[length(currentFit2$values)]
      next
    }
    startingValues2 <- currentFit$pars
    lik[which(value == values)] <- currentFit$values[length(currentFit$values)]
  }

  return(data.frame("parameterValue" = values, "Likelihood" = lik))

}

#' profileLikelihoodFromMx
#'
#' computes the profile likelihood for a selected parameter based on an mxModel
#'
#' @param mxModelObject model from OpenMx
#' @param parameter label of the parameter for which the profile likelihood should be computed
#' @param values values for the parameter at which the profile likelihood will be computed
#' @param startingValues starting values
#' @export
profileLikelihoodFromMx <- function(mxModelObject,
                                    parameter,
                                    values,
                                    startingValues = omxGetParameters(mxModelObject),
                                    tryHard = FALSE
){

  lik <- rep(NA, length(values))
  startingValues2 <- startingValues[] # startingValues2 will be used as an alternative set of starting values based on the previous iteration
  free <- rep(TRUE, length(startingValues))
  free[names(startingValues) == parameter] <- FALSE
  for(value in values){
    cat(paste0(which(values == value), " of ", length(values)), "\n")
    # set parameter
    startingValues[parameter] <- value
    mxModelObject1 <- omxSetParameters(mxModelObject, values = startingValues, labels = names(startingValues), free = free)

    # optimize model
    invisible(capture.output(mxModelObjectFit <- try(mxTryHardctsem(mxModelObject1, extraTries = ifelse(tryHard,10,0), silent = TRUE),
                                                     silent = TRUE), type = c("output", "message")))

    if(!is.null(startingValues2)){
      startingValues2[parameter] <- value
      mxModelObject2 <- omxSetParameters(mxModelObject, values = startingValues, labels = names(startingValues), free = free)

      # optimize model
      invisible(capture.output(mxModelObjectFit2 <- try(mxTryHardctsem(mxModelObject2, extraTries = ifelse(tryHard,10,0), silent = TRUE),
                                                        silent = TRUE), type = c("output", "message")))

    }
    if(any(class(mxModelObjectFit) == "try-error") && any(class(mxModelObjectFit2) == "try-error")){
      warning(paste0("Optimization for ", parameter, " = ", value, " failed."))
      next}
    if((any(class(mxModelObjectFit) == "try-error" || mxModelObjectFit$output$status$code > 0) && (any(class(mxModelObjectFit2) == "try-error") || mxModelObjectFit2$output$status$code > 0))){
      warning(paste0("Optimization resulted in errors for ", parameter, " = ", value, ". Skipping this value."))
      next
    }
    if(any(class(mxModelObjectFit) == "try-error") || mxModelObjectFit$output$status$code > 0 || mxModelObjectFit$fitfunction$result[[1]] > mxModelObjectFit2$fitfunction$result[[1]]){
      startingValues2 <- omxGetParameters(mxModelObjectFit2)
      lik[which(value == values)] <- mxModelObjectFit2$fitfunction$result[[1]]
      next
    }
    startingValues2 <- omxGetParameters(mxModelObjectFit)
    lik[which(value == values)] <- mxModelObjectFit$fitfunction$result[[1]]
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
      invisible(capture.output(o <- try(Rcpp_cpptsemKalmanModel$computeAndFitKalman(0), silent = T), type = "message"))
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
      invisible(capture.output(o <- try(Rcpp_cpptsemKalmanModel$computeAndFitKalman(0), silent = T), type = "message"))
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


#' startFromSparse
#'
#' alpha status function. Start regularized model from a setting where all parameters are already at their target values
#'
#' NOTE: Function located in file utils.R
#'
#' @param ctsemObject Fitted object of class ctsemFit
#' @param dataset Data set in wide format compatible with ctsemOMX
#' @param regIndicators Labels of the regularized parameters (e.g. drift_eta1_eta2).
#' @param targetVector named vector with values towards which the parameters are regularized (Standard is regularization towards zero)
#' @param lambdas vector of penalty values (tuning parameter). E.g., seq(0,1,.01). Alternatively, lambdas can be set to "auto". regCtsem will then compute an upper limit for lambda and test lambdasAutoLength increasing lambda values
#' @param lambdasAutoLength if lambdas == "auto", lambdasAutoLength will determine the number of lambdas tested.
#' @param lambdasAutoCurve It is often a good idea to have unequally spaced lambda steps (e.g., .01,.02,.05,1,5,20). If lambdasAutoCurve is close to 1 lambda values will be equally spaced, if lambdasAutoCurve is large lambda values will be more concentrated close to 0. See ?getCurvedLambda for more informations.
#' @param penalty Currently supported are lasso, ridge and adaptiveLasso
#' @param adaptiveLassoWeights weights for the adaptive lasso. Defaults to 1/(|theta|^adaptiveLassoPower), where theta is the maximum likelihood estimate of the regularized parameters.
#' @param adaptiveLassoPower power for the adaptive lasso weights. The weights will be set to 1/(|theta|^adaptiveLassoPower).
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
#' @param nMultistart number of additional tries when optimizing the models
#' @param fitFull boolean: Only used for adaptiveLASSO weights. Should the full model (model without regularization) be fitted (TRUE) or only approximated (FALSE). Approximation might work sometimes, but not always.
#' @param optimizeRegCtsem if set to false, the function will only return the lambda_max, a vector with sparse parameter values and a vector for the full, unregularized model parameters
#' @examples
#'
#' @author Jannik Orzek
#' @import numDeriv
#' @export
startFromSparse <- function(ctsemObject,
                            dataset,
                            regIndicators,
                            targetVector = NULL,
                            lambdasAutoLength = 50,
                            lambdasAutoCurve = 10,
                            penalty = "lasso",
                            adaptiveLassoWeights = NULL,
                            adaptiveLassoPower = -1,

                            cvSample = NULL,
                            autoCV = "No",
                            k = 5,
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
                            trainingWheels = TRUE,

                            nMultistart = 3,
                            fitFull = TRUE,
                            optimizeRegCtsem = TRUE){
  if(autoCV!= "No") stop("Cross-Validation not yet implemented")
  if(!is.null(subjectSpecificParameters)) stop("Subject specific parameters not yet implemented")
  if(!fitFull & (penalty == "adaptiveLASSO" | standardizeDrift != "No")){
    warning("When using fitFull = FALSE the standardization and the adaptive LASSO weights are based on a (possibly very poor) approximation of the full model. Be wary of that!")
  }

  if(is.null(targetVector)){
    targetVector <- rep(0, length(regIndicators))
    names(targetVector) <- regIndicators
    targetVector <- targetVector
  }

  objective <- ifelse(ctsemObject$ctfitargs$objective == "Kalman", "Kalman", "ML")

  cpptsemObject <- try(regCtsem::cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = dataset))
  parameters <- cpptsemObject$getParameterValues()

  # change parameters to target values:
  sparseParameters <- parameters

  # step 1: set the regularized parameters to their respective target and estimate the model:

  sparseParameters[regIndicators] <- targetVector[regIndicators]
  freeParam <- rep(TRUE, length(sparseParameters))
  names(freeParam) <- names(sparseParameters)
  freeParam[regIndicators] <- FALSE

  cpptsemObject$setParameterValues(sparseParameters, names(sparseParameters))

  # check if start values are feasible
  tryFit <- fitCpptsem(parameterValues = sparseParameters[freeParam],
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
          sparseParameters[paramterTable$label[p]] <- -.45
        }else{
          sparseParameters[paramterTable$label[p]] <- -.05
        }
      }
      if(paramterTable$matrix[p] == "DIFFUSIONbase"){
        if(paramterTable$row[p] == paramterTable$col[p]){
          sparseParameters[paramterTable$label[p]] <- log(10)
        }else{
          sparseParameters[paramterTable$label[p]] <- 0
        }
      }
      if(paramterTable$matrix[p] == "T0VARbase"){
        if(paramterTable$row[p] == paramterTable$col[p]){
          sparseParameters[paramterTable$label[p]] <- log(1)
        }else{
          sparseParameters[paramterTable$label[p]] <- 0
        }
      }
      if(paramterTable$matrix[p] == "MANIFESTMEANS"){
        sparseParameters[paramterTable$label[p]] <- 0
      }
    }
    sparseParameters[regIndicators] <- targetVector[regIndicators]
    tryFit <- fitCpptsem(parameterValues = sparseParameters[freeParam],
                         cpptsemObject = cpptsemObject,
                         objective = objective,
                         free = freeParam,
                         failureReturns = NA)
    if(is.na(tryFit)){
      stop("Impossible to start model with tested starting values. Try setting the starting values manually.")
    }
  }

  # optimize
  sparseModel <- try(optimizeCpptsem(cpptsemObject = cpptsemObject,
                                     free = freeParam,
                                     nMultistart = nMultistart),
                     silent = TRUE)

  if(any(class(sparseModel) == "try-error")){
    stop("Error when optimizing sparse model. Try different starting values.")
  }

  nonZeroParam <- sparseModel$pars
  namesNonZeroParam <- names(nonZeroParam)

  # step 2: compute gradients with regularized parameters set to target and unregularized parameters set to nonZeroParam estimates
  sparseParameters[namesNonZeroParam] <- nonZeroParam
  cpptsemObject$setParameterValues(sparseParameters, names(sparseParameters))
  grad <- regCtsem::exact_getCppGradients(cpptsemObject, objective = objective)

  if(tolower(penalty) == "lasso" && standardizeDrift != "No"){
    warning("Standardizing with parameters of the sparse model")
    cpptsemObject$setParameterValues(sparseParameters, names(sparseParameters))
    tryFit <- fitCpptsem(parameterValues = sparseParameters[freeParam],
                         cpptsemObject = cpptsemObject,
                         objective = objective,
                         free = freeParam,
                         failureReturns = NA)
  }
  if(tolower(penalty) == "adaptivelasso"){
    # Problem: We need the full model to define the adaptive LASSO weights
    if(fitFull){
      fullModel <- try(optimizeCpptsem(cpptsemObject = cpptsemObject,
                                       nMultistart = nMultistart),
                       silent = TRUE)
      fullParameters <- fullModel$par
    }else{
      warning("Full model will only be approximated!")
      # step 3: compute hessian
      hess <- numDeriv::hessian(func = fitCpptsem,
                                x = cpptsemObject$getParameterValues(),
                                cpptsemObject = cpptsemObject,
                                objective = "ML",
                                failureReturns = .Machine$double.xmax/2)

      fn <- function(x, sparseParameters, grad, hess){
        diff <- matrix((x - sparseParameters), nrow = 1)
        return(diff%*%matrix(grad, ncol = 1) + .5*diff%*%hess%*%t(diff))
      }
      param <- sparseParameters
      approximatedFull <- optim(par = param,
                                fn = fn,
                                sparseParameters = param,
                                grad = grad,
                                hess = hess)
      fullParameters <- approximatedFull$par
    }

    cpptsemObject$setParameterValues(fullParameters, names(fullParameters))
  }
  tryFit <- fitCpptsem(parameterValues = cpptsemObject$getParameterValues(),
                       cpptsemObject = cpptsemObject,
                       objective = objective,
                       failureReturns = NA)
  adaptiveLassoWeights <- getAdaptiveLassoWeights(cpptsemObject = cpptsemObject,
                                                  penalty = penalty,
                                                  adaptiveLassoWeights = adaptiveLassoWeights,
                                                  adaptiveLassoPower =  adaptiveLassoPower,
                                                  standardizeDrift = standardizeDrift)
  # define maxLambda as the maximal gradient of the regularized parameters
  maxLambda <- max(abs(grad[regIndicators]) * adaptiveLassoWeights[regIndicators]^(-1))/ifelse(scaleLambdaWithN, nrow(dataset),1)

  lambdas <- rev(regCtsem::getCurvedLambda(maxLambda = maxLambda + maxLambda/25, # adding some wiggle room as there will always be some deviations
                                           lambdasAutoCurve = lambdasAutoCurve,
                                           lambdasAutoLength = lambdasAutoLength))

  # fit ctsem model to have correct starting values:
  cpptsemObject$setParameterValues(sparseParameters, names(sparseParameters))
  ctsemObject$mxobj <- omxSetParameters(ctsemObject$mxobj,
                                        labels = names(cpptsemObject$getParameterValues()),
                                        values = cpptsemObject$getParameterValues())
  ctsemObject$mxobj <- mxRun(ctsemObject$mxobj, useOptimizer = FALSE)

  if(optimizeRegCtsem)  {
    tryFit <- try(regCtsem::regCtsem(ctsemObject = ctsemObject,
                                     dataset = dataset,
                                     regIndicators = regIndicators,
                                     targetVector = targetVector,
                                     lambdas = lambdas,
                                     lambdasAutoLength = lambdasAutoLength,
                                     lambdasAutoCurve = lambdasAutoCurve,
                                     penalty = penalty,
                                     adaptiveLassoWeights = adaptiveLassoWeights,
                                     adaptiveLassoPower = adaptiveLassoPower,
                                     cvSample = cvSample,
                                     autoCV = autoCV,
                                     k = k,
                                     subjectSpecificParameters = subjectSpecificParameters,
                                     standardizeDrift = "No",
                                     scaleLambdaWithN = scaleLambdaWithN,
                                     returnFitIndices = returnFitIndices,
                                     BICWithNAndT = BICWithNAndT,
                                     optimization = optimization,
                                     optimizer = optimizer,
                                     control = control,
                                     verbose = verbose,
                                     trainingWheels = trainingWheels)
    )
    if(any(class(tryFit) == "try-error")) {warning("Error while fitting regularized model.")}else{return(tryFit)}
  }

  return(list("maxLambda" = maxLambda, "sparseParameters" = sparseParameters, "fullParameters" = cpptsemObject$getParameterValues()))

}

approximateLOOCV <- function(cpptsemObject, wideData, Hessian = NULL, eps = 1e-4){
  warning("WORK IN PROGRESS")
  parameterValues <- cpptsemObject$getParameterValues()

  if(class(cpptsemObject) == "Rcpp_cpptsemRAMmodel") {
    objective <- "ML"
  }else if(class(cpptsemObject) == "Rcpp_cpptsemKalmanModel"){
    objective <- "Kalman"
  }else{stop("Model of unknown class passed to approximateLOOCV")}

  # step 1: compute Hessian (if not provided)
  if(is.null(Hessian)){
    Hessian <- optimHess(par = parameterValues,
                         fn = regCtsem::fitCpptsem,
                         cpptsemObject = cpptsemObject,
                         objective = objective,
                         failureReturns = .5*.Machine$double.xmax)
  }
  if(any(eigen(Hessian)$values < 0)) stop("Hessian is not positive definite.")

  # step 2: compute scores (individual gradients)
  if(objective == "ML"){
    scores <- computeScores.Rcpp_cpptsemRAMmodel(cpptsemObject = cpptsemObject, wideData = wideData, eps = eps)
  }
  if(objective == "Kalman"){
    scores <- computeScores.Rcpp_cpptsemKalmanModel(cpptsemObject = cpptsemObject, wideData = wideData, eps = eps)
  }

  # approximate parameters for Leave One Out Training sets
  trainingParameters <- matrix(NA, nrow = length(parameterValues), ncol = nrow(wideData))
  rownames(trainingParameters) <- names(parameterValues)

  invHessian <- solve(Hessian)

  for(i in 1:nrow(wideData)){
    trainingParameters[,i] <- parameterValues + (1/(nrow(wideData)-1))*invHessian%*%matrix(scores[names(parameterValues),i], ncol = 1)
  }
}


individualMinus2LogLikelihoods.Rcpp_cpptsemRAMmodel <- function(cpptsemObject, wideData){
  dataset <- subset(wideData, select = !grepl("dT", colnames(wideData)) & !grepl("intervalID", colnames(wideData)) )

  expectedMeans <- cpptsemObject$expectedMeans
  expectedCovariance <- cpptsemObject$expectedCovariance

  inidividualLikelihoods <- rep(0, nrow(wideData))
  for(i in 1:nrow(wideData)){
    missings <- is.na(dataset[i,])

    # check if all missing
    if(all(missings)) next

    expectedMeans_i <- expectedMeans[!missings]
    expectedCovariance_i <- expectedCovariance[!missings , !missings]
    observed_i <- as.matrix(dataset[i,!missings])
    inidividualLikelihoods[i] <- regCtsem:::computeIndividualM2LL(nObservedVariables = length(observed_i),
                                                                  rawData = observed_i,
                                                                  expectedMeans = expectedMeans_i,
                                                                  expectedCovariance = expectedCovariance_i)
  }
  return(inidividualLikelihoods)
}


individualMinus2LogLikelihoods.Rcpp_cpptsemKalmanModel <- function(cpptsemObject, wideData){
  return(cpptsemObject$indM2LL)
}

computeScores.Rcpp_cpptsemRAMmodel <- function(cpptsemObject, wideData, eps = 1e-4){
  parameterValues <- cpptsemObject$getParameterValues()

  m2LLs <- matrix(NA, nrow = nrow(wideData), ncol = 2)
  scores <- matrix(NA, nrow = length(parameterValues), ncol = nrow(wideData))
  rownames(scores) <- names(parameterValues)

  parameterValues_i <- parameterValues

  for(i in 1:length(parameterValues)){

    # step backward
    parameterValues_i[i] <- parameterValues_i[i] - eps
    cpptsemObject$setParameterValues(parameterValues_i, names(parameterValues_i))

    # fit
    cpptsemObject$computeRAM()
    cpptsemObject$fitRAM()

    # get individual fits
    m2LLs[,2] <- individualMinus2LogLikelihoods.Rcpp_cpptsemRAMmodel(cpptsemObject, wideData)

    # step forward
    parameterValues_i[i] <- parameterValues_i[i] + 2*eps
    cpptsemObject$setParameterValues(parameterValues_i, names(parameterValues_i))

    # fit
    cpptsemObject$computeRAM()
    cpptsemObject$fitRAM()

    # get individual fits
    m2LLs[,1] <- individualMinus2LogLikelihoods.Rcpp_cpptsemRAMmodel(cpptsemObject, wideData)

    # reset
    parameterValues_i[i] <- parameterValues_i[i] - eps
    cpptsemObject$setParameterValues(parameterValues_i, names(parameterValues_i))

    # compute gradients
    scores[names(parameterValues_i[i]), ] <- (m2LLs[,1] - m2LLs[,2])/2
  }
  return(scores)
}


computeScores.Rcpp_cpptsemKalmanModel <- function(cpptsemObject, wideData, eps = 1e-4){
  parameterValues <- cpptsemObject$getParameterValues()

  m2LLs <- matrix(NA, nrow = nrow(wideData), ncol = 2)
  scores <- matrix(NA, nrow = length(parameterValues), ncol = nrow(wideData))
  rownames(scores) <- names(parameterValues)

  parameterValues_i <- parameterValues

  for(i in 1:length(parameterValues)){

    # step backward
    parameterValues_i[i] <- parameterValues_i[i] - eps
    cpptsemObject$setParameterValues(parameterValues_i, names(parameterValues_i))

    # fit
    cpptsemObject$computeAndFitKalman(0)

    # get individual fits
    m2LLs[,2] <- individualMinus2LogLikelihoods.Rcpp_cpptsemKalmanModel(cpptsemObject, wideData)

    # step forward
    parameterValues_i[i] <- parameterValues_i[i] + 2*eps
    cpptsemObject$setParameterValues(parameterValues_i, names(parameterValues_i))

    # fit
    cpptsemObject$computeAndFitKalman(0)

    # get individual fits
    m2LLs[,1] <- individualMinus2LogLikelihoods.Rcpp_cpptsemKalmanModel(cpptsemObject, wideData)

    # reset
    parameterValues_i[i] <- parameterValues_i[i] - eps
    cpptsemObject$setParameterValues(parameterValues_i, names(parameterValues_i))

    # compute gradients
    scores[names(parameterValues_i[i]), ] <- (m2LLs[,1] - m2LLs[,2])/2
  }

  return(scores)
}
