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
#' Importantly, the untransformed parameters are returned.
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
  return(extractParameterTableFromMx(ctsemObject$mxobj))
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
#' @param lambdasAutoLength lambdasAutoLength will determine the number of lambdas tested.
#' @param lambdasAutoCurve It is often a good idea to have unequally spaced lambda steps (e.g., .01,.02,.05,1,5,20). If lambdasAutoCurve is close to 1 lambda values will be equally spaced, if lambdasAutoCurve is large lambda values will be more concentrated close to 0. See ?getCurvedLambda for more informations.
#' @param penalty Currently supported are lasso, ridge and adaptiveLasso
#' @param adaptiveLassoWeights weights for the adaptive lasso. Defaults to 1/(|theta|^adaptiveLassoPower), where theta is the maximum likelihood estimate of the regularized parameters.
#' @param adaptiveLassoPower power for the adaptive lasso weights. The weights will be set to 1/(|theta|^adaptiveLassoPower).
#' @param cvSample cross-validation sample. Has to be in wide format and compatible with ctsemOMX
#' @param autoCV Should automatic cross-validation be used? Possible are "No", "kFold" or "Blocked". kFold splits the dataset in k groups by selecting independent units from the rows. Blocked is a within-unit split, where for each person blocks of observations are deleted. See Bulteel, K., Mestdagh, M., Tuerlinckx, F., & Ceulemans, E. (2018). VAR(1) based models do not always outpredict AR(1) models in typical psychological applications. Psychological Methods, 23(4), 740–756. https://doi.org/10.1037/met0000178 for a more detailed explanation
#' @param k number of cross-validation folds if autoCV = "kFold" or autoCV = "Blocked"
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
    message("When using fitFull = FALSE the standardization and the adaptive LASSO weights are based on a (possibly very poor) approximation of the full model. Be wary of that!")
  }

  if(is.null(targetVector)){
    targetVector <- rep(0, length(regIndicators))
    names(targetVector) <- regIndicators
    targetVector <- targetVector
  }

  objective <- ifelse(ctsemObject$ctfitargs$objective == "Kalman", "Kalman", "ML")

  cpptsemObject <- try(cpptsemFromCtsem(ctsemModel = ctsemObject, wideData = dataset))
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
  grad <- exact_getCppGradients(cpptsemObject, objective = objective)

  if(tolower(penalty) == "lasso" && standardizeDrift != "No"){
    message("Standardizing with parameters of the sparse model")
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
      message("Full model will only be approximated!")
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
      approximatedFull <- stats::optim(par = param,
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

  lambdas <- rev(getCurvedLambda(maxLambda = maxLambda + maxLambda/25, # adding some wiggle room as there will always be some deviations
                                           lambdasAutoCurve = lambdasAutoCurve,
                                           lambdasAutoLength = lambdasAutoLength))

  # fit ctsem model to have correct starting values:
  cpptsemObject$setParameterValues(sparseParameters, names(sparseParameters))
  ctsemObject$mxobj <- omxSetParameters(ctsemObject$mxobj,
                                        labels = names(cpptsemObject$getParameterValues()),
                                        values = cpptsemObject$getParameterValues())
  ctsemObject$mxobj <- mxRun(ctsemObject$mxobj, useOptimizer = FALSE)

  if(optimizeRegCtsem)  {
    tryFit <- try(regCtsem(ctsemObject = ctsemObject,
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

