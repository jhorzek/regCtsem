#' getMaxRegValue
#'
#' computes an approximation of the lowest regValue which will set all regularized parameters to zero. This function is adapted from Murphy (2012) Machine learning: a probabilistic perspective. See p. 434 for more details.
#'
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param regIndicators Labels of the regularized parameters (e.g. drift_eta1_eta2)
#' @param differenceApprox Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#' @param adaptiveLassoWeights weights for the adaptive lasso. If auto, defaults to the unregularized parameter estimates.
#' @author Jannik Orzek
#' @import OpenMx
#' @export
getMaxRegValue <- function(mxObject, regIndicators, differenceApprox, adaptiveLassoWeights){
  # This function is adapted from Murphy (2012) Machine learning: a probabilistic perspective. See p. 434 for more details.
  cat("Computing upper bound for regValues ... ")

  converged <- FALSE
  # warning("automatically determining the maximal lambda with getLambdaMax is still experimental! It only produces an approximation of the required lambdaMax which can be too large or too small.")

  it <- 0

  while(!converged){
    if(it == length(regIndicators)){
      stop("Error while automatically setting the regValues: The models did not converge. Try setting the regValues manually.")
    }
    # extract parameter vector
    param <- OpenMx::omxGetParameters(mxObject)

    numberRegularized <- length(regIndicators) - it

    # Problem: setting all regularized parameters to zero might result in an impossible model
    # if this error occurs, regCtsem will iteratively try to set a subset of the parameters to zero
    # the size of the parameters determines the order in which they are set to zero
    # however, this is only a rough approximation and might result in an unsatisfactory maxRegValue
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
        warning("Error when determining the regValues automatically: Setting all regularized parameters to zero resulted in an impossible model. regCtsem will try to at least set a subset of the regularized parameters to zero; however, this might result in a wrong upper bound for the regValues! Consider setting the regValues manually.")
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
        warning("Error when determining the regValues automatically: Setting all regularized parameters to zero resulted in an impossible model. regCtsem will try to at least set a subset of the regularized parameters to zero; however, this might result in a wrong upper bound for the regValues! Consider setting the regValues manually.")
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
    warning(paste0("regCtsem did set ", numberRegularized, " of the ", length(regIndicators), " regularized parameters to zero when determining the maximal regValue."))
  }
  # define maxRegValue as the maximal gradient of the regularized parameters
  maxRegValue <- max(abs(grad[regIndicators]) * adaptiveLassoWeights[regIndicators]^(-1))

  cat("DONE \n")

  return(list("maxRegValue" = maxRegValue, "sparseParameters" = param))
}




