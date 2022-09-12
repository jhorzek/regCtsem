## Functions to compute the fit for an individual person

#' individualM2LL
#'
#' compute the -2 log Likelihood for a singe person
#' @param regCtsemObject regularized model
#' @param parameterValues values for the parameters
#' @param i indicator for person(s)
#' @export
individualM2LL <- function(regCtsemObject, parameterValues, i){

  if(is.null(regCtsemObject$setup$cpptsemObject)) regCtsemObject <- restore(regCtsemObject)
  regCtsemObject$setup$cpptsemObject$setParameterValues(parameterValues, names(parameterValues))

  indM2LLs <- rep(NA,length(i))

  if(regCtsemObject$setup$objective == "ML"){
    regCtsemObject$setup$cpptsemObject$computeRAM()

    for(j in 1:length(i)){
      indM2LLs[j] <- individualM2LLML(regCtsemObject$setup$cpptsemObject,
                                      wideData = regCtsemObject$setup$dataset,
                                      i = i[j])
    }

  }else if(regCtsemObject$setup$objective == "Kalman"){
    for(j in 1:length(i)){
      indM2LLs[j] <- individualM2LLKalman(regCtsemObject$setup$cpptsemObject,
                                          i = i[j])
    }

  }else{
    stop("Model not fitted with ML or Kalman")
  }

  return(indM2LLs)
}

#' individualM2LLML
#'
#' Compute the -2log-Likelihood for a single individual
#' @param cpptsemObject cpptsemObject of class Rcpp_cpptsemRAMmodel
#' @param wideData dataset in wide format
#' @param i index for person
#' @export
individualM2LLML <- function(cpptsemObject, wideData, i){

  # extract data only
  dataset <- subset(wideData, select = !grepl("dT", colnames(wideData)) & !grepl("intervalID", colnames(wideData)) )
  expectedMeans <- cpptsemObject$expectedMeans
  expectedCovariance <- cpptsemObject$expectedCovariance

  missings <- is.na(dataset[i,])

  # check if all missing
  if(all(missings)) return(0)

  expectedMeans_i <- expectedMeans[!missings]
  expectedCovariance_i <- expectedCovariance[!missings , !missings]
  observed_i <- as.matrix(dataset[i,!missings])
  inidividualm2LL <- computeIndividualM2LL(nObservedVariables = length(observed_i),
                                                      rawData = observed_i,
                                                      expectedMeans = expectedMeans[!missings],
                                                      expectedCovariance = expectedCovariance[!missings , !missings])
  return(inidividualm2LL)
}

#' individualM2LLKalman
#'
#' Compute the -2log-Likelihood for a single individual
#' @param cpptsemObject cpptsemObject of class Rcpp_cpptsemKalmanModel
#' @param i index of individual
#' @export
individualM2LLKalman <- function(cpptsemObject, i){
  cpptsemObject$computeAndFitKalman(i)
  return(cpptsemObject$indM2LL[i,1])
}

#' individualGradient
#'
#' compute the gradients for a singe person
#' @param regCtsemObject regularized model
#' @param parameterValues values for the parameters
#' @param i indicator for person(s)
#' @param eps value for numerical approximation of gradients
#' @export
individualGradient <- function(regCtsemObject, parameterValues, i, eps = 1e-4){
  gradients <- matrix(NA, nrow = length(i), ncol = length(parameterValues))
  colnames(gradients) <- names(parameterValues)
  rownames(gradients) <- paste0("person_", i)

  for(par in 1:length(parameterValues)){
    # step forward
    parameterValues_ <- parameterValues
    parameterValues_[par] <- parameterValues_[par] + eps
    gradients[,par] <- individualM2LL(regCtsemObject = regCtsemObject, parameterValues = parameterValues_, i = i)

    # step back
    parameterValues_ <- parameterValues
    parameterValues_[par] <- parameterValues[par] - eps
    gradients[,par] <- gradients[,par] - individualM2LL(regCtsemObject = regCtsemObject, parameterValues = parameterValues_, i = i)

    # gradients
    gradients[,par] <- gradients[,par]/(2*eps)
  }

  return(gradients)

}


#' approximateLOOCV
#'
#' @param regCtsemObject regularized model
#' @param breakCondition break condition for the coordinate descent algorithm
#' @param maxIterations maximal number of iterations for the coordinate descent algorithm
#' @param gradientStepSize step size used for the gradient approximation
#' @export
approximateLOOCV <- function(regCtsemObject,
                             breakCondition = 1e-6,
                             maxIterations = 100,
                             gradientStepSize = 1e-4){
  stop("This function is under development.")
  N <- nrow(regCtsemObject$setup$dataset)
  parameterValues <- regCtsemObject$parameterEstimatesRaw
  nParameters <- nrow(parameterValues)

  LOOM2LL <- matrix(NA, nrow = N, ncol = length(regCtsemObject$setup$lambdas))
  rownames(LOOM2LL) <- paste0("person_", 1:N)
  colnames(LOOM2LL) <- paste0("lambda_", regCtsemObject$setup$lambdas)

  ipcs <- array(NA, dim = c(nrow(regCtsemObject$setup$dataset),
                            nParameters,
                            ncol(parameterValues)),
                dimnames = list("x" = paste0("person_", 1:N),
                                "y" = rownames(parameterValues),
                                "z" = paste0("lambda_", regCtsemObject$setup$lambdas))
  )

  for(lambda in 1:ncol(parameterValues)){
    currentParameterValues <- parameterValues[,lambda]
    # compute person-specific gradients
    gradients <- individualGradient(regCtsemObject = regCtsemObject, parameterValues = currentParameterValues, i = 1:N, eps = gradientStepSize)

      # we will use the hessian provided by the optimizer
      hessian <- regCtsemObject$misc$Hessians[,,lambda]

      approximateHessian <- FALSE
      NHessian <- N # sample size for the Hessian

    # Compute IPCs

    delta <- matrix(0, nrow = nParameters, ncol = 1)
    rownames(delta) <- names(currentParameterValues)

    z <- matrix(0, nrow = nParameters, ncol = 1)
    rownames(z) <- names(currentParameterValues)

    e <- matrix(0, nrow = nParameters, ncol = 1)

    for(i in 1:N){
      it <- 1
      delta[] <- 0
      delta_old <- delta

      while(TRUE){
        z[] <- 0
        for(par in sample(1:nParameters, nParameters)){
          e[] <- 0
          e[par] <- 1

            hessianXe <- hessian%*%e
            hessianXDelta <- hessian%*%delta

          if(names(currentParameterValues[par]) %in% regIndicators){
            testValue <- (t(e)%*%hessianXe)*(currentParameterValues[par] + delta[par])

            if(testValue <= NHessian*t(e)%*%matrix(gradients[i,], ncol = 1) + t(e)%*%hessianXDelta + NHessian*lambda){

              z[par] <- -(NHessian*t(e)%*%matrix(gradients[i,], ncol = 1) + t(e)%*%hessianXDelta + NHessian*lambda)/(t(e)%*%hessianXe)
              delta[par] <- delta[par] + z[par]

            }else if(testValue >= NHessian*t(e)%*%matrix(gradients[i,], ncol = 1) + t(e)%*%hessianXDelta - NHessian*lambda){

              z[par] <- -(NHessian*t(e)%*%matrix(gradients[i,], ncol = 1) + t(e)%*%hessianXDelta - NHessian*lambda)/(t(e)%*%hessianXe)
              delta[par] <- delta[par] + z[par]

            }else{

              z[par] <- -(currentParameterValues[par] + delta[par])
              delta[par] <- delta[par] + z[par]

            }

          }else{

            # unregularized parameters
            z[par] <- -(NHessian*t(e)%*%matrix(gradients[i,], ncol = 1) + t(e)%*%hessianXDelta)/(t(e)%*%hessianXe)
            delta[par] <- delta[par] + z[par]

          }
        }# end for parameter

        # check inner stopping criterion:
        breakInner <- max(abs(delta_old-delta)/(max(abs(currentParameterValues)))) < .0001
        delta_old <- delta

        if(breakInner){

          break

        }

        if(it >= maxIterations){

          break

        }
        it <- it + 1

      }# end while

      # set individual parameters
      ipcs[i,,lambda] <- currentParameterValues + delta

    }# end for i in 1:N

    for(i in 1:N){
      parametersWithoutI <- apply(ipcs[-i,,lambda], 2, mean)

      loocv_i <- individualM2LL(regCtsemObject = regCtsemObject, parameterValues = parametersWithoutI, i = i)

      if(any(class(loocv_i) == "try-error")) next
      LOOM2LL[i,lambda] <- loocv_i

    }

  }

  loocvMeans <- apply(LOOM2LL,2,mean)

  finalParameters <- parameterValues[,which(loocvMeans == min(loocvMeans))]

  return(list("IPCs" = ipcs,
              "LOOM2LL" = LOOM2LL,
              "finalParameters" = finalParameters)
  )
}
