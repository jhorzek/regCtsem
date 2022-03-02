#' ipcGlmnet
#'
#' @param parameterValues labeled vector with parameter values
#' @param fitfunction function which takes the data of a single person (vector), the parameter values and the additional arguments specified with ... as input and returns the fit value for this person
#' @param data data set with persons in rows and variables in columns
#' @param regIndicators vector with labels of regularized parameters
#' @param lambda size of the penalty value
#' @param hessian matrix with hessian
#' @param breakCondition break condition for the coordinate descent algorithm
#' @param maxIterations maximal number of iterations for the coordinate descent algorithm
#' @param gradientStepSize step size used for the gradient approximation
#' @param ... additional arguments passed to fitfunction
#' @export
ipcGlmnet <- function(parameterValues,
                      fitfunction,
                      data,
                      regIndicators,
                      lambda,
                      hessian = NULL,
                      breakCondition = 1e-6,
                      maxIterations = 100,
                      gradientStepSize = 1e-4,
                      ...){

  N <- nrow(data)
  nParameters <- length(parameterValues)

  if(is.null(names(parameterValues))) stop("parameterValues must have names")

  if(!all(regIndicators %in% names(parameterValues))) stop("Some or all regIndicators are not in the names of parameterValues.")

  # compute hessian if necessary
  if(is.null(hessian)){

    hessianFunction <- function(parameterValues, fitfunction, data, ...){
      m2LL <- sum(apply(data, 1, fitfunction, parameterValues, ...))
    }

    hessian <- numDeriv::hessian(func = hessianFunction, x = parameterValues, ...)

    if(any(eigen(hessian)$values) < 0) stop("Hessian not positive definite.")

  }

  # compute person-specific gradients
  gradients <- matrix(NA, nrow = N, ncol = nParameters)
  colnames(gradients) <- names(parameterValues)

  for(par in 1:nParameters){
    # step back
    pars <- parameterValues
    pars[par] <- pars[par] - gradientStepSize
    stepBack <- apply(data, 1, fitfunction, pars, ...)

    # step forward
    pars <- parameterValues
    pars[par] <- pars[par] + gradientStepSize
    stepForward <- apply(data, 1, fitfunction, pars, ...)

    # compute gradients
    gradients[,par] <- (stepForward - stepBack)/(2*gradientStepSize)
  }

  # Compute IPCs
  IPCs <- matrix(NA, nrow = N, ncol = nParameters)
  colnames(IPCs) <- names(parameterValues)

  delta <- matrix(0, nrow = nParameters, ncol = 1)
  rownames(delta) <- names(parameterValues)

  z <- matrix(0, nrow = nParameters, ncol = 1)
  rownames(z) <- names(parameterValues)

  e <- matrix(0, nrow = nParameters, ncol = 1)

  for(i in 1:N){
    it <- 1
    delta[] <- 0
    while(TRUE){
      z[] <- 0
      for(par in sample(1:nParameters, nParameters)){
        e[] <- 0
        e[par] <- 1

        if(names(parameterValues[par]) %in% regIndicators){
          testValue <- (t(e)%*%hessian%*%e)*(parameterValues[par] + delta[par])

          if(testValue <= N*t(e)%*%matrix(gradients[i,], ncol = 1) + t(e)%*%hessian%*%delta + N*lambda){

            z[par] <- -(N*t(e)%*%matrix(gradients[i,], ncol = 1) + t(e)%*%hessian%*%delta + N*lambda)/(t(e)%*%hessian%*%e)
            delta[par] <- delta[par] + z[par]

          }else if(testValue >= N*t(e)%*%matrix(gradients[i,], ncol = 1) + t(e)%*%hessian%*%delta - N*lambda){

            z[par] <- -(N*t(e)%*%matrix(gradients[i,], ncol = 1) + t(e)%*%hessian%*%delta - N*lambda)/(t(e)%*%hessian%*%e)
            delta[par] <- delta[par] + z[par]

          }else{

            z[par] <- -(parameterValues[par] + delta[par])
            delta[par] <- delta[par] + z[par]

          }

        }else{

          # unregularized parameters
          z[par] <- -(N*t(e)%*%matrix(gradients[i,], ncol = 1) + t(e)%*%hessian%*%delta)/(t(e)%*%hessian%*%e)
          delta[par] <- delta[par] + z[par]

        }
      }# end for parameter

      # check inner stopping criterion:
      hessTimesD <- diag(diag((1/N)*hessian))%*%z^2

      breakInnner <- max(hessTimesD) < gradientStepSize

      if(is.na(breakInnner)){

        stop("The direction z appears to go to infinity.")

      }
      if(breakInnner){

        break

      }

      if(it >= maxIterations){

        break

      }
      it <- it + 1

    }# end while

    # set individual parameters
    IPCs[i,] <- parameterValues + delta

  }# end for i in 1:N

  loocv <- rep(NA, N)
  for(i in 1:N){
    parametersWithoutI <- apply(IPCs[-i,], 2, mean)

    loocv_i <- try(fitfunction(data[i,], parametersWithoutI, ...), silent = TRUE)
    if(any(class(loocv_i) == "try-error")) next
    loocv[i] <- loocv_i

  }

  return(list("IPCs" = IPCs,
              "loocv" = loocv
  )
  )
}

#' returnIndividualFit
#'
#' This function returns the fit of the model for the vector individualData. It is currently not very efficient
#'
#' @param individualData vector with data of one individual
#' @param parameterValues parameter values
#' @param fullData matrix with entire data set
#' @param cpptsemObject object created with cpptsemFromCtsem
#' @export
returnIndividualFit <- function(individualData, parameterValues, fullData, cpptsemObject){

  # determine if Kalman or ML
  if(any(class(cpptsemObject) == "Rcpp_cpptsemKalmanModel")){
    objective <- "Kalman"
  }else{
    objective <- "ML"
  }

  if(objective == "ML"){
    cpptsemObject$setParameterValues(parameterValues, names(parameterValues))
    cpptsemObject$computeRAM()

    # in case of ML, we have to calculate the individual fit by hand
    m2LL_i <- regCtsem:::computeIndividualM2LL(nObservedVariables = sum(!is.na(individualData)),
                                               rawData = individualData[!is.na(individualData)],
                                               expectedMeans = cpptsemObject$expectedMeans[!is.na(individualData)],
                                               expectedCovariance = cpptsemObject$expectedCovariance[!is.na(individualData),
                                                                                                     !is.na(individualData)]
    )
    return(m2LL_i)
  }

  # in case of Kalman filter, the individual likelihoods are already computed; however we have to fit the full model, which
  # is quite inefficient; this has to be changed later on...
  if(objective == "Kalman"){
    i <- which(apply(fullData, 1, function(x) identical(x, individualData)))[1]

    return(cpptsemObject$m2LLs[i])
  }


}
