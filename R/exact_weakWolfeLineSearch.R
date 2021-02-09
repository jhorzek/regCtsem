#' exact_weakWolfeLineSearch
#'
#' @param gradientModel mxObject for computing the derivarive of the likelihood with respect to the parameter estimates
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param thetaNames names of the parameter estimates
#' @param regIndicators vector with names of parameters to regularize
#' @param lambda penalty value
#' @param theta_kp1 parameter values of iteration k plus 1
#' @param m2LL_kp1 -2 log likelihood of iteration k plus 1
#' @param g_kp1 gradients of iteration k plus 1
#' @param d vector with updates to parameter estimates
#' @param c1 tuning parameter for Armijo condition
#' @param c2 c2 constant for lineSearch. This constant controls the Curvature condition
#' @param stepSize Initial stepsize of the outer iteration (theta_{k+1} = theta_k + Stepsize \* Stepdirection)
#' @param differenceApprox which approximation for the gradients should be used? Recommended is central
#' @param eps_numericDerivative controls the precision of the central gradient approximation. The default (1.1 * 10^(-16))^(1/3) is derived in Nocedal, J., & Wright, S. J. (2006). Numerical optimization (2nd ed), p. 197
#' @param maxIter_line maximal number of iterations for line search
#' @param verbose 0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress
#' @param epsWW additional parameter to shorten the search if the upper and lower bound are already very close to each other
#' @export
exact_weakWolfeLineSearch <- function(gradientModel, objective, gradientModelcpp = NULL, adaptiveLassoWeights, thetaNames, regIndicators, lambda,
                                      theta_kp1, m2LL_kp1, g_kp1, d, c1 = 0.0001, c2 = .9,
                                      stepSize = 1, differenceApprox = "central",
                                      maxIter_line = 100, verbose = 0, eps_WW = .0001, eps_numericDerivative = (1.1 * 10^(-16))^(1/3)){
  # Adapted from HANSO (see LEWIS, D.S. &  OVERTONNON, M. L., (2013). SMOOTH OPTIMIZATION VIA BFGS)

  if(!is.null(gradientModelcpp)){
    # deep clone of gradientModelcpp
    gradientModelcpp_i <- rlang::duplicate(gradientModelcpp, shallow =FALSE)
  }

  if(!is.null(adaptiveLassoWeights)){
    adaptiveLassoWeightsMatrix <- diag(adaptiveLassoWeights)
  }else{
    adaptiveLassoWeightsMatrix <- diag(length(thetaNames))
  }
  rownames(adaptiveLassoWeightsMatrix) <- thetaNames
  colnames(adaptiveLassoWeightsMatrix) <- thetaNames

  stepSize_L <- 0 # lower bound
  stepSize_U <- Inf # upper bound

  # get penalized M2LL for step size 0:

  h_0 <- m2LL_kp1 + lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1)[regIndicators,]))

  # get (sub-)gradients for step size 0:
  g_0 <- regCtsem::exact_getSubgradients(theta = theta_kp1,
                                          jacobian = g_kp1,
                                          regIndicators = regIndicators,
                                          lambda = lambda,
                                          lineSearch = "wolfe",
                                          adaptiveLassoWeightsMatrix = adaptiveLassoWeightsMatrix)

  # Inexact Line Search
  i <- 0

  while(i<maxIter_line){
    i <- i+1 # update iterations

    if(verbose==2){
      cat(paste0("\r",
                 "## Linesearch: [Lower step size: ", sprintf('%.3f', stepSize_L, 3),
                 ", current step size: ", sprintf('%.3f', stepSize, 3),
                 ", upper step size: ", sprintf('%.3f', stepSize_U, 3),
                 "] ##"
      )
      )
      flush.console()
    }

    if(stepSize_U-stepSize_L < eps_WW){
      # The following part is adapted from https://github.com/ceres-solver/ceres-solver/blob/master/internal/ceres/line_search.cc
      # and is a workaround for when the line search algorithm does not find
      # an exact point that satisfies both Wolfe conditions and results in an upper and lower bound which are
      # approximately identical. In this case, the iterations are stopped and the current lower bound is returned
      if(verbose>1){
        warning("Could not find an exact point satisfying both Wolfe conditions.")
      }
      # an additional problem when the hessian is not well approximated is that the search direction might be wrong.
      # this would result in an increasing Likelihood which is indicated by a lower bound of 0.
      if(stepSize_L == 0){warning("Searchdirection might be problematic. Line search returned step size 0.")}
      return(stepSize_L)
    }

    # new theta
    theta_kp1_td <- rlang::duplicate(theta_kp1)+stepSize*d
    # get L(x+stepSize*d) and L'(x+stepSize*d)
    if(!is.null(gradientModelcpp)){
      gradientModelcpp_i$setParameterValues(theta_kp1_td, thetaNames)
      if(tolower(objective) == "ml"){
      invisible(capture.output(out1 <- try(gradientModelcpp_i$computeRAM(), silent = TRUE), type = "message"))
      invisible(capture.output(out2 <- try(gradientModelcpp_i$fitRAM(), silent = TRUE), type = "message"))
      }else{
        invisible(capture.output(out1 <- try(gradientModelcpp_i$computeAndFitKalman(), silent = TRUE), type = "message"))
        out2 <- NA
      }
      if(any(class(out1)== "try-error") | any(class(out2)== "try-error") | !is.finite(gradientModelcpp_i$m2LL)){
        # if the starting values are not feasible
        stepSize <- .5*stepSize
        next
        }

      if(any(class(gradientModelcpp_i)=="try-error")){
        # if the starting values are far off, the gradients can be very large and testing for the initial
        # step size might result in infeasible parameter values. In this case: try smaller step size
        stepSize <- .5*stepSize
        next
      }
    }else{
      gradientModel_i <- OpenMx::omxSetParameters(model = gradientModel, labels = thetaNames, values = theta_kp1_td)
      gradientModel_i <- suppressWarnings(try(OpenMx::mxRun(gradientModel_i, silent = TRUE), silent = TRUE))
      if(any(class(gradientModel_i)=="try-error")){
        # if the starting values are far off, the gradients can be very large and testing for the initial
        # step size might result in infeasible parameter values. In this case: try smaller step size
        stepSize <- .5*stepSize
        next
      }
      if(!is.na(gradientModel_i$output$status$code)){
        if(gradientModel_i$output$status$code == 10){
          # if the starting values are not feasible
          stepSize <- .5*stepSize
          next}
      }

    }


    if(!is.null(gradientModelcpp)){
      m2LL_kp1_td <- gradientModelcpp_i$m2LL
      if(tolower(objective) == "ml"){
      invisible(capture.output(g_kp1_td <- try(gradientModelcpp_i$approxRAMGradients(eps_numericDerivative)[thetaNames], silent = TRUE), type = "message"))
      }else{
        invisible(capture.output(g_kp1_td <- try(gradientModelcpp_i$approxKalmanGradients(eps_numericDerivative)[thetaNames], silent = TRUE), type = "message"))
        }
      if(any(class(g_kp1_td)=="try-error")){
        # sometimes the stepsize of 1 results in NA; in this case: use smaller step size. Also, remove unfeasible starting values
        stepSize <- .9*stepSize
        next
      }
      g_kp1_td <- matrix(g_kp1_td, nrow = length(g_kp1_td), ncol = 1)
      rownames(g_kp1_td) <- thetaNames

      if(any(is.na(g_kp1_td))){
        # sometimes the stepsize of 1 results in NA; in this case: use smaller step size. Also, remove unfeasible starting values
        stepSize <- .9*stepSize
        next
      }
    }else{
      m2LL_kp1_td <- gradientModel_i$fitfunction$result[[1]]
      g_kp1_td <- gradientModel_i$compute$steps[[1]]$output[["gradient"]]
      g_kp1_td <- g_kp1_td[,differenceApprox] # use specified gradient approximation
      g_kp1_td <- matrix(g_kp1_td, nrow = length(g_kp1_td), ncol = 1)
      rownames(g_kp1_td) <- thetaNames

      if((!is.na(gradientModel_i$output$status$code) & gradientModel_i$output$status$code == 10) | any(is.na(g_kp1_td))){
        # sometimes the stepsize of 1 results in NA; in this case: use smaller step size. Also, remove unfeasible starting values
        stepSize <- .9*stepSize
        next
      }
    }




    # compute h(stepSize) = L(x+td) + p(x+td) - L(x) - p(x), where p(x) is the penalty function
    h_t <- m2LL_kp1_td +  lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1_td)[regIndicators,]))
    # compute h'(stepSize)
    g_t <- regCtsem::exact_getSubgradients(theta = theta_kp1_td,
                                            jacobian = g_kp1_td,
                                            regIndicators = regIndicators,
                                            lambda = lambda,
                                            lineSearch = "wolfe",
                                            adaptiveLassoWeightsMatrix = adaptiveLassoWeightsMatrix)

    # Check Armijo and Curvature Condition
    Armijo <- h_t <= h_0 + c1*stepSize*t(d)%*%g_0 # for first wolfe condition (Armijo condition)
    Curvature <- t(g_t)%*%d >= c2*t(d)%*%g_0 # for second wolfe condition (curvature condition)

    if (!Armijo){
      # if the first condition is not satisfied at stepSize (i.e., stepSize is not resulting in a sufficient decrease of the function), set the upper limit to stepSize
      stepSize_U <- stepSize
    }else if(!Curvature){
      # if stepSize results in a sufficient decrease, but the second wolfe condition is not satisfied at stepSize (i.e. the step size is too small), set the lower limit to stepSize
      stepSize_L <- stepSize
    }else{
      # if both conditions are satisfied: end line search
      return(stepSize)
    }
    if(stepSize_U < Inf){
      # if the upper limit stepSize_U has been updated, set stepSize to the mean of the interval (lower, upper)
      stepSize <- (stepSize_L+stepSize_U)/2
    }else{
      # if the upper limit has not been updated, set t to 2 times the lower limit
      stepSize <- 2*stepSize_L
    }
  }
  return(stepSize)
}



