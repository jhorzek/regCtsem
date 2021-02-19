#' exact_GLMNETLineSearch
#'
#' performs the line search procedure described by Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421 Equation 20.
#'
#' @param gradientModel Object of Type MxModel which specifies how the gradients of the likelihood-function are computed (the jacobian)
#' @param objective which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)
#' @param gradientModelcpp cpptsem object which specifies how the gradients of the likelihood-function are computed (the jacobian)
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param thetaNames names of the parameter estimates
#' @param regIndicators vector with names of parameters to regularize
#' @param lambda penalty value
#' @param theta_kp1 parameter values of iteration k plus 1
#' @param m2LL_kp1 -2 log likelihood of iteration k plus 1
#' @param g_kp1 gradients of iteration k plus 1
#' @param d vector with updates to parameter estimates
#' @param differenceApprox which approximation for the gradients should be used? Recommended is central
#' @param eps_numericDerivative controls the precision of the central gradient approximation. The default (1.1 * 10^(-16))^(1/3) is derived in Nocedal, J., & Wright, S. J. (2006). Numerical optimization (2nd ed), p. 197
#' @param stepSize Initial stepsize of the outer iteration (theta_{k+1} = theta_k + Stepsize \* Stepdirection)
#' @param sig GLMNET & GIST: GLMNET: only relevant when lineSearch = 'GLMNET' | GIST: sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
#' @param gam GLMNET when lineSearch = 'GLMNET'. Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.
#' @param maxIter_line maximal number of iterations for line search
#' @export
exact_GLMNETLineSearch <- function(gradientModel, objective, gradientModelcpp,
                                   adaptiveLassoWeights, thetaNames, regIndicators,
                                   lambda, theta_kp1, m2LL_kp1, g_kp1, H_k,
                                   d, differenceApprox, eps_numericDerivative,
                                   stepSize, sig, gam = 0, maxIter_line){

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

  # get penalized M2LL for step size 0:

  f_0 <- m2LL_kp1 + lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1)[regIndicators,]))
  pen_0 <- lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1)[regIndicators,]))

  i <- 0

  if(stepSize >= 1){
    stepSizeInit <- .9
  }else{
    stepSizeInit <- stepSize
  }
  while(TRUE){
    stepSize <- stepSizeInit^i
    theta_kp1_td <- theta_kp1+stepSize*d

    # compute new fitfunction value
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
        i <- i+1
        next
      }

      if(any(class(gradientModelcpp_i)=="try-error")){
        # if the starting values are far off, the gradients can be very large and testing for the initial
        # step size might result in infeasible parameter values. In this case: try smaller step size
        i <- i+1
        next
      }
    }else{
      gradientModel_i <- OpenMx::omxSetParameters(model = gradientModel, labels = thetaNames, values = theta_kp1_td)
      gradientModel_i <- suppressWarnings(try(OpenMx::mxRun(gradientModel_i, silent = TRUE), silent = TRUE))
      if(any(class(gradientModel_i)=="try-error")){
        # if the starting values are far off, the gradients can be very large and testing for the initial
        # step size might result in infeasible parameter values. In this case: try smaller step size
        i <- i+1
        next
      }
      if(!is.na(gradientModel_i$output$status$code)){
        if(gradientModel_i$output$status$code == 10){
          # if the starting values are not feasible
          i <- i+1
          next}
      }

    }

    if(!is.null(gradientModelcpp)){
      m2LL_kp1_td <- gradientModelcpp_i$m2LL
    }else{
      m2LL_kp1_td <- gradientModel_i$fitfunction$result[[1]]

      if((!is.na(gradientModel_i$output$status$code) & gradientModel_i$output$status$code == 10)){
        # sometimes the stepSize of 1 results in NA; in this case: use smaller step size. Also, remove unfeasible starting values
        i <- i+1
        next
      }
    }

    # compute h(stepSize) = L(x+td) + p(x+td) - L(x) - p(x), where p(x) is the penalty function
    f_new <- m2LL_kp1_td +  lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1_td)[regIndicators,]))
    p_new <- lambda*sum(abs((adaptiveLassoWeightsMatrix%*%theta_kp1_td)[regIndicators,]))

    # test line search criterion
    lineCriterion <- f_new - f_0 <= sig*stepSize*(t(g_kp1)%*%d + gam*t(d)%*%H_k%*%d + p_new - pen_0)
    if(lineCriterion){
      break
    }
    i <- i+1
    if(i >= maxIter_line){
      warning("Line search found no stepSize within the maximal number of line search iterations.")
      break
    }
  }
  return(stepSize)
}





