#' exact_armijoLineSearch
#'
#' performs Armijo line search
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
#' @param stepSize Initial stepsize of the outer iteration (theta_{k+1} = theta_k + Stepsize \* Stepdirection)
#' @param differenceApprox which approximation for the gradients should be used? Recommended is central
#' @param maxIter_line maximal number of iterations for line search
#'
#' @export
exact_armijoLineSearch <- function(gradientModel, adaptiveLassoWeights, thetaNames, regIndicators, lambda,
                             theta_kp1, m2LL_kp1, g_kp1, d, c1 = 0.0001, stepSize = .99, differenceApprox = "central", maxIter_line = 100){
  stop("lineSearch = 'armijo' is deprecated. Use lineSearch = 'wolfe'.")
  if(!is.null(adaptiveLassoWeights)){stop("not implemented for adaptive lasso")}
  # Adapted from Huang lslx

  # get penalized M2LL for step size 0:
  h_0 <- m2LL_kp1 + lambda*sum(abs((theta_kp1)[regIndicators,]))

  # get (sub-)gradients for step size 0:
  g_0 <- regCtsem::exact_getSubgradients(theta = theta_kp1, jacobian = g_kp1, regIndicators = regIndicators, lambda = lambda, lineSearch = "armijo")

  # Inexact Line Search
  i <- 0

  while(i<maxIter_line){
    i <- i+1 # update iterations
    if(stepSize <.00001){
      return(stepSize)
    }
    # new theta
    theta_kp1_td <- theta_kp1+stepSize*d
    # get L(x+stepSize*d) and L'(x+stepSize*d)
    gradientModel <- OpenMx::omxSetParameters(model = gradientModel, labels = thetaNames, values = theta_kp1_td)
    gradientModel <- OpenMx::mxRun(gradientModel, silent = TRUE)
    m2LL_kp1_td <- gradientModel$fitfunction$result[[1]]
    g_kp1_td <- gradientModel$compute$steps[[1]]$output[["gradient"]]
    g_kp1_td <- g_kp1_td[,differenceApprox] # use specified gradient approximation
    g_kp1_td <- matrix(g_kp1_td, nrow = length(g_kp1_td), ncol = 1)
    rownames(g_kp1_td) <- thetaNames

    # compute h(stepSize) = L(x+td) + p(x+td) - L(x) - p(x), where p(x) is the penalty function
    h_t <- m2LL_kp1_td + lambda*sum(abs((theta_kp1_td)[regIndicators,]))
    # compute h'(stepSize)
    g_t <- regCtsem::exact_getSubgradients(theta = theta_kp1_td, jacobian = g_kp1_td, regIndicators = regIndicators, lambda = lambda, lineSearch = "armijo")

    # Check Armijo
    if(h_t-h_0 <= c1*stepSize*(t(g_kp1)%*%d+lambda*sum(abs((theta_kp1_td)[regIndicators,]))-lambda*sum(abs((theta_kp1)[regIndicators,])))){
      return(stepSize)
    }else{
      stepSize <- stepSize^2
    }
  }
}



