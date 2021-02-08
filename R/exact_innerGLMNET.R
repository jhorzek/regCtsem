#' exact_innerGLMNET
#'
#' performs the inner optimization routine of GLMNET. exact_innerGLMNET returns a direction vector for the next step in the outer optimization routine.
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param thetaNames Vector with names of theta-parameters
#' @param regIndicators Names of regularized parameters
#' @param lambda Penaltiy value
#' @param theta_kp1 Theta at iteration k+1
#' @param g_kp1 Gradients of the likelihood function at iteration k+1
#' @param H_kp1 Hessian of the likelihood function at iteration k+1
#' @param maxIter_in Maximal number of iterations of the inner optimization algorithm
#' @param eps_in Stopping criterion for the inner iterations
#' @export
exact_innerGLMNET <- function(adaptiveLassoWeights, thetaNames, regIndicators, lambda, theta_kp1, g_kp1, H_kp1, maxIter_in, eps_in){
  ## inner loop: optimize directions

  d <- matrix(0,nrow = length(thetaNames), ncol = 1) # initialize step vector (direction d)
  rownames(d) = thetaNames

  iter_in <- 0

  while(iter_in < maxIter_in){
    iter_in <- iter_in+1

    # initialize direction z
    z <- matrix(0,nrow = length(thetaNames), ncol = 1)

    # random order of updates
    d_order <- sample(1:length(d),length(d))

    # iterate over parameters
    for(d_i in d_order){

      # compute derivative elements:
      dp_k <- g_kp1[d_i]+(H_kp1%*%d)[d_i]
      d2p_k <- H_kp1[d_i,d_i]

      # if the parameter is regularized:
      if(names(d[d_i,]) %in% regIndicators){
        theta_i <- theta_kp1[d_i]
        adaptiveWeight <- ifelse(is.null(adaptiveLassoWeights), 1, adaptiveLassoWeights[d_i])

        # adjust d for regularized parameters
        if((dp_k-lambda*adaptiveWeight)>=(d2p_k*(theta_i+d[d_i]))){
          # condition 1
          z_j <- -(dp_k-lambda*adaptiveWeight)/(d2p_k)
          z[d_i] <- z_j
          d[d_i] <- d[d_i] + z_j
        }else if((dp_k+lambda*adaptiveWeight)<=(d2p_k*(theta_i+d[d_i]))){
          # condition 2
          z_j <- -(dp_k+lambda*adaptiveWeight)/(d2p_k)
          z[d_i] <- z_j
          d[d_i] <- d[d_i] + z_j
        }else{
          # condition 3
          z_j <- -(theta_i+d[d_i])
          z[d_i] <- z_j
          d[d_i] <- d[d_i]+z_j
        }
      }else{
        # if not regularized: coordinate descent with newton direction
        z_j <- -dp_k/d2p_k
        z[d_i] <- z_j
        d[d_i] <- d[d_i]+z_j
      }
    }

    # check inner stopping criterion:
    HessTimesD <- diag(diag(H_kp1))%*%z^2

    breakInnner <- max(HessTimesD)<eps_in
    if(is.na(breakInnner)){
      stop("Error in inner exact_innerGLMNET: The direction z appears to go to infinity.")
    }
    if(breakInnner){
      break
    }
  }

  # return step direction
  return(d)
}


