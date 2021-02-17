#' exact_getBFGS
#'
#' computes the BFGS Hessian approximation
#' @param theta_k Theta at iteration k
#' @param g_k Gradients of the likelihood function at iteration k
#' @param H_k Hessian of the likelihood function at iteration k
#' @param theta_kp1 Theta at iteration k+1
#' @param g_kp1 Gradients of the likelihood function at iteration k+1
#' @param cautios boolen: should the update be skipped if it would result in a non positive definite Hessian?
#' @param hessianEps controls when the update of the Hessian approximation is skipped
#' @export
exact_getBFGS <- function(theta_k, g_k, H_k, theta_kp1, g_kp1, cautious = TRUE, hessianEps = .001){

  y <- g_kp1-g_k
  d <- theta_kp1-theta_k

  # test if positive definiteness is ensured
  skipUpdate <- try((t(y)%*%d < hessianEps) && cautious, silent = TRUE)
  if(any(class(skipUpdate) == "try-error") || skipUpdate){
    # Hessian might become non-positive definite. Return without update
    return(H_k)
  }
  if(t(y)%*%d < 0){
    warning("Hessian update possibly non-positive definite.")
  }

  H_kp1 <- H_k - (H_k%*%d%*%t(d)%*%H_k)/as.numeric(t(d)%*%H_k%*%d) + (y%*%t(y))/as.numeric(t(y)%*%d)

  if(anyNA(H_kp1)){
    warning("Invalid Hessian. Returning previous Hessian")
    return(H_k)
  }
  HessianEigen <- eigen(H_kp1)
  iscomplex <- any(!Im(HessianEigen$values) == 0)
  if(iscomplex){
    eigenVectors <- Re(HessianEigen$vectors)
    eigenValues <- Re(HessianEigen$values)
    H_kp1 <- eigenVectors%*%diag(eigenValues)%*%t(eigenVectors)
  }
  if(any(HessianEigen$values < 0)){
    while(any(eigen(H_kp1)$values < 0)){
      H_kp1 <- H_kp1 + .01*diag(nrow(H_kp1))
    }
  }

  return(H_kp1)

}



