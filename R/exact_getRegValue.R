#' exact_getRegValue
#'
#' computes sum(lambda*abs(regularized Values))
#'
#' @param regIndicators Names of regularized parameters
#' @param lambda Penaltiy value
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @export
exact_getRegValue <- function(lambda, theta, regIndicators, adaptiveLassoWeights = NULL){

  regVal <- 0
  if(!is.vector(theta)){
    thetaNames <- rownames(theta)
    # iterate over theta
    for(th in thetaNames){
      # if theta is regularized
      if(th %in% regIndicators){
        regVal <- regVal + lambda*abs(theta[th,])*ifelse(is.null(adaptiveLassoWeights),1,adaptiveLassoWeights[th])
      }
    }
    names(regVal) <- NULL
    return(regVal)
  }

  thetaNames <- names(theta)
  # iterate over theta
  for(th in thetaNames){
    # if theta is regularized
    if(th %in% regIndicators){
      regVal <- regVal + lambda*abs(theta[th])*ifelse(is.null(adaptiveLassoWeights),1,adaptiveLassoWeights[th])
    }
  }
  names(regVal) <- NULL
  return(regVal)
}


