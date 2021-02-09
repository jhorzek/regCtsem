#' exact_getSubgradients
#'
#' Computes the subgradients of a lasso penalized likelihood f(theta) = L(theta)+lambda*p(theta)
#'
#' @param theta vector with named parameters
#' @param jacobian derivative of L(theta)
#' @param regIndicators names of regularized parameters
#' @param lambda lambda value
#' @export
exact_getSubgradients <- function(theta, jacobian, regIndicators, lambda, lineSearch, adaptiveLassoWeightsMatrix){

  # first part: derivative of Likelihood
  subgradient <- jacobian

  # second part: derivative of penalty term
  ## Note: if adaptiveLassoWeightsMatrix*parameter != 0, the derivative is: lambda*adaptiveLassoWeight*sign(parameter)
  ## if adaptiveLassoWeightsMatrix*parameter == 0, the derivative is in [-lambda*adaptiveLassoWeight, +lambda*adaptiveLassoWeight]
  for(regularizedParameterLabel in regIndicators){
    absoluteValueOf <- theta[regularizedParameterLabel,]
    if(absoluteValueOf != 0){
      penaltyGradient <- adaptiveLassoWeightsMatrix[regularizedParameterLabel,regularizedParameterLabel]*lambda*sign(absoluteValueOf)
      subgradient[regularizedParameterLabel,] <- subgradient[regularizedParameterLabel,] + penaltyGradient
    }else{
      penaltyGradient <- c(-adaptiveLassoWeightsMatrix[regularizedParameterLabel,regularizedParameterLabel]*lambda, adaptiveLassoWeightsMatrix[regularizedParameterLabel,regularizedParameterLabel]*lambda)
      # check if likelihood gradient is within interval:
      setZero <- (subgradient[regularizedParameterLabel,] > penaltyGradient[1]) && (subgradient[regularizedParameterLabel,] < penaltyGradient[2])
      subgradient[regularizedParameterLabel,] <- ifelse(setZero,
                                                        0,
                                                        sign(subgradient[regularizedParameterLabel,])*(abs(subgradient[regularizedParameterLabel,]) -
                                                                                                         adaptiveLassoWeightsMatrix[regularizedParameterLabel,
                                                                                                                                    regularizedParameterLabel]*lambda))
    }
  }
  return(subgradient)
}



