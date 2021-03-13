#' controlApprox
#'
#' The following arguments can be used to adjust the approximate optimization
#'
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
#' @param zeroThresh threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
#' @export
controlApprox <- function(){
  return(controlApprox <- list(
    "epsilon" = .001, # epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
    "zeroThresh" = .04 # threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
  )
  )
}

#' controlGIST
#'
#' The following arguments can be used to adjust the GIST optimization
#'
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
#' @param tryCpptsem should regCtsem try to translate the model to cpptsem? This can speed up the computation considerably but might fail for some models
#' @param forceCpptsem should cpptsem be enforced even if results differ from ctsem? Sometimes differences between cpptsem and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
#' @param stepSize initial step size of the outer iteration
#' @param sig sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
#' @param differenceApprox Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#' @param maxIter_out Maximal number of outer iterations
#' @param maxIter_in Maximal number of inner iterations
#' @param break_outer Stopping criterion for outer iterations. It has to be a named value. By default (name: gradient), a relative first-order condition is checked, where the maximum absolute value of the gradients is compared to break_outer (see https://de.mathworks.com/help/optim/ug/first-order-optimality-measure.html). Alternatively, an absolute tolerance can be passed to the function (e.g., break_outer = c("gradient" = .0001)). Instead of relative gradients, the change in parameters can used as breaking criterion. To this end, use c("parameterChange" = .00001)
#' @param eta if the current step size fails, eta will decrease the step size. Must be > 1
#' @param stepsizeMin Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param stepsizeMax Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param GISTLinesearchCriterion criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
#' @param GISTNonMonotoneNBack in case of non-monotone line search: Number of preceding regM2LL values to consider
#' @param approxFirst Should approximate optimization be used first to obtain start values for exact optimization? 1 and 2 are using OpenMx with 1 = optimization only for first lambda, 2 = optimization for all lambdas. 3 ensures that the fit will not be worse than in the sparse model if lambdas = "auto" or sparseParameters are provided. To this end, 10 models between the current parameter estimates and the sparse parameter estimates are tested and the one with the lowest regM2LL is used for starting values. 4 = optimizing using optim or OpenMx if cpptsem is not available, 5 = optimizing using Rsolnp or OpenMx if cpptsem is not available (requires installation of Rsolnp). "auto" will default to 3 if lambdas = "auto" and 4 otherwise
#' @param numStart Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param approxOpt Used if approxFirst = 3. Should each of the generated starting values be optimized slightly? This can substantially improve the fit of the generated starting values. 1 = optimization with optim, 2 = optimization with Rsolnp
#' @param approxMaxIt Used if approxFirst = 3 and approxOpt > 1. How many outer iterations should be given to each starting values vector? More will improve the selected starting values but slow down the computation. If approxFirst =  4, or approxFirst = 5 this will control the number of outer iteration in optim or solnp .
#' @export
controlGIST <- function(){
  return(controlGIST <- list(
    "tryCpptsem" = TRUE, # should regCtsem try to translate the model to C++? This can speed up the computation considerably but might fail for some models
    "forceCpptsem" = FALSE, # should the C++ translation be enforced even if results differ from ctsem? Sometimes differences between the C++ implementation and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
    "stepSize" =  1, # initial step size of the outer iteration
    "sig" = 10^(-5), # sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
    "maxIter_out" = 100, # Maximal number of outer iterations
    "maxIter_in" = 1000, # Maximal number of inner iterations
    "break_outer" = c("gradient" = "max(max(abs(startingValues))*.001, .001)"), # Stopping criterion for outer iterations. It has to be a named value. By default (name: gradient), a relative first-order condition is checked, where the maximum absolute value of the gradients is compared to break_outer (see https://de.mathworks.com/help/optim/ug/first-order-optimality-measure.html). Alternatively, an absolute tolerance can be passed to the function (e.g., break_outer = c("gradient" = .0001)). Instead of relative gradients, the change in parameters can used as breaking criterion. To this end, use c("parameterChange" = .00001)
    "eta" = 2, # if the current step size fails, eta will decrease the step size. Must be > 1
    "stepsizeMin" = 1/(10^30), # Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
    "stepsizeMax" = 10^30, # Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
    "GISTLinesearchCriterion" = "monotone", # criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
    "GISTNonMonotoneNBack" = 5,# in case of non-monotone line search: Number of preceding regM2LL values to consider
    "approxFirst" = "auto", # Should approximate optimization be used first to obtain start values for exact optimization? 1 and 2 are using OpenMx with 1 = optimization only for first lambda, 2 = optimization for all lambdas. 3 ensures that the fit will not be worse than in the sparse model if lambdas = "auto" or sparseParameters are provided. To this end, numStart + 2 models between the current parameter estimates and the sparse parameter estimates are tested and the one with the lowest regM2LL is used for starting values. 4 = optimizing using optim or OpenMx if cpptsem is not available, 5 = optimizing using Rsolnp or OpenMx if cpptsem is not available (requires installation of Rsolnp)
    "numStart" = 0, # Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
    "approxOpt" = 1, # Used if approxFirst = 3. Should each of the generated starting values be optimized slightly? This can substantially improve the fit of the generated starting values. 1 = optimization with optim, 2 = optimization with Rsolnp
    "approxMaxIt" = "auto", # Used if approxFirst = 3 and approxOpt > 1. How many outer iterations should be given to each starting values vector? More will improve the selected starting values but slow down the computation. If approxFirst =  4, or approxFirst = 5 this will control the number of outer iteration in optim or solnp .
    "differenceApprox" = "central" # Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
  )
  )
}

#' controlGLMNET
#'
#' The following arguments can be used to adjust the GLMNET optimization
#'
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
#' @param tryCpptsem should regCtsem try to translate the model to cpptsem? This can speed up the computation considerably but might fail for some models
#' @param forceCpptsem should cpptsem be enforced even if results differ from ctsem? Sometimes differences between cpptsem and ctsem can result from problems with numerical precision which will lead to the m,atrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
#' @param stepSize initial step size of the outer iteration
#' @param lineSearch String indicating which linesearch should be used. Defaults to the one described in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Alternatively (not recommended) Wolfe conditions (lineSearch = "Wolfe") can be used in the outer iteration. Setting to "none" is also not recommended!.
#' @param c1 c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
#' @param c2 c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
#' @param sig for lineSearch = 'GLMNET': Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 < sigma < 1
#' @param gam for lineSearch = 'GLMNET': Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 <= gamma < 1
#' @param differenceApprox Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
#' @param initialHessianApproximation Which initial hessian approximation should be used? Possible are: 'ident' for an identity matrix and 'OpenMx' (here the hessian approxmiation from the mxObject is used). If the Hessian from 'OpenMx' is not positive definite, the negative Eigenvalues will be 'flipped' to positive Eigenvalues. This works sometimes, but not always. Alternatively, a matrix can be provided which will be used as initial Hessian
#' @param maxIter_out Maximal number of outer iterations
#' @param maxIter_in Maximal number of inner iterations
#' @param maxIter_line Maximal number of iterations for the lineSearch procedure
#' @param eps_out Stopping criterion for outer iterations
#' @param eps_in Stopping criterion for inner iterations
#' @param eps_WW Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned
#' @param approxFirst Should approximate optimization be used first to obtain start values for exact optimization? 1 and 2 are using OpenMx with 1 = optimization only for first lambda, 2 = optimization for all lambdas. 3 ensures that the fit will not be worse than in the sparse model if lambdas = "auto" or sparseParameters are provided. To this end, 10 models between the current parameter estimates and the sparse parameter estimates are tested and the one with the lowest regM2LL is used for starting values. 4 = optimizing using optim or Opemx if cpptsem is not available, 5 = optimizing using optim or Opemx if cpptsem is not available (requires installation of Rsolnp)
#' @param numStart Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param approxOpt Used if approxFirst = 3. Should each of the generated starting values be optimized slightly? This can substantially improve the fit of the generated starting values. 1 = optimization with optim, 2 = optimization with Rsolnp
#' @param approxMaxIt Used if approxFirst = 3 and approxOpt > 1. How many outer iterations should be given to each starting values vector? More will improve the selected starting values but slow down the computation. If approxFirst =  4, or approxFirst = 5 this will control the number of outer iteration in optim or solnp .
#' @export
controlGLMNET <- function(){
  return(
    controlGLMNET <- list(
      "tryCpptsem" = TRUE, # should regCtsem try to translate the model to C++? This can speed up the computation considerably but might fail for some models
      "forceCpptsem" = FALSE, # should the C++ translation be enforced even if results differ from ctsem? Sometimes differences between the C++ implementation and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
      "stepSize" =  1, # initial step size of the outer iteration
      "lineSearch" = "GLMNET", # String indicating which linesearch should be used. Defaults to the one described in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Alternatively (not recommended) Wolfe conditions (lineSearch = "Wolfe") can be used in the outer iteration. Setting to "none" is also not recommended!.
      "c1" = .0001, # c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
      "c2" = .9, # c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
      "sig" = 10^(-5), # for lineSearch = 'GLMNET': Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 < sigma < 1
      "gam" = 0, # for lineSearch = 'GLMNET': Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 <= gamma < 1
      "initialHessianApproximation" = "OpenMx", # Which initial hessian approximation should be used? Possible are: 'ident' for an identity matrix and 'OpenMx' (here the hessian approxmiation from the mxObject is used). If the Hessian from 'OpenMx' is not positive definite, the negative Eigenvalues will be 'flipped' to positive Eigenvalues. This works sometimes, but not always. Alternatively, a matrix can be provided which will be used as initial Hessian
      "maxIter_out" = 100, # Maximal number of outer iterations
      "maxIter_in" = 1000, # Maximal number of inner iterations
      "maxIter_line" = 500, # Maximal number of iterations for the lineSearch procedure
      "eps_out" = .0000000001, # Stopping criterion for outer iterations
      "eps_in" = .0000000001, # Stopping criterion for inner iterations
      "eps_WW" = .0001, #Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned
      "approxFirst" = "auto", # Should approximate optimization be used first to obtain start values for exact optimization? 1 and 2 are using OpenMx with 1 = optimization only for first lambda, 2 = optimization for all lambdas. 3 ensures that the fit will not be worse than in the sparse model if lambdas = "auto" or sparseParameters are provided. To this end, 10 models between the current parameter estimates and the sparse parameter estimates are tested and the one with the lowest regM2LL is used for starting values. 4 = optimizing using optim or OpenMx if cpptsem is not available, 5 = optimizing using Rsolnp or OpenMx if cpptsem is not available (requires installation of Rsolnp). "auto" will default to 3 if lambdas = "auto" and 4 otherwise
      "numStart" = 0, # Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
      "approxOpt" = 1, # Used if approxFirst = 3. Should each of the generated starting values be optimized slightly? This can substantially improve the fit of the generated starting values. 1 = optimization with optim, 2 = optimization with Rsolnp
      "approxMaxIt" = "auto", # Used if approxFirst = 3 and approxOpt > 1. How many outer iterations should be given to each starting values vector? More will improve the selected starting values but slow down the computation. If approxFirst =  4, or approxFirst = 5 this will control the number of outer iteration in optim or solnp .
      "differenceApprox" = "central" # Which approximation should be used for calculating the gradients in the gradientModel. central is recommended
    )
  )
}



