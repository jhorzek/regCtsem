#' controlApprox
#'
#' The following arguments can be used to adjust the approximate optimization
#'
#' @param forceCpptsem should cpptsem be enforced even if results differ from ctsem? Sometimes differences between cpptsem and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
#' @param epsilon epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
#' @param zeroThresh threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
#' @param controlApproxOptimizer settings passed to the optimizer in approximate optimization. Currently, Rsolnp and optimx are supported. See ?controlOptimx and ?controlSolnp for details on the lists passed to controlApprox
#' @export
controlApprox <- function(forceCpptsem = FALSE, # should the C++ translation be enforced even if results differ from ctsem? Sometimes differences between the C++ implementation and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
                          epsilon = .001, # epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
                          zeroThresh = .04, # threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
                          controlApproxOptimizer = controlRsolnp(control = list("outer.iter" = 500, "trace" = 0))){
  return(controlApprox <- list(
    "forceCpptsem" = forceCpptsem, # should the C++ translation be enforced even if results differ from ctsem? Sometimes differences between the C++ implementation and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
    "epsilon" = epsilon, # epsilon is used to transform the non-differentiable lasso penalty to a differentiable one if optimization = approx
    "zeroThresh" = zeroThresh, # threshold below which parameters will be evaluated as == 0 in lasso regularization if optimization = approx
    "controlApproxOptimizer" = controlApproxOptimizer
  )
  )
}

#' controlOptimx
#'
#' list with settings used for optimization with optimx
#' @param package set to "optimx"
#' @param failureReturns what will the fitting function return if the current points are impossible? Depends on the method used
#' @param hess Should the Hessian be computed at the solution
#' @param lower lower bounds for paramters
#' @param upper upper bounds for parameters
#' @param method see ?optimx
#' @param itnmax maximal number of iterations
#' @param control control passed to optimx
#' @export
controlOptimx <- function(package = "optimx",
                          failureReturns = .Machine$double.xmax/2,
                          hess = NULL,
                          lower = -Inf,
                          upper = Inf,
                          method = "L-BFGS-B",
                          hessian = FALSE,
                          itnmax = 200,
                          control = list("dowarn" = FALSE,
                                         "kkt" = TRUE,
                                         "maxit" = 200)){
  return(list(
    "package" = package, # package from where the optimizer is taken
    "failureReturns" = failureReturns, # not part of the optimx documentation; controls the value that will be returned if the current paramter values are impossible
    "hess" = hess,
    "lower" = lower,
    "upper" = upper,
    "method" = method,
    "hessian" = hessian,
    "itnmax" = itnmax,
    "control" = control)
  )
}

#' controlRsolnp
#'
#' list with settings used for optimization with Rsolnp
#' @param package set to "Rsolnp"
#' @param failureReturns what will the fitting function return if the current points are impossible?
#' @param eqfun Equality constraints function. See ?Rsolnp::solnp
#' @param eqB Equality constraints. See ?Rsolnp::solnp
#' @param ineqfun Inequality constraints function. See ?Rsolnp::solnp
#' @param ineqLB Inequality constraints lower bound. See ?Rsolnp::solnp
#' @param ineqUB Inequality constraints upper bound. See ?Rsolnp::solnp
#' @param LB Lower bound. See ?Rsolnp::solnp
#' @param UB Upper bound. See ?Rsolnp::solnp
#' @param control control passed to Rsolnp
#' @export
controlRsolnp <- function(package = "Rsolnp",
                          failureReturns = .Machine$double.xmax/2,
                          eqfun = NULL, eqB = NULL, ineqfun = NULL, ineqLB = NULL,
                          ineqUB = NULL, LB = NULL, UB = NULL,
                          control = list("outer.iter" = 50, "trace" = 0)
){
  return(list(
    "package" = package, # package from where the optimizer is taken
    "failureReturns" = failureReturns, # not part of the optimx documentation; controls the value that will be returned if the current paramter values are impossible
    "eqfun" = eqfun,
    "eqB" = eqB,
    "ineqfun" = ineqfun,
    "ineqLB" = ineqLB,
    "LB" = LB,
    "UB" = UB,
    "control" = control)
  )
}


#' controlGIST
#'
#' The following arguments can be used to adjust the GIST optimization
#'
#' @param forceCpptsem should cpptsem be enforced even if results differ from ctsem? Sometimes differences between cpptsem and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
#' @param stepSize initial step size of the outer iteration
#' @param sig sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
#' @param maxIter_out Maximal number of outer iterations
#' @param maxIter_in Maximal number of inner iterations
#' @param break_outer Stopping criterion for outer iterations; has to be a named value. By default the change in fit is used, with c("fitChange" = 1e-5) meaning that this change should be smaller than 1e-5. Additionally, a relative change in parameters is used as breaking criterion c("parameterChange" = .00001). Alternatively (name: gradient), a relative first-order condition is checked, where the maximum absolute value of the gradients is compared to break_outer (see https://de.mathworks.com/help/optim/ug/first-order-optimality-measure.html). Example: c("gradient" = "max(max(abs(startingValues))*.001, .001)") . Alternatively, an absolute tolerance can be passed to the function (e.g., break_outer = c("gradient" = .0001)).
#' @param eta if the current step size fails, eta will decrease the step size. Must be > 1
#' @param stepsizeMin Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param stepsizeMax Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
#' @param GISTLinesearchCriterion criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
#' @param GISTNonMonotoneNBack in case of non-monotone line search: Number of preceding regM2LL values to consider
#' @param approxFirst Should approximate optimization be used first to obtain start values for exact optimization?
#' @param numStart Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param controlApproxOptimizer settings passed to the optimizer in approximate optimization. Currently, Rsolnp and optimx are supported. See ?controlOptimx and ?controlSolnp for details on the lists passed to controlApprox
#' @export
controlGIST <- function(forceCpptsem = FALSE, # should the C++ translation be enforced even if results differ from ctsem? Sometimes differences between the C++ implementation and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
                        stepSize =  1, # initial step size of the outer iteration
                        sig = 10^(-5), # sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
                        maxIter_out = 100, # Maximal number of outer iterations
                        maxIter_in = 1000, # Maximal number of inner iterations
                        break_outer = c("fitChange" = 10^(-5)), # Stopping criterion for outer iterations; has to be a named value. By default the change in fit is used, with c("fitChange" = 1e-5) meaning that this change should be smaller than 1e-5. Additionally, a relative change in parameters is used as breaking criterion c("parameterChange" = .00001). Alternatively (name: gradient), a relative first-order condition is checked, where the maximum absolute value of the gradients is compared to break_outer (see https://de.mathworks.com/help/optim/ug/first-order-optimality-measure.html). Example: c("gradient" = "max(max(abs(startingValues))*.001, .001)") . Alternatively, an absolute tolerance can be passed to the function (e.g., break_outer = c("gradient" = .0001)).
                        eta = 2, # if the current step size fails, eta will decrease the step size. Must be > 1
                        stepsizeMin = 1/(10^30), # Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
                        stepsizeMax = 10^30, # Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
                        GISTLinesearchCriterion = "monotone", # criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
                        GISTNonMonotoneNBack = 5,# in case of non-monotone line search: Number of preceding regM2LL values to consider
                        approxFirst = TRUE, # Should approximate optimization be used first to obtain start values for exact optimization?
                        numStart = 0, # Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
                        controlApproxOptimizer = controlRsolnp(control = list("outer.iter" = 50, "trace" = 0))){
  return(controlGIST <- list(
    "forceCpptsem" = forceCpptsem, # should the C++ translation be enforced even if results differ from ctsem? Sometimes differences between the C++ implementation and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
    "stepSize" =  stepSize, # initial step size of the outer iteration
    "sig" = sig, # sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.
    "maxIter_out" = maxIter_out, # Maximal number of outer iterations
    "maxIter_in" = maxIter_in, # Maximal number of inner iterations
    "break_outer" = break_outer, # Stopping criterion for outer iterations; has to be a named value. By default the change in fit is used, with c("fitChange" = 1e-5) meaning that this change should be smaller than 1e-5. Additionally, a relative change in parameters is used as breaking criterion c("parameterChange" = .00001). Alternatively (name: gradient), a relative first-order condition is checked, where the maximum absolute value of the gradients is compared to break_outer (see https://de.mathworks.com/help/optim/ug/first-order-optimality-measure.html). Example: c("gradient" = "max(max(abs(startingValues))*.001, .001)") . Alternatively, an absolute tolerance can be passed to the function (e.g., break_outer = c("gradient" = .0001)).
    "eta" = eta, # if the current step size fails, eta will decrease the step size. Must be > 1
    "stepsizeMin" = stepsizeMin, # Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
    "stepsizeMax" = stepsizeMax, # Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)
    "GISTLinesearchCriterion" = GISTLinesearchCriterion, # criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.
    "GISTNonMonotoneNBack" = GISTNonMonotoneNBack,# in case of non-monotone line search: Number of preceding regM2LL values to consider
    "approxFirst" = approxFirst, # Should approximate optimization be used first to obtain start values for exact optimization?
    "numStart" = numStart, # Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
    "controlApproxOptimizer" = controlApproxOptimizer
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
#' @param initialHessianApproximation Which initial hessian approximation should be used? Possible are: 'ident' for an identity matrix and 'OpenMx' (here the hessian approxmiation from the mxObject is used). If the Hessian from 'OpenMx' is not positive definite, the negative Eigenvalues will be 'flipped' to positive Eigenvalues. This works sometimes, but not always. Alternatively, a matrix can be provided which will be used as initial Hessian
#' @param maxIter_out Maximal number of outer iterations
#' @param maxIter_in Maximal number of inner iterations
#' @param maxIter_line Maximal number of iterations for the lineSearch procedure
#' @param eps_out Stopping criterion for outer iterations
#' @param eps_in Stopping criterion for inner iterations
#' @param eps_WW Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned
#' @param approxFirst Should approximate optimization be used first to obtain start values for exact optimization?
#' @param numStart Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
#' @param controlApproxOptimizer settings passed to the optimizer in approximate optimization. Currently, Rsolnp and optimx are supported. See ?controlOptimx and ?controlSolnp for details on the lists passed to controlApprox
#' @export
controlGLMNET <- function(      tryCpptsem = TRUE, # should regCtsem try to translate the model to C++? This can speed up the computation considerably but might fail for some models
                                forceCpptsem = FALSE, # should the C++ translation be enforced even if results differ from ctsem? Sometimes differences between the C++ implementation and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
                                stepSize =  1, # initial step size of the outer iteration
                                lineSearch = "GLMNET", # String indicating which linesearch should be used. Defaults to the one described in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Alternatively (not recommended) Wolfe conditions (lineSearch = "Wolfe") can be used in the outer iteration. Setting to "none" is also not recommended!.
                                c1 = .0001, # c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
                                c2 = .9, # c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
                                sig = 10^(-5), # for lineSearch = 'GLMNET': Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 < sigma < 1
                                gam = 0, # for lineSearch = 'GLMNET': Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 <= gamma < 1
                                initialHessianApproximation = "OpenMx", # Which initial hessian approximation should be used? Possible are: 'ident' for an identity matrix and 'OpenMx' (here the hessian approxmiation from the mxObject is used). If the Hessian from 'OpenMx' is not positive definite, the negative Eigenvalues will be 'flipped' to positive Eigenvalues. This works sometimes, but not always. Alternatively, a matrix can be provided which will be used as initial Hessian
                                maxIter_out = 100, # Maximal number of outer iterations
                                maxIter_in = 1000, # Maximal number of inner iterations
                                maxIter_line = 500, # Maximal number of iterations for the lineSearch procedure
                                eps_out = .0000000001, # Stopping criterion for outer iterations
                                eps_in = .0000000001, # Stopping criterion for inner iterations
                                eps_WW = .0001, #Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned
                                approxFirst = TRUE, # Should approximate optimization be used first to obtain start values for exact optimization?
                                numStart = 0, # Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
                                controlApproxOptimizer = controlRsolnp(control = list("outer.iter" = 50, "trace" = 0))){
  return(
    controlGLMNET <- list(
      "tryCpptsem" = tryCpptsem, # should regCtsem try to translate the model to C++? This can speed up the computation considerably but might fail for some models
      "forceCpptsem" = forceCpptsem, # should the C++ translation be enforced even if results differ from ctsem? Sometimes differences between the C++ implementation and ctsem can result from problems with numerical precision which will lead to the matrix exponential of RcppArmadillo differing from the OpenMx matrix exponential. If you want to ensure the faster optimization, set to TRUE. See vignette("MatrixExponential", package = "regCtsem") for more details
      "stepSize" =  stepSize, # initial step size of the outer iteration
      "lineSearch" = "lineSearch", # String indicating which linesearch should be used. Defaults to the one described in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421. Alternatively (not recommended) Wolfe conditions (lineSearch = "Wolfe") can be used in the outer iteration. Setting to "none" is also not recommended!.
      "c1" = c1, # c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"
      "c2" = c2, # c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"
      "sig" = sig, # for lineSearch = 'GLMNET': Controls the sigma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 < sigma < 1
      "gam" = gam, # for lineSearch = 'GLMNET': Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999–2030. https://doi.org/10.1145/2020408.2020421, Equation 20. Defaults to 0. Has to be in 0 <= gamma < 1
      "initialHessianApproximation" = initialHessianApproximation, # Which initial hessian approximation should be used? Possible are: 'ident' for an identity matrix and 'OpenMx' (here the hessian approxmiation from the mxObject is used). If the Hessian from 'OpenMx' is not positive definite, the negative Eigenvalues will be 'flipped' to positive Eigenvalues. This works sometimes, but not always. Alternatively, a matrix can be provided which will be used as initial Hessian
      "maxIter_out" = maxIter_out, # Maximal number of outer iterations
      "maxIter_in" = maxIter_in, # Maximal number of inner iterations
      "maxIter_line" = maxIter_line, # Maximal number of iterations for the lineSearch procedure
      "eps_out" = eps_out, # Stopping criterion for outer iterations
      "eps_in" = eps_in, # Stopping criterion for inner iterations
      "eps_WW" = eps_WW, #Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned
      "approxFirst" = approxFirst, # Should approximate optimization be used first to obtain start values for exact optimization?
      "numStart" = numStart, # Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)
      "controlApproxOptimizer" = controlApproxOptimizer
    )
  )
}



