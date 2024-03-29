% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regCtsem.R
\name{exact_regCtsem}
\alias{exact_regCtsem}
\title{exact_regCtsem}
\usage{
exact_regCtsem(
  cpptsemObject,
  dataset,
  regIndicators,
  targetVector,
  lambdas,
  lambdasAutoLength,
  lambdasAutoCurve,
  penalty,
  adaptiveLassoWeights,
  returnFitIndices,
  BICWithNAndT,
  Tpoints,
  cvSampleCpptsemObject,
  optimizer,
  objective,
  sparseParameters,
  stepSize,
  lineSearch,
  c1,
  c2,
  sig,
  gam,
  initialHessianApproximation,
  maxIter_out,
  maxIter_in,
  maxIter_line,
  eps_out,
  eps_in,
  eps_WW,
  eta,
  stepsizeMin,
  stepsizeMax,
  GISTLinesearchCriterion,
  GISTNonMonotoneNBack,
  break_outer,
  scaleLambdaWithN,
  approxFirst,
  numStart,
  nMultistart,
  controlApproxOptimizer,
  verbose
)
}
\arguments{
\item{cpptsemObject}{Object of type cpptsem}

\item{dataset}{Please provide a data set in wide format compatible to ctsemOMX}

\item{regIndicators}{Labels for the regularized parameters (e.g. drift_eta1_eta2)}

\item{targetVector}{named vector with values towards which the parameters are regularized}

\item{lambdas}{vector of penalty values (tuning parameter). E.g., seq(0,1,.01). Alternatively, lambdas can be set to "auto". regCtsem will then compute an upper limit for lambda and test lambdasAutoLength increasing lambda values}

\item{lambdasAutoLength}{if lambdas == "auto", lambdasAutoLength will determine the number of lambdas tested.}

\item{lambdasAutoCurve}{It is often a good idea to have unequally spaced lambda steps (e.g., .01,.02,.05,1,5,20). If lambdasAutoCurve is close to 1 lambda values will be equally spaced, if lambdasAutoCurve is large lambda values will be more concentrated close to 0. See ?getCurvedLambda for more informations.}

\item{penalty}{Currently supported are ridge, lasso, and adaptiveLasso}

\item{adaptiveLassoWeights}{weights for the adaptive lasso. If auto, defaults to the inverse of unregularized parameter estimates.}

\item{returnFitIndices}{Boolean: should fit indices be returned?}

\item{BICWithNAndT}{Boolean: TRUE = Use N and T in the formula for the BIC (-2log L + log(N+T)*k, where k is the number of parameters in the model). FALSE = Use both N in the formula for the BIC (-2log L + log(N))}

\item{Tpoints}{Number of time points (used for BICWithNAndT)}

\item{cvSampleCpptsemObject}{cppstem for cross-validation}

\item{optimizer}{Either "GIST" or "GLMNET"}

\item{objective}{which objective should be used? Possible are "ML" (Maximum Likelihood) or "Kalman" (Kalman Filter)}

\item{sparseParameters}{labeled vector with parameter estimates of the most sparse model. Required for approxFirst = 3. If regValues = "auto" the sparse parameters will be computed automatically.}

\item{stepSize}{GLMNET & GIST: initial step size of the outer iteration}

\item{lineSearch}{GLMNET: String indicating which linesearch should be used. Defaults to the one described in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999-2030. https://doi.org/10.1145/2020408.2020421. Alternatively (not recommended) Wolfe conditions (lineSearch = "Wolfe") can be used in the outer iteration. Setting to "none" is also not recommended!}

\item{c1}{GLMNET: c1 constant for lineSearch. This constant controls the Armijo condition in lineSearch if lineSearch = "Wolfe"}

\item{c2}{GLMNET: c2 constant for lineSearch. This constant controls the Curvature condition in lineSearch if lineSearch = "Wolfe"}

\item{sig}{GLMNET & GIST: GLMNET: only relevant when lineSearch = 'GLMNET' | GIST: sigma value in Gong et al. (2013). Sigma controls the inner stopping criterion and must be in (0,1). Generally, a larger sigma enforce a steeper decrease in the regularized likelihood while a smaller sigma will result in faster acceptance of the inner iteration.}

\item{gam}{GLMNET when lineSearch = 'GLMNET'. Controls the gamma parameter in Yuan, G.-X., Ho, C.-H., & Lin, C.-J. (2012). An improved GLMNET for l1-regularized logistic regression. The Journal of Machine Learning Research, 13, 1999-2030. https://doi.org/10.1145/2020408.2020421. Defaults to 0.}

\item{initialHessianApproximation}{GLMNET: Which initial hessian approximation should be used? Possible are: 'identity' for an identity matrix and 'OpenMx' (here the hessian approxmiation from the mxObject is used). If the Hessian from 'OpenMx' is not positive definite, the negative Eigenvalues will be 'flipped' to positive Eigenvalues. "estimate" will estimate the Hessian using optimHess. All of these approaches work most of the time, but not always. Alternatively, a matrix can be provided which will be used as initial Hessian}

\item{maxIter_out}{GLMNET & GIST: Maximal number of outer iterations}

\item{maxIter_in}{GLMNET & GIST: Maximal number of inner iterations}

\item{maxIter_line}{GLMNET: Maximal number of iterations for the lineSearch procedure}

\item{eps_out}{GLMNET: Stopping criterion for outer iterations}

\item{eps_in}{GLMNET: Stopping criterion for inner iterations}

\item{eps_WW}{GLMNET: Stopping criterion for weak Wolfe line search. If the upper - lower bound of the interval is < epsWW, line search will be stopped and stepSize will be returned}

\item{eta}{GIST: if the current step size fails, eta will decrease the step size. Must be > 1}

\item{stepsizeMin}{GIST: Minimal acceptable step size. Must be > 0. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)}

\item{stepsizeMax}{GIST: Maximal acceptable step size. Must be > stepsizeMin. A larger number corresponds to a smaller step from one to the next iteration. All step sizes will be computed as described by Gong et al. (2013)}

\item{GISTLinesearchCriterion}{criterion for accepting a step. Possible are 'monotone' which enforces a monotone decrease in the objective function or 'non-monotone' which also accepts some increase.}

\item{GISTNonMonotoneNBack}{in case of non-monotone line search: Number of preceding regM2LL values to consider}

\item{break_outer}{criterion for breaking outer iterations of GIST. See ?controlGIST for more information}

\item{scaleLambdaWithN}{Boolean: Should the penalty value be scaled with the sample size? True is recommended as the likelihood is also sample size dependent}

\item{approxFirst}{Should approximate optimization be used first to obtain start values for exact optimization?}

\item{numStart}{Used if approxFirst = 3. regCtsem will try numStart+2 starting values (+2 because it will always try the current best and the parameters provided in sparseParameters)}

\item{nMultistart}{controls how many different starting values are tried when estimating lambda_max}

\item{controlApproxOptimizer}{settings passed to optimx or Rsolnp}

\item{verbose}{0 (default), 1 for convergence plot, 2 for parameter convergence plot and line search progress}
}
\description{
creates a regCtsem object for exact optimization
}
\details{
NOTE: Function located in file regCtsem.R
}
\author{
Jannik Orzek
}
\keyword{internal}
