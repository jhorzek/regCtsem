#' checkSetup
#'
#' internal checks
#'
#' @param argsIn list with arguments
#' @author Jannik Orzek
#' @export
checkSetup <- function(argsIn){
  if(is.null(argsIn$ctsemObject) & is.null(argsIn$mxObject)){
    stop("Both ctsemObject and mxObect are missing. You need to provide at least one")
  }

  if(any(class(argsIn$ctsemObject)=="ctsemFit")){
    # check if ctsemObject is estimated with Kalman
    if(argsIn$ctsemObject$ctfitargs$objective == "Kalman"){
      stop("It seems like the provided ctsemObject was fitted with Kalman filter. To use the Kalman filter, provide the object of type ctsemInit from ctModel instead of the fitted model. Set the objective = 'Kalman' and provide a dataset")
    }
  }
  if(!is.null(argsIn$mxObject) & any(class(argsIn$mxObject)=="MxModel")){
    # check if ctsemObject is estimated with Kalman
    if(any(class(argsIn$mxObject$expectation) == "MxExpectationStateSpace") &  !any(class(argsIn$ctsemObject)=="ctsemInit")){
      stop("It seems like the provided mxObject was fitted with Kalman filter. To use the Kalman filter, provide the object of type ctsemInit from ctModel instead of the fitted model. Set the objective = 'Kalman' and provide a dataset")
    }
  }

  if(!(argsIn$optimization == "exact" | argsIn$optimization == "approx")){
    stop(paste("Optimization was set to", optimization, "however only exact or approx are supported."))
  }
  if(argsIn$optimization == "approx" & argsIn$extraTries == 1){
    warning(paste("Approximate optimization often requires multiple tries to find the optimal solution. Use extraTries to automatically try different statring values (e.g., extraTries = 5)."))
  }
  if(argsIn$cores > 1 & argsIn$verbose>0){
    stop("verbose > 0 only possible for single core execution. Set cores = 1 or verbose = 0.")
  }

  if(tolower(argsIn$lineSearch) == "armijo"){
    stop("lineSearch = 'armijo' is deprecated. Use lineSearch = 'Wolfe'.")
  }

  if(tolower(argsIn$objective) == "kalman"){
    if(is.null(argsIn$ctsemObject)){stop("Kalman filter requires a ctsemObject. Set up the model in ctsemOMX first")}
    if(!any(class(argsIn$ctsemObject) == "ctsemInit")){stop("ctsemObject has to be of class ctsemInit. You can get the ctsemInit object from ctModel.")}
    if(is.null(argsIn$dataset)){stop("No dataset was provided")}
  }

  if(tolower(argsIn$penalty) == "ridge"){
    if(argsIn$standardizeDrift){
      stop("Automatic drift standardization not supported for ridge regularization")
    }
    if(argsIn$optimization == "exact"){
      stop("Use optimization = approx for ridge regularization. This is confusing, but approx is the right optimization for ridge regularized SEM. The approx and exact differentiation is only relevant for lasso and adaptivelasso.")
    }
  }

  if(tolower(argsIn$penalty) == "adaptivelasso"){
    if(is.null(argsIn$adaptiveLassoWeights)){
      stop("Adaptive Lasso requires the adaptiveLassoWeights to be specified. auto uses the inverse of the absolute values of unregularized parameter estimates. This is only recommended with standardizeDrift = FALSE")
    }
    if(argsIn$standardizeDrift){
      stop("Combination of automatic standardization and adaptive lasso is not implemented. Standardization is a special variant of the adaptive lasso. Use penalty = 'lasso' or standardizeDrift = FALSE.")
    }
  }
  #if(any(argsIn$regValues == "auto") & (tolower(argsIn$penalty) == "adaptivelasso" | argsIn$standardizeDrift)){
  #  stop("regValues = 'auto' currently not supported for adative lasso or automatic standardization of drift parameters.")
  #}
}


