#' exact_getT0VAR
#'
#' computes the updated T0VAR given the old parameter values in an mxObject and the updates to the parameters in a vector d
#'
#' @param mxObject mxObject with old parameter values
#' @param d vector with updates to parameter values
#' @param stationarity set to TRUE, if stationaritiy is assumed for the T0VAR
#' @export
exact_getT0VAR <- function(mxObject, d){
  stationarity <- !("T0VAR" %in% names(mxObject$algebras))
  if(stationarity){
    stop("Not yet implemented for stationarity = 'T0VAR'. Changes in derivative necessary")

    T0Varmodel <- mxModel(
      # with free parameters
      mxObject$DRIFT,
      mxObject$DIFFUSIONbase,

      # algebras
      mxObject$T0VAR,
      mxObject$asymDIFFUSIONalg,
      mxObject$DRIFTHATCH,
      mxObject$DIFFUSION,
      mxObject$DIFFUSIONchol,
      mxObject$II)

    # set Parameters
    d_subset <- names(OpenMx::omxGetParameters(T0Varmodel))
    T0Varmodel_new <- OpenMx::omxSetParameters(T0Varmodel, labels = d_subset, values = OpenMx::omxGetParameters(T0Varmodel)+d[d_subset,])

    T0Varmodel_new.fit <- OpenMx::mxRun(T0Varmodel_new, silent = TRUE)

    return(T0Varmodel_new.fit$T0VAR$values)
  }else{
    if(!any(mxObject$T0VARbase$free)){
      return(mxObject$T0VAR$result)
    }
    d_subset <- unique(OpenMx::cvectorize(mxObject$T0VARbase$labels[!is.na(mxObject$T0VARbase$labels)])) # rownames(d)[grepl("T0var", rownames(d))]
    T0VARbaseVal <- mxObject$T0VARbase$values
    T0VARbaseLab <- mxObject$T0VARbase$labels
    for(d_sub in d_subset){
      T0VARbaseVal[T0VARbaseLab == d_sub & !is.na(T0VARbaseLab)] <- T0VARbaseVal[T0VARbaseLab == d_sub & !is.na(T0VARbaseLab)] + d[d_sub,]
    }
    T0VARchol <- OpenMx::vec2diag(exp(OpenMx::diag2vec(T0VARbaseVal))) + T0VARbaseVal - OpenMx::vec2diag(OpenMx::diag2vec(T0VARbaseVal))
    T0VAR <- T0VARchol %*% t(T0VARchol)

    return(T0VAR)

  }

}


