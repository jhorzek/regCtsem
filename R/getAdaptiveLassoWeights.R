#' getAdaptiveLassoWeights
#'
#' Computes the adaptive lasso weights
#'
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject. Provide either ctsemObject or mxObject
#' @param penalty type. Currently supported are lasso and ridge for optimization = approx and lasso for optimization = exact
#' @param adaptiveLassoWeights weights for the adaptive lasso.
#' @param standardizeDrift Boolean: Should Drift parameters be standardized automatically using T0VAR?
#' @author Jannik Orzek
#' @export
getAdaptiveLassoWeights <- function(mxObject, penalty, adaptiveLassoWeights, standardizeDrift){

  # if adaptiveLassoWeights were provided and no standardization is requested:
  if(tolower(penalty) == "adaptivelasso" & is.numeric(adaptiveLassoWeights) & (!standardizeDrift)){
    return(adaptiveLassoWeights)
  }

  # otherwise: set adaptiveLassoWeights
  ## if automatic drift standardization was requested
  if(standardizeDrift){
    # check if mxObject was optimized
    ## TODO

    # check if adaptiveLassoWeights were provided
    if(is.numeric(adaptiveLassoWeights)){ stop("standardizeDrift and provided adaptiveLassoWeights can not be combined automatically. Consider setting adaptiveLassoWeights = 'auto'") }

    # check if lasso or adaptiveLasso
    if(tolower(penalty) == "lasso"){
      thetaNames <- names(OpenMx::omxGetParameters(mxObject))
      adaptiveLassoWeights <- rep(1,length(thetaNames))
      names(adaptiveLassoWeights) <- thetaNames

      # compute standardizers
      T0VAR <- getT0VAR(mxObject)
      driftLabels <- mxObject$DRIFT$labels
      if(anyNA(driftLabels)){
        autoDriftLabels <- matrix(paste0("_autoDriftLabel_", rep(seq_len(nrow(driftLabels)), each = ncol(driftLabels)), "_", seq_len(ncol(driftLabels))),
                                  nrow = nrow(driftLabels), ncol = ncol(driftLabels), byrow = T)
        driftLabels[is.na(driftLabels)] <- autoDriftLabels[is.na(driftLabels)]
      }
      flatStandardizers <- regCtsem::getFlatStdizer(T0VAR = T0VAR, driftLabels = driftLabels)
      flatStandardizers <- flatStandardizers[rownames(flatStandardizers) %in% thetaNames,]

      for(thetaName in names(flatStandardizers)){
        adaptiveLassoWeights[thetaName] <- flatStandardizers[thetaName]*adaptiveLassoWeights[thetaName]
      }

      return(adaptiveLassoWeights)
    }

    if(tolower(penalty) == "adaptivelasso" && adaptiveLassoWeights == "auto"){
      stop("Automatic standardization and adaptive lasso weights can not be combined.")
      thetaNames <- names(OpenMx::omxGetParameters(mxObject))
      adaptiveLassoWeights <- abs(OpenMx::omxGetParameters(mxObject))^(-1)
      # compute standardizers
      T0VAR <- getT0VAR(mxObject)
      driftLabels <- mxObject$DRIFT$labels
      if(anyNA(driftLabels)){
        autoDriftLabels <- matrix(paste0("_autoDriftLabel_", rep(seq_len(nrow(driftLabels)), each = ncol(driftLabels)), "_", seq_len(ncol(driftLabels))),
                                  nrow = nrow(driftLabels), ncol = ncol(driftLabels), byrow = T)
        driftLabels[is.na(driftLabels)] <- autoDriftLabels[is.na(driftLabels)]
      }

      flatStandardizers <- getFlatStdizer(T0VAR = T0VAR, driftLabels = driftLabels)
      flatStandardizers <- flatStandardizers[rownames(flatStandardizers) %in% thetaNames,]

      for(thetaName in rownames(flatStandardizers)){
        adaptiveLassoWeights[thetaName] <- flatStandardizers[thetaName]*adaptiveLassoWeights[thetaName]
      }

      return(adaptiveLassoWeights)
    }
  }

  ## if lasso was requested, but no automatic standardization
  if(tolower(penalty) == "lasso"){
    thetaNames <- names(OpenMx::omxGetParameters(mxObject))
    adaptiveLassoWeights <- rep(1,length(thetaNames))
    names(adaptiveLassoWeights) <- thetaNames
    return(adaptiveLassoWeights)
  }

  if(tolower(penalty) == "ridge"){
    thetaNames <- names(OpenMx::omxGetParameters(mxObject))
    adaptiveLassoWeights <- rep(1,length(thetaNames))
    names(adaptiveLassoWeights) <- thetaNames
    return(adaptiveLassoWeights)
  }

  ## if adatpive lasso with automatic weigths was requested, but no automatic standardization
  if(tolower(penalty) == "adaptivelasso" && adaptiveLassoWeights == "auto"){
    adaptiveLassoWeights <- abs(OpenMx::omxGetParameters(mxObject))^(-1)
    return(adaptiveLassoWeights)
  }

  # else
  stop("Error while computing adaptive lasso weigths")

}




