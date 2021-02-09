#' getVariancesInParameterEstimates
#'
#' computes the variances from the parameterEstimates matrix
#' @param mxObject mxObject
#' @param parameterEstimates parameter estimates from regularized ctsem
#' @export
getVariancesInParameterEstimates <- function(mxObject, parameterEstimates){
  if(!is.matrix(parameterEstimates)){
    stop("parameterEstimates has to be of class matrix")
  }

  latentNames <- diag(mxObject$DRIFT$labels)
  latentNames <- sub(x = latentNames, pattern = "drift_", replacement = "")
  manifestNames <- rownames(mxObject$LAMBDA$values)

  for(regValue in 1:ncol(parameterEstimates)){
    tempMxObject <- mxObject
    tempMxObject <- OpenMx::omxSetParameters(model = tempMxObject,
                                             labels = rownames(parameterEstimates),
                                             values = parameterEstimates[,regValue]
                                             )

    # T0VAR
    if(any(mxObject$T0VARbase$free)){
      T0VAR <- regCtsem::getVarianceFromVarianceBase2(varianceBaseValues = tempMxObject$T0VARbase$values)

      if(regValue == 1){
        T0VARBaseLabels <- tempMxObject$T0VARbase$labels[!is.na(tempMxObject$T0VARbase$labels)]
        T0VARLabels <- paste0("T0VAR_",rep(latentNames, each = length(latentNames)),
                              "_",
                              rep(latentNames, length(latentNames)))
        T0VARLabels <- matrix(T0VARLabels,
                              nrow = length(latentNames),
                              ncol = length(latentNames),
                              byrow = TRUE)

        T0VARs <- OpenMx::cvectorize(T0VAR)
      }else{
        T0VARs <- cbind(T0VARs, OpenMx::cvectorize(T0VAR))
      }

    }

    # DIFFUSION

    if(any(mxObject$DIFFUSIONbase$free)){
      DIFFUSION <- regCtsem::getVarianceFromVarianceBase2(varianceBaseValues = tempMxObject$DIFFUSIONbase$values)

      if(regValue == 1){
        DIFFUSIONBaseLabels <- tempMxObject$DIFFUSIONbase$labels[!is.na(tempMxObject$DIFFUSIONbase$labels)]
        DIFFUSIONLabels <- paste0("DIFFUSION_",rep(latentNames, each = length(latentNames)),
                              "_",
                              rep(latentNames, length(latentNames)))
        DIFFUSIONLabels <- matrix(DIFFUSIONLabels,
                              nrow = length(latentNames),
                              ncol = length(latentNames),
                              byrow = TRUE)

        DIFFUSIONs <- OpenMx::cvectorize(DIFFUSION)
      }else{
        DIFFUSIONs <- cbind(DIFFUSIONs, OpenMx::cvectorize(DIFFUSION))
      }

    }

    # MANIFESTVAR

    if(any(mxObject$MANIFESTVARbase$free)){
      MANIFESTVAR <- regCtsem::getVarianceFromVarianceBase2(varianceBaseValues = tempMxObject$MANIFESTVARbase$values)

      if(regValue == 1){
        MANIFESTVARBaseLabels <- tempMxObject$MANIFESTVARbase$labels[!is.na(tempMxObject$MANIFESTVARbase$labels)]
        MANIFESTVARLabels <- paste0("MANIFESTVAR_",rep(manifestNames, each = length(manifestNames)),
                                  "_",
                                  rep(manifestNames, length(manifestNames)))
        MANIFESTVARLabels <- matrix(MANIFESTVARLabels,
                                  nrow = length(manifestNames),
                                  ncol = length(manifestNames),
                                  byrow = TRUE)

        MANIFESTVARs <- OpenMx::cvectorize(MANIFESTVAR)
      }else{
        MANIFESTVARs <- cbind(MANIFESTVARs, OpenMx::cvectorize(MANIFESTVAR))
      }

    }

    if(any(grepl("asymDIFFUSIONalg", mxObject$T0VAR$labels))){
      # in this case: T0VAR is set to stationarity
      DRIFTHATCH <- tempMxObject$DRIFT$values %x% tempMxObject$II$values + tempMxObject$II$values %x% tempMxObject$DRIFT$values
      asymDIFFUSIONalg <- -solve(DRIFTHATCH) %*% cvectorize(DIFFUSION)
      T0VAR <- matrix(asymDIFFUSIONalg, nrow = length(latentNames), byrow = F)

      if(regValue == 1){
        T0VARBaseLabels <- c()
        T0VARLabels <- paste0("T0VAR_",rep(latentNames, each = length(latentNames)),
                                    "_",
                                    rep(latentNames, length(latentNames)))
        T0VARLabels <- matrix(T0VARLabels,
                                    nrow = length(latentNames),
                                    ncol = length(latentNames),
                                    byrow = TRUE)

        T0VARs <- OpenMx::cvectorize(T0VAR)
      }else{
        T0VARs <- cbind(T0VARs, OpenMx::cvectorize(T0VAR))
      }

    }

  }
  # replace values in parameterEstimates
  # T0VAR
  if(any(mxObject$T0VARbase$free)){
    rownames(T0VARs) <- OpenMx::cvectorize(T0VARLabels)
    parameterEstimates <- parameterEstimates[!(rownames(parameterEstimates) %in% T0VARBaseLabels),]
    if(!is.matrix(parameterEstimates)){
      parLabels <- names(parameterEstimates)
      parameterEstimates <- matrix(parameterEstimates, nrow = length(parLabels))
      rownames(parameterEstimates) <- parLabels
    }
    parameterEstimates <- rbind(parameterEstimates, T0VARs)
  }
  if(any(grepl("asymDIFFUSIONalg", mxObject$T0VAR$labels))){
    rownames(T0VARs) <- OpenMx::cvectorize(T0VARLabels)
    if(!is.matrix(parameterEstimates)){
      parLabels <- names(parameterEstimates)
      parameterEstimates <- matrix(parameterEstimates, nrow = length(parLabels))
      rownames(parameterEstimates) <- parLabels
    }
    parameterEstimates <- rbind(parameterEstimates, T0VARs)
  }

  # DIFFUSION
  if(any(mxObject$DIFFUSIONbase$free)){
    rownames(DIFFUSIONs) <- OpenMx::cvectorize(DIFFUSIONLabels)
    parameterEstimates <- parameterEstimates[!(rownames(parameterEstimates) %in% DIFFUSIONBaseLabels),]
    if(!is.matrix(parameterEstimates)){
      parLabels <- names(parameterEstimates)
      parameterEstimates <- matrix(parameterEstimates, nrow = length(parLabels))
      rownames(parameterEstimates) <- parLabels
    }
    parameterEstimates <- rbind(parameterEstimates, DIFFUSIONs)
  }

  # MANIFESTVAR
  if(any(mxObject$MANIFESTVARbase$free)){
    rownames(MANIFESTVARs) <- OpenMx::cvectorize(MANIFESTVARLabels)
    parameterEstimates <- parameterEstimates[!(rownames(parameterEstimates) %in% MANIFESTVARBaseLabels),]
    if(!is.matrix(parameterEstimates)){
      parLabels <- names(parameterEstimates)
      parameterEstimates <- matrix(parameterEstimates, nrow = length(parLabels))
      rownames(parameterEstimates) <- parLabels
    }
    parameterEstimates <- rbind(parameterEstimates, MANIFESTVARs)
  }

  return(parameterEstimates)

}



