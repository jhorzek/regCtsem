prepareDiscreteElements <- function(mxObject, ctMatrices,
                                    nlatent, nmanifest, Tpoints, dT){

  retList <- list()

  # identify unique dTs. This will later on be used to compute the corresponding discrete time values just once
  dTUnique <- unique(dT)
  uniqueNumbers <- seq_len(length(dTUnique))
  discreteNumbering <- rep(NA, length(dT))
  for(i in seq_len(length(dTUnique))){
    discreteNumbering[dT == dTUnique[i]] <- uniqueNumbers[i]
  }

  ctParameterNames <- names(ctMatrices)

  if("DRIFT" %in% ctParameterNames){
    discreteDRIFTNames <- paste0("discreteDRIFT_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteDRIFTUnique <- vector("list", length = length(unique(discreteDRIFTNames))+2)
    names(discreteDRIFTUnique) <- c("labels", "dT", unique(discreteDRIFTNames))
    discreteDRIFTUnique$labels <- unique(discreteDRIFTNames)
    discreteDRIFTUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteDRIFTNames)))){
      discreteDRIFTUnique[[unique(discreteDRIFTNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }
    retList[["discreteDRIFTUnique"]] <- discreteDRIFTUnique
  }

  if("T0TRAITEFFECT" %in% ctParameterNames){
    T0TRAITEFFECTName <- "T0TRAITEFFECT"

    discreteTRAITNames <- paste0("discreteTRAIT_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteTRAITUnique <- vector("list", length = length(unique(discreteTRAITNames))+2)
    names(discreteTRAITUnique) <- c("labels", "dT", unique(discreteTRAITNames))
    discreteTRAITUnique$labels <- unique(discreteTRAITNames)
    discreteTRAITUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteTRAITNames)))){
      discreteTRAITUnique[[unique(discreteTRAITNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }
    retList[["discreteTRAITUnique"]] <- discreteTRAITUnique
  }

  if("DIFFUSIONbase" %in% ctParameterNames){
    discreteDIFFUSIONNames <- paste0("discreteDIFFUSION_", discreteNumbering)
    DRIFTHASHExponentialNames <- paste0("DRIFTHASHExponential_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteDIFFUSIONUnique <- vector("list", length = length(unique(discreteDIFFUSIONNames))+2)
    names(discreteDIFFUSIONUnique) <- c("labels", "dT", unique(discreteDIFFUSIONNames))
    discreteDIFFUSIONUnique$labels <- unique(discreteDIFFUSIONNames)
    discreteDIFFUSIONUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteDIFFUSIONNames)))){
      discreteDIFFUSIONUnique[[unique(discreteDIFFUSIONNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }

    DRIFTHASHExponentialUnique <- vector("list", length = length(unique(DRIFTHASHExponentialNames))+2)
    names(DRIFTHASHExponentialUnique) <- c("labels", "dT", unique(DRIFTHASHExponentialNames))
    DRIFTHASHExponentialUnique$labels <- unique(DRIFTHASHExponentialNames)
    DRIFTHASHExponentialUnique$dT <- dTUnique
    for(i in seq_len(length(unique(DRIFTHASHExponentialNames)))){
      DRIFTHASHExponentialUnique[[unique(DRIFTHASHExponentialNames)[i]]] <- matrix(0, nrow = nlatent, ncol = nlatent)
    }

    retList[["discreteDIFFUSIONUnique"]] = discreteDIFFUSIONUnique
    retList[["DRIFTHASHExponentialUnique"]] = DRIFTHASHExponentialUnique
  }

  if("CINT" %in% ctParameterNames){
    discreteCINTNames <- paste0("discreteCINT_", discreteNumbering)
    # combine numbered discrete time parameters and time intervals
    discreteCINTUnique <- vector("list", length = length(unique(discreteCINTNames))+2)
    names(discreteCINTUnique) <- c("labels", "dT", unique(discreteCINTNames))
    discreteCINTUnique$labels <- unique(discreteCINTNames)
    discreteCINTUnique$dT <- dTUnique
    for(i in seq_len(length(unique(discreteCINTNames)))){
      discreteCINTUnique[[unique(discreteCINTNames)[i]]] <- matrix(0, nrow = nlatent, ncol = 1)
    }

    retList[["discreteCINTUnique"]] <- discreteCINTUnique
  }


  return(retList)
}

