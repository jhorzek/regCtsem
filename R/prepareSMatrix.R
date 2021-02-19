prepareSMatrix <- function(mxObject, ctMatrices, nlatent, nmanifest, Tpoints, dT, stationaryT0VAR){
  retList <- list()
  # SParameterIndicators will be used to identify the segments in the
  # S matrix which belong to the discreteDIFFUSION MANIFESVAR and T0VAR

  SParameterIndicators <- data.frame()
  SParameterIndicatorsLabels <- c()

  # identify unique dTs. This will later on be used to compute the corresponding discrete time values just once
  dTUnique <- unique(dT)
  uniqueNumbers <- seq_len(length(dTUnique))
  discreteNumbering <- rep(NA, length(dT))
  for(i in seq_len(length(dTUnique))){
    discreteNumbering[dT == dTUnique[i]] <- uniqueNumbers[i]
  }

  ctParameterNames <- names(ctMatrices)

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

    SParameterIndicatorsLabels <- c(SParameterIndicatorsLabels, discreteDIFFUSIONNames)
  }

  if("MANIFESTVARbase"%in% ctParameterNames){
    MANIFESTVARNames <- rep("MANIFESTVAR", Tpoints)
    SParameterIndicatorsLabels <- c(SParameterIndicatorsLabels, MANIFESTVARNames)
  }
  if("T0VARbase" %in% ctParameterNames | stationaryT0VAR){
    T0VARNames <- "T0VAR"
    SParameterIndicatorsLabels <- c(SParameterIndicatorsLabels, T0VARNames)
  }

  if("TRAITVARbase" %in% ctParameterNames){
    TRAITVARNames <- "TRAITVAR"
    SParameterIndicatorsLabels <- c(SParameterIndicatorsLabels, TRAITVARNames)
  }

  SParameterIndicators <- data.frame("label" = SParameterIndicatorsLabels,
                                     "row_start" = NA,
                                     "row_end" = NA,
                                     "col_start" = NA,
                                     "col_end" = NA
  )

  ## Initialize S

  #SInitializer <- matrix(0, nrow = (nlatent)*Tpoints +nlatent # DRIFTPart and TRAITPart
  #                       + nmanifest*Tpoints # MANIFESTPart
  #                       ,
  #                       ncol = (nlatent)*Tpoints + nlatent # eta <-> eta and eta <- trait
  #                       + nmanifest*Tpoints # y <-> eta
  #)
  SInitializer <- mxObject$S$values

  ## segment S

  counterOffset <- 0

  if("DIFFUSIONbase" %in% ctParameterNames){
    # segment discreteDIFFUSION
    discreteDIFFUSIONsRowStart <- seq(from = nlatent+1, # skip initial time point
                                      to = (Tpoints-1)*nlatent+1,
                                      by = nlatent)
    discreteDIFFUSIONsColStart <- seq(from = nlatent+1, # skip initial time point
                                      to = (Tpoints-1)*nlatent+1,
                                      by = nlatent)

    for(discreteDIFFUSIONs in seq_len(length(discreteDIFFUSIONNames)) + counterOffset){
      if(!(SParameterIndicators[discreteDIFFUSIONs,"label"] == discreteDIFFUSIONNames[discreteDIFFUSIONs-counterOffset])){
        stop("Wrong element selected in segment discreteDIFFUSIONs")
      }
      SParameterIndicators[discreteDIFFUSIONs,"row_start"] <- discreteDIFFUSIONsRowStart[discreteDIFFUSIONs-counterOffset]
      SParameterIndicators[discreteDIFFUSIONs,"row_end"] <- discreteDIFFUSIONsRowStart[discreteDIFFUSIONs-counterOffset] + nlatent -1
      SParameterIndicators[discreteDIFFUSIONs,"col_start"] <- discreteDIFFUSIONsColStart[discreteDIFFUSIONs-counterOffset]
      SParameterIndicators[discreteDIFFUSIONs,"col_end"] <- discreteDIFFUSIONsColStart[discreteDIFFUSIONs-counterOffset] + nlatent -1
    }
    counterOffset <- length(discreteDIFFUSIONNames)

    retList[["discreteDIFFUSIONUnique"]] = discreteDIFFUSIONUnique
    retList[["DRIFTHASHExponentialUnique"]] = DRIFTHASHExponentialUnique
  }

  if("MANIFESTVARbase" %in% ctParameterNames){
    if("TRAITVARbase" %in% ctParameterNames){
      traitVarOffset <- nlatent
    }else{
        traitVarOffset <- 0
        }

    # segment MANIFESTVARs
    MANIFESTVARRowStart <- seq(from = Tpoints*(nlatent) + traitVarOffset + 1,
                               to = Tpoints*(nlatent) + traitVarOffset + (Tpoints-1)*nmanifest + 1,
                               by = nmanifest)
    MANIFESTVARColStart <- seq(from = Tpoints*(nlatent) + traitVarOffset + 1,
                               to = Tpoints*(nlatent) + traitVarOffset + (Tpoints-1)*nmanifest + 1,
                               by = nmanifest)

    for(MANIFESTVAR in seq_len(length(MANIFESTVARNames)) + counterOffset){
      if(!(SParameterIndicators[MANIFESTVAR,"label"] == MANIFESTVARNames[MANIFESTVAR-counterOffset])){
        stop("Wrong element selected in segment TRAITs")
      }
      SParameterIndicators[MANIFESTVAR,"row_start"] <- MANIFESTVARRowStart[MANIFESTVAR-counterOffset]
      SParameterIndicators[MANIFESTVAR,"row_end"] <- MANIFESTVARRowStart[MANIFESTVAR-counterOffset] + nmanifest -1
      SParameterIndicators[MANIFESTVAR,"col_start"] <- MANIFESTVARColStart[MANIFESTVAR-counterOffset]
      SParameterIndicators[MANIFESTVAR,"col_end"] <- MANIFESTVARColStart[MANIFESTVAR-counterOffset] + nmanifest -1
    }

    counterOffset <- counterOffset + length(MANIFESTVARNames)
  }

  if("T0VARbase" %in% ctParameterNames | stationaryT0VAR){
    # segment T0VAR
    SParameterIndicators[counterOffset+1,"row_start"] <- 1
    SParameterIndicators[counterOffset+1,"row_end"] <- nlatent
    SParameterIndicators[counterOffset+1,"col_start"] <- 1
    SParameterIndicators[counterOffset+1,"col_end"] <- nlatent
  }

  if("TRAITVARbase" %in% ctParameterNames){
    # segment TRAITVAR
    SParameterIndicators[counterOffset+2,"row_start"] <- nlatent*Tpoints + 1
    SParameterIndicators[counterOffset+2,"row_end"] <- nlatent*Tpoints + nlatent
    SParameterIndicators[counterOffset+2,"col_start"] <- nlatent*Tpoints + 1
    SParameterIndicators[counterOffset+2,"col_end"] <- nlatent*Tpoints + nlatent
  }

  # cpp indices start at 0
  cppSParameterIndicators <- SParameterIndicators
  cppSParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] <- cppSParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] - 1

  retList[["SParameterIndicators"]] = SParameterIndicators
  retList[["cppSParameterIndicators"]] = cppSParameterIndicators
  retList[["SInitializer"]] = SInitializer

  return(retList)

}

