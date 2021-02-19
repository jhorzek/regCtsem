prepareAMatrix <- function(mxObject, ctMatrices, nlatent, nmanifest, Tpoints, dT){

  retList <- list()
  # AParameterIndicators will be used to identify the segments in the
  # A matrix which belong to the discreteDRIFTs discreteTRAITs and LAMBDAS

  AParameterIndicators <- data.frame()
  AParameterIndicatorsLabels <- c()

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

    AParameterIndicatorsLabels <- c(AParameterIndicatorsLabels, discreteDRIFTNames)
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

    AParameterIndicatorsLabels <- c(AParameterIndicatorsLabels, discreteTRAITNames)
  }

  if("LAMBDA" %in% ctParameterNames){
    LAMBDANames <- rep("LAMBDA", Tpoints)

    AParameterIndicatorsLabels <- c(AParameterIndicatorsLabels, LAMBDANames)
  }



  AParameterIndicators <- data.frame("label" = AParameterIndicatorsLabels,
                                     "row_start" = NA,
                                     "row_end" = NA,
                                     "col_start" = NA,
                                     "col_end" = NA
  )

  # A matrix
  #AInitializer <- matrix(0, nrow = (nlatent)*Tpoints +nlatent# DRIFTPart and TRAITPart
  #                       + nmanifest*Tpoints # MANIFESTPart
  #                       ,
  #                       ncol = (nlatent)*Tpoints +nlatent # eta <- eta and eta <- trait
  #                       + nmanifest*Tpoints # y <- eta
  #)

  AInitializer <- mxObject$A$values

  ## segment A

  counterOffset <- 0

  # segment DRIFTs
  if("DRIFT" %in% ctParameterNames){
    discreteDRIFTRowStart <- seq(from = nlatent+1, # skip initial time point
                                 to = (Tpoints-1)*nlatent+1,
                                 by = nlatent)
    discreteDRIFTColStart <- seq(from = 1,
                                 to = (Tpoints-2)*nlatent+1,
                                 by = nlatent)

    for(discreteDRIFT in seq_len(length(discreteDRIFTNames)) + counterOffset){
      if(!(AParameterIndicators[discreteDRIFT,"label"] == discreteDRIFTNames[discreteDRIFT-counterOffset])){
        stop("Wrong element selected in segment DRIFTs")
      }
      AParameterIndicators[discreteDRIFT,"row_start"] <- discreteDRIFTRowStart[discreteDRIFT-counterOffset]
      AParameterIndicators[discreteDRIFT,"row_end"] <- discreteDRIFTRowStart[discreteDRIFT-counterOffset] + nlatent -1
      AParameterIndicators[discreteDRIFT,"col_start"] <- discreteDRIFTColStart[discreteDRIFT-counterOffset]
      AParameterIndicators[discreteDRIFT,"col_end"] <- discreteDRIFTColStart[discreteDRIFT-counterOffset] + nlatent -1
    }
    counterOffset <- length(discreteDRIFTNames)

    retList[["discreteDRIFTUnique"]] <- discreteDRIFTUnique
  }

  # segment TRAITs
  if("T0TRAITEFFECT" %in% ctParameterNames){
    addTrait <- nlatent
    if(any(mxObject$T0TRAITEFFECT$free)){
      stop("free T0TRAITEFFECT detected. This is not implemented.")
    }
    AInitializer[1:nlatent, (nlatent*Tpoints+1) : (nlatent*Tpoints+nlatent)] <- ctMatrices$T0TRAITEFFECT$values

    discreteTRAITRowStart <- seq(from = nlatent+1, # skip initial time point
                                 to = (Tpoints-1)*nlatent+1,
                                 by = nlatent)
    discreteTRAITColStart <- rep(nlatent*Tpoints+1, Tpoints-1)

    for(discreteTRAIT in seq_len(length(discreteTRAITNames)) + counterOffset){
      if(!(AParameterIndicators[discreteTRAIT,"label"] == discreteTRAITNames[discreteTRAIT-counterOffset])){
        stop("Wrong element selected in segment TRAITs")
      }
      AParameterIndicators[discreteTRAIT,"row_start"] <- discreteTRAITRowStart[discreteTRAIT-counterOffset]
      AParameterIndicators[discreteTRAIT,"row_end"] <- discreteTRAITRowStart[discreteTRAIT-counterOffset] + nlatent -1
      AParameterIndicators[discreteTRAIT,"col_start"] <- discreteTRAITColStart[discreteTRAIT-counterOffset]
      AParameterIndicators[discreteTRAIT,"col_end"] <- discreteTRAITColStart[discreteTRAIT-counterOffset] + nlatent -1
    }
    counterOffset <- counterOffset + length(discreteTRAITNames)

    retList[["discreteTRAITUnique"]] <- discreteTRAITUnique
  }else{
    addTrait <- 0
  }

  if("LAMBDA" %in% ctParameterNames){
    LAMBDARowStart <- seq(from = Tpoints*(nlatent) + addTrait + 1,
                          to = Tpoints*(nlatent) + addTrait + (Tpoints-1)*nmanifest + 1,
                          by = nmanifest)
    LAMBDAColStart <- seq(from = 1,
                          to = (Tpoints-1)*(nlatent) + 1,
                          by = nlatent)

    for(LAMBDA in seq_len(length(LAMBDANames)) + counterOffset){
      if(!(AParameterIndicators[LAMBDA,"label"] == LAMBDANames[LAMBDA-counterOffset])){
        stop("Wrong element selected in segment TRAITs")
      }
      AParameterIndicators[LAMBDA,"row_start"] <- LAMBDARowStart[LAMBDA-counterOffset]
      AParameterIndicators[LAMBDA,"row_end"] <- LAMBDARowStart[LAMBDA-counterOffset] + nmanifest -1
      AParameterIndicators[LAMBDA,"col_start"] <- LAMBDAColStart[LAMBDA-counterOffset]
      AParameterIndicators[LAMBDA,"col_end"] <- LAMBDAColStart[LAMBDA-counterOffset] + nlatent -1
    }
  }

  # cpp indices start at 0
  cppAParameterIndicators <- AParameterIndicators
  cppAParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] <- cppAParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] - 1

  retList[["AParameterIndicators"]] = AParameterIndicators
  retList[["cppAParameterIndicators"]] = cppAParameterIndicators
  retList[["AInitializer"]] = AInitializer

  return(retList)
}

