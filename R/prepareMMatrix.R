prepareMMatrix <- function(mxObject, ctMatrices, nlatent, nmanifest, Tpoints, dT, stationaryT0MEANS){

  retList <- list()
  # MParameterIndicators will be used to identify the segments in the
  # M matrix which belong to the CINT and MANIFESTMEANS
  MParameterIndicators <- data.frame()
  MParameterIndicatorsLabels <- c()

  # identify unique dTs. This will later on be used to compute the corresponding discrete time values just once
  dTUnique <- unique(dT)
  uniqueNumbers <- seq_len(length(dTUnique))
  discreteNumbering <- rep(NA, length(dT))
  for(i in seq_len(length(dTUnique))){
    discreteNumbering[dT == dTUnique[i]] <- uniqueNumbers[i]
  }

  ctParameterNames <- names(ctMatrices)

  if("T0MEANS" %in% ctParameterNames | stationaryT0MEANS){
    MParameterIndicatorsLabels <- c(MParameterIndicatorsLabels, "T0MEANS")
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
    MParameterIndicatorsLabels <- c(MParameterIndicatorsLabels, discreteCINTNames)
  }

  if("MANIFESTMEANS" %in% ctParameterNames){
    MANIFESTMEANSNames <- rep("MANIFESTMEANS", Tpoints)
    MParameterIndicatorsLabels <- c(MParameterIndicatorsLabels, MANIFESTMEANSNames)
  }


  ## initialize M

  #MInitializer <- matrix(0, nrow = (nlatent)*Tpoints + nlatent + nmanifest*Tpoints, ncol = 1)
  MInitializer  <- matrix(mxObject$M$values, ncol = 1)
  rownames(MInitializer) <- names(mxObject$M$values)
  # create indicators for T0MEANS and CINT
  MParameterIndicators <- data.frame("label" = MParameterIndicatorsLabels,
                                     "row_start" = NA,
                                     "row_end" = NA,
                                     "col_start" = NA,
                                     "col_end" = NA
  )

  if("T0MEANS" %in% ctParameterNames | stationaryT0MEANS){

    # T0MEANS
    MParameterIndicators[1,"row_start"] <- 1
    MParameterIndicators[1,"row_end"] <- nlatent
    MParameterIndicators[1,"col_start"] <- 1
    MParameterIndicators[1,"col_end"] <- 1

    offset <- 1
  }else{
    offset <- 0
  }

  if("CINT" %in% ctParameterNames){

    # discreteCINT
    MParameterIndicators[1:length(discreteCINTNames)+offset,"row_start"] <- seq(nlatent+1, Tpoints*nlatent, nlatent)
    MParameterIndicators[1:length(discreteCINTNames)+offset,"row_end"] <- seq(nlatent+1, Tpoints*nlatent, nlatent)+nlatent-1
    MParameterIndicators[1:length(discreteCINTNames)+offset,"col_start"] <- 1
    MParameterIndicators[1:length(discreteCINTNames)+offset,"col_end"] <- 1

    retList[["discreteCINTUnique"]] <- discreteCINTUnique
    offset <- offset + length(discreteCINTNames)
  }

  if("MANIFESTMEANS" %in% ctParameterNames){
    if(!is.null(mxObject$T0TRAITEFFECT)){
      traitOffset <- nlatent
    }else{
      traitOffset <- 0
    }
    # MANIFESTMEANS
    MParameterIndicators[(offset+1) : (offset + length(MANIFESTMEANSNames) ),"row_start"] <- seq(from = Tpoints * nlatent + traitOffset + 1,
                                                                                                 to = Tpoints * nlatent + traitOffset + nmanifest*(Tpoints-1) + 1,
                                                                                                 by = nmanifest)
    MParameterIndicators[(offset+1) : (offset + length(MANIFESTMEANSNames) ),"row_end"] <- seq(from = Tpoints * nlatent + traitOffset + 1,
                                                                                               to = Tpoints * nlatent + traitOffset + nmanifest*(Tpoints-1) + 1,
                                                                                               by = nmanifest) + nmanifest - 1
    MParameterIndicators[(offset+1) : (offset + length(MANIFESTMEANSNames) ),"col_start"] <- 1
    MParameterIndicators[(offset+1) : (offset + length(MANIFESTMEANSNames) ),"col_end"] <- 1

  }

  # cpp indices start at 0
  cppMParameterIndicators <- MParameterIndicators
  cppMParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] <- cppMParameterIndicators[,c("row_start", "row_end", "col_start", "col_end")] - 1

  retList[["MParameterIndicators"]] = MParameterIndicators
  retList[["cppMParameterIndicators"]] = cppMParameterIndicators
  retList[["MInitializer"]] = MInitializer

  return(retList)

}


