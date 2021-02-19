prepareDiscreteElementNames <- function(ctMatrices, Tpoints, dT){

  retList <- list()

  # identify unique dTs. This will later on be used to compute the corresponding discrete time values just once
  dTUnique <- unique(as.vector(dT))
  uniqueNumbers <- seq_len(length(dTUnique))
  discreteNumbering <- matrix(NA, nrow = nrow(dT), ncol = ncol(dT))
  for(i in seq_len(length(dTUnique))){
    discreteNumbering[dT == dTUnique[i]] <- uniqueNumbers[i]
  }

  ctParameterNames <- names(ctMatrices)

  if("DRIFT" %in% ctParameterNames){
    discreteDRIFTNames <- matrix(paste0("discreteDRIFT_", discreteNumbering),
                                 nrow = nrow(discreteNumbering),
                                 ncol = ncol(discreteNumbering))
    retList[["discreteDRIFTNames"]] <- discreteDRIFTNames
  }

  if("T0TRAITEFFECT" %in% ctParameterNames){
    discreteTRAITNames <- matrix(paste0("discreteTRAIT_", discreteNumbering),
                                 nrow = nrow(discreteNumbering),
                                 ncol = ncol(discreteNumbering))
    retList[["discreteTRAITNames"]] <- discreteTRAITNames
  }

  if("DIFFUSIONbase" %in% ctParameterNames){
    discreteDIFFUSIONNames <- matrix(paste0("discreteDIFFUSION_", discreteNumbering),
                                     nrow = nrow(discreteNumbering),
                                     ncol = ncol(discreteNumbering))
    DRIFTHASHExponentialNames <- matrix(paste0("DRIFTHASHExponential_", discreteNumbering),
                                        nrow = nrow(discreteNumbering),
                                        ncol = ncol(discreteNumbering))

    retList[["discreteDIFFUSIONNames"]] = discreteDIFFUSIONNames
    retList[["DRIFTHASHExponentialNames"]] = DRIFTHASHExponentialNames
  }

  if("CINT" %in% ctParameterNames){
    discreteCINTNames <- matrix(paste0("discreteCINT_", discreteNumbering),
                                nrow = nrow(discreteNumbering),
                                ncol = ncol(discreteNumbering))

    retList[["discreteCINTNames"]] <- discreteCINTNames
  }


  return(retList)
}

