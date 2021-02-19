prepareRAMData <- function(dataset, individualMissingPatternID, uniqueMissingPatterns){

  if(!is.matrix(uniqueMissingPatterns)){
    cols = length(uniqueMissingPatterns)
    uniqueMissingPatterns <- matrix(uniqueMissingPatterns, nrow = 1, ncol = cols)
  }

  uniqueMissingPatternIDs <- unique(individualMissingPatternID)
  patternNames <- paste0("missingnessPattern_", seq_len(length(uniqueMissingPatternIDs)))
  datasetRAM <- vector("list", length = length(uniqueMissingPatternIDs)+1)
  names(datasetRAM) <- c("names", patternNames)
  datasetRAM[["names"]] <- patternNames

  for(id in uniqueMissingPatternIDs){
    patternName <- patternNames[id]
    datasetRAM[[patternName]] <- list(
                             "sampleSize" = NA,
                             "nObservedVariables" = NA,
                             "rawData" = NA,
                             "observedMeans" = NA,
                             "observedCov" = NA,
                             "cppExpectedSelector" = NA,
                             "expectedMeans" = NA,
                             "expectedCovariance" = NA,
                             "m2LL" = -999999.9)
    datasetRAM[[patternName]] [["rawData"]] <- matrix(dataset[individualMissingPatternID == id,!uniqueMissingPatterns[id,]], nrow = sum(individualMissingPatternID == id))
    datasetRAM[[patternName]] [["sampleSize"]] <- nrow(datasetRAM[[patternName]] [["rawData"]])
    datasetRAM[[patternName]] [["nObservedVariables"]] <- sum(!uniqueMissingPatterns[id,])

    if(datasetRAM[[patternName]] [["sampleSize"]] > 1){
      datasetRAM[[patternName]] [["observedMeans"]] <- matrix(apply(datasetRAM[[patternName]] [["rawData"]], 2, mean), ncol = 1)
      datasetRAM[[patternName]] [["observedCov"]] <- ((datasetRAM[[patternName]] [["sampleSize"]]-1)/(datasetRAM[[patternName]] [["sampleSize"]]))*cov(datasetRAM[[patternName]] [["rawData"]])
    }else{
      datasetRAM[[patternName]] [["rawData"]] <- t(datasetRAM[[patternName]] [["rawData"]])
      datasetRAM[[patternName]] [["observedMeans"]] <- matrix(-999999.9, nrow = datasetRAM[[patternName]] [["nObservedVariables"]], ncol = 1)
      datasetRAM[[patternName]] [["observedCov"]] <- matrix(-999999.9, nrow = datasetRAM[[patternName]] [["nObservedVariables"]], ncol = datasetRAM[[patternName]] [["nObservedVariables"]])
    }

    datasetRAM[[patternName]] [["cppExpectedSelector"]] <- regCtsem::getRowsAndColsOfNonMissing(uniqueMissingPatterns[id,])

    datasetRAM[[patternName]] [["expectedMeans"]] <- matrix(-999999.9, nrow = datasetRAM[[patternName]] [["nObservedVariables"]], ncol = 1)
    datasetRAM[[patternName]] [["expectedCovariance"]] <- matrix(-999999.9, nrow = datasetRAM[[patternName]] [["nObservedVariables"]], ncol = datasetRAM[[patternName]] [["nObservedVariables"]])
  }

  return(datasetRAM)
}

