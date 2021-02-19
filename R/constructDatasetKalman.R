#' constructDatasetKalman
#'
#' Separates data and time intervals. Computes the number of unique time intervals and the unqiue missingness patterns
#' @param wideData dataset in wide format compatible with ctsem
#' @import mgcv
constructDatasetKalman <- function(wideData){
  wideData <- as.matrix(wideData)
  # separate data from time intervals
  dT <- subset(wideData, select = grepl("dT", colnames(wideData)))

  dataset <- subset(wideData, select = !grepl("dT", colnames(wideData)) & !grepl("intervalID", colnames(wideData)) )

  # extract unique dts
  dTUnique <- unique(as.vector(dT))

  # make dt indicators
  dtIndicators <- matrix(NA, nrow = nrow(dT), ncol = ncol(dT))
  for(i in seq_len(length(dTUnique))){
    dtIndicators[dT == dTUnique[i]] <- i-1
  }

  # extract unique missingness patterns

  isMissing <- is.na(dataset)
  uniqueMissingPatterns <- mgcv::uniquecombs(isMissing, ordered=FALSE)
  individualMissingPatternID <- attr(uniqueMissingPatterns,"index")

  dataList <- list("dataset" = as.matrix(dataset),
                   "dT" = dT,
                   "dTUnique" = dTUnique,
                   "dtIndicators" = dtIndicators,
                   "uniqueMissingPatterns" = uniqueMissingPatterns,
                   "individualMissingPatternID" = individualMissingPatternID
  )

  return(dataList)
}

