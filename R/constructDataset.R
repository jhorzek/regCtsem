#' constructDataset
#'
#' Separates data and time intervals. Computes the number of unique time intervals and the unqiue missingness patterns
#' @param wideData dataset in wide format compatible with ctsem
#' @import mgcv
constructDataset <- function(wideData){
  # separate data from time intervals
  dT <- subset(wideData, select = grepl("dT", colnames(wideData)))

  if(any(apply(dT, 2, function(x) length(unique(x)))>1)){
    stop("Observed columns with more than one unqiue dT. This is not possible")
  }
  dT <- apply(dT, 2, unique)

  dataset <- subset(wideData, select = !grepl("dT", colnames(wideData)) & !grepl("intervalID", colnames(wideData)) )

  # extract unique dts
  dTUnique <- unique(dT)

  # extract unique missingness patterns

  isMissing <- is.na(dataset)
  uniqueMissingPatterns <- mgcv::uniquecombs(isMissing, ordered=FALSE)
  individualMissingPatternID <- attr(uniqueMissingPatterns,"index")

  dataList <- list("dataset" = as.matrix(dataset),
                   "dT" = dT,
                   "dTUnique" = dTUnique,
                   "uniqueMissingPatterns" = uniqueMissingPatterns,
                   "individualMissingPatternID" = individualMissingPatternID
                   )

  return(dataList)
}

