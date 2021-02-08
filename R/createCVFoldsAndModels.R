#' createCVFoldsAndModels
#'
#' @param mxObject Fitted object of class MxObject extracted from ctsemObject.
#' @param dataset only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param k number of cross-validation folds if autoCV = TRUE (k-fold cross-validation)
#' @author Jannik Orzek
#' @export
createCVFoldsAndModels <- function(mxObject, dataset = NULL, k){
  if(is.null(dataset)){
  fullData <- mxObject$data$observed}else{
    fullData <- dataset
  }

  cvFolds <- regCtsem::createFolds(mxObject = mxObject, dataset = dataset, k = k)

  testSets <- vector("list", length = k)
  names(testSets) <- names(cvFolds)

  trainSets <- vector("list", length = k)
  names(trainSets) <- names(cvFolds)

  cvModels <- vector("list", length = k)
  names(cvModels) <- names(cvFolds)

  return(list("fullData" = fullData, "cvFolds" = cvFolds, "testSets" = testSets, "trainSets" = trainSets, "cvModels" = cvModels))
}


