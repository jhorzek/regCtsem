#' createFolds
#'
#' created folds for automatic cross-validation
#'
#' @param mxObject Object of class mxModel
#' @param dataset only required if objective = "Kalman" and ctsemObject ist of type ctsemInit. Please provide a data set in wide format compatible to ctsemOMX
#' @param k number of cross-validation folds (k-fold cross-validation)
#' @author Jannik Orzek
#' @export
createFolds = function(mxObject, dataset = NULL, k){
  if(is.null(dataset)){
    sampleSize <- mxObject$data$numObs
    }else{
      sampleSize <- nrow(dataset)
    }

  # shuffle
  Folds <- sample(1:sampleSize, sampleSize)

  Folds <- split(Folds, cut(x = 1:length(Folds), breaks = k, labels = paste0("fold", 1:k)))

  return(Folds)
}



