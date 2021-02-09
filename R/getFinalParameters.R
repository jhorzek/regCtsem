#' getFinalParameters
#'
#' returns the final parameters for a
#' @param regCtsemObject fitted regularized continuous time model
#' @param criterion select a criterion. Possible are AIC, BIC, cvM2LL
#' @author Jannik Orzek
#' @import OpenMx
#' @export
getFinalParameters <- function(regCtsemObject, criterion = NULL){
  if(!regCtsemObject$setup$autoCV){
  minCriterionValue <- max(which(regCtsemObject$fit[criterion,] == min(regCtsemObject$fit[criterion,], na.rm = TRUE)))
  regValues <- regCtsemObject$setup$regValues
  bestRegValue <- regValues[minCriterionValue]
  return(list("criterion" = criterion,
              "lambda" = bestRegValue,
              "parameters" = regCtsemObject$parameters[,as.character(bestRegValue)]))
  }
  minCriterionValue <- max(which(regCtsemObject$fit["mean CV fit",] == min(regCtsemObject$fit["mean CV fit",], na.rm = TRUE)))
  regValues <- regCtsemObject$setup$regValues
  bestRegValue <- regValues[minCriterionValue]
  return(list("criterion" = "mean CV fit",
              "lambda" = bestRegValue))
}



