#' summary.regCtsem
#'
#' @param criterion select a criterion. Possible are AIC, BIC, cvM2LL
#' @param returnSummary Boolean: should the summary results be returned
#' @author Jannik Orzek
#' @import OpenMx
#' @export
summary.regCtsem <- function(regCtsemObject, criterion = NULL){
  cat("\n")
  consoleWidth <- getOption("width")
  sum1 <- paste0("|---", "regCtsem results", paste0(rep("-", consoleWidth - 16-5), collapse = ""), "|\n")
  # without automatic cross-validation
  if(!regCtsemObject$setup$autoCV){

    if(!is.null(criterion)){
      if(!criterion %in% rownames(regCtsemObject$fit)){
        removeEstimatedParameters <- !grepl("estimatedParameters", rownames(regCtsemObject$fit))
        stop(paste("Unknown criterion selected. Possible are: ", paste0(fitLabels[removeEstimatedParameters], collapse = ", ")))
      }
      if(any(is.na(regCtsemObject$fit[criterion,]))){
        warning("NAs in fit. Only using non-NA results")
      }
      minCriterionValue <- max(which(regCtsemObject$fit[criterion,] == min(regCtsemObject$fit[criterion,], na.rm = TRUE)))
      regValues <- as.numeric(colnames(regCtsemObject$fit))
      bestRegValue <- regValues[minCriterionValue]
      sum2 <- paste("\nThe best ", criterion, " value was observed for lambda = ", bestRegValue,".\n", sep = "")
      sum3 <- paste("\n --- Parameter estimates: --- \n")
      sum4 <- round(regCtsemObject$parameters[,as.character(bestRegValue)],3)
      sum5 <- paste("\n --- Fit: --- \n")
      sum6 <- round(regCtsemObject$fit[,as.character(bestRegValue)],3)
      sum7 <- paste0("|",paste0(rep("-", consoleWidth - 2), collapse = ""),"|\n")
      cat(paste0(sum1,sum2,sum3))
      print(sum4)
      cat(sum5)
      print(sum6)
      cat(sum7)
    }else{

      sum3 <- paste("\n --- Parameter estimates: --- \n")
      sum4 <- round(regCtsemObject$parameters,3)
      sum5 <- paste("\n --- Fit: --- \n")
      sum6 <- round(regCtsemObject$fit,3)
      sum7 <- paste0("|",paste0(rep("-", consoleWidth - 2), collapse = ""),"|\n")
      cat(paste0(sum1,sum3))
      print(sum4)
      cat(sum5)
      print(sum6)
      cat(sum7)
    }
  }else{

    # with automatic cross-validation

    fitLabels <- c(paste0("fold", 1:regCtsemObject$setup$k), "mean CV fit")

    if(any(is.na(regCtsemObject$fit))){
      warning("NAs in fit. Only using non-NA results")
    }

    minCriterionValue <- max(which(regCtsemObject$fit["mean CV fit",] == min(regCtsemObject$fit["mean CV fit",], na.rm = TRUE)))
    regValues <- as.numeric(colnames(regCtsemObject$fit))
    bestRegValue <- regValues[minCriterionValue]
    sum2 <- paste("\nThe best ", "mean CV fit", " value was observed for lambda = ", bestRegValue,".\n", sep = "")
    sum3 <- paste("To obtain the final parameter estimates, rerun the model with the full data set and regValue = ", bestRegValue,".\n")
    sum5 <- paste("\n --- Fit: --- \n")
    sum6 <- round(regCtsemObject$fit, 3)
    sum7 <- paste0("|",paste0(rep("-", consoleWidth - 2), collapse = ""),"|\n")
    cat(paste0(sum1,sum2))
    message(sum3)
    cat(sum5)
    print(sum6)
    cat(sum7)

  }
}


