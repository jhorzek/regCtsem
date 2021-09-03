#' summary.regCtsem
#'
#' @param regCtsemObject Object of type regCtsem
#' @param criterion select a criterion. Possible are AIC, BIC, cvM2LL
#' @author Jannik Orzek
#' @export
summary.regCtsem <- function(regCtsemObject, criterion = NULL){
  cat("\n")
  consoleWidth <- getOption("width")
  ret <- list()
  ret <- append(ret, paste0("|---", "regCtsem results", paste0(rep("-", consoleWidth - 30), collapse = ""), "|"))

  if(!is.null(criterion)){
    if(!criterion %in% rownames(regCtsemObject$fit)){
      removeEstimatedParameters <- !grepl("estimatedParameters", rownames(regCtsemObject$fit))
      stop(paste("Unknown criterion selected. Possible are: ", paste0(fitLabels[removeEstimatedParameters], collapse = ", ")))
    }
    if(any(is.na(regCtsemObject$fit[criterion,]))){
      warning("NAs in fit. Only using non-NA results")
    }
    minCriterionValue <- max(which(regCtsemObject$fit[criterion,] == min(regCtsemObject$fit[criterion,], na.rm = TRUE)))
    lambdas <- regCtsemObject$setup$lambdas
    bestLambda <- lambdas[minCriterionValue]
    ret <- append(ret, paste("The best ", criterion, " value was observed for:", sep = ""))
    ret <- append(ret, list("lambda" = bestLambda))
    ret <- append(ret, paste(" --- Parameter estimates: --- "))
    ret <- append(ret, list("ParameterEstimates" = regCtsemObject$parameters[,minCriterionValue]))
    ret <- append(ret, paste(" --- Fit: --- "))
    ret <- append(ret, list("fit" = regCtsemObject$fit[,minCriterionValue]))
    ret <- append(ret, paste0("|",paste0(rep("-", consoleWidth - 10), collapse = ""),"|"))
    return (ret)
  }else{

    ret <- append(ret, paste(" --- Parameter estimates: --- "))
    ret <- append(ret, list("ParameterEstimates" = regCtsemObject$parameters))
    ret <- append(ret, paste(" --- Fit: --- "))
    ret <- append(ret, list("fit" = regCtsemObject$fit))
    ret <- append(ret, paste0("|",paste0(rep("-", consoleWidth - 10), collapse = ""),"|"))
    return (ret)
  }
}

#' summary.regCtsemCV
#'
#' @param regCtsemObject Object of type regCtsemCV
#' @author Jannik Orzek
#' @export
summary.regCtsemCV <- function(regCtsemObject){
  cat("\n")
  consoleWidth <- getOption("width")
  ret <- list()
  ret <- append(ret, paste0("|---", "regCtsem results", paste0(rep("-", consoleWidth - 30), collapse = ""), "|"))

  fitLabels <- c(paste0("fold", 1:regCtsemObject$setup$k), "mean CV fit")

  if(any(is.na(regCtsemObject$fit))){
    warning("NAs in fit. Only using non-NA results")
  }

  minCriterionValue <- max(which(regCtsemObject$fit["mean",] == min(regCtsemObject$fit["mean",], na.rm = TRUE)))

  lambdas <- regCtsemObject$setup$lambdas
  bestLambda <- lambdas[minCriterionValue]
  ret <- append(ret, paste("The best CV fit value was observed for:", sep = ""))
  ret <- append(ret, list("lambda" = bestLambda))
  ret <- append(ret, paste("To compute the final parameter estimates, rerun the model with the full data set and lambda = ", bestLambda,".", sep = ""))
  ret <- append(ret, paste(" --- Fit: --- "))
  ret <- append(ret, list("fit" = regCtsemObject$fit))
  ret <- append(ret, paste0("|",paste0(rep("-", consoleWidth - 10), collapse = ""),"|"))
  return (ret)

}

#' summary.regCtsemMultiSubject
#'
#' @param regCtsemObject Object of type regCtsemMultiSubject
#' @param criterion select a criterion. Possible are AIC or BIC
#' @author Jannik Orzek
#' @export
summary.regCtsemMultiSubject <- function(regCtsemObject, criterion = NULL){
  cat("\n")
  consoleWidth <- getOption("width")
  ret <- list()
  ret <- append(ret, paste0("|---", "regCtsem results", paste0(rep("-", consoleWidth - 30), collapse = ""), "|"))

  if(!is.null(criterion)){
    if(!criterion %in% rownames(regCtsemObject$fit)){
      removeEstimatedParameters <- !grepl("estimatedParameters", rownames(regCtsemObject$fit))
      stop(paste("Unknown criterion selected. Possible are: ", paste0(fitLabels[removeEstimatedParameters], collapse = ", ")))
    }
    if(any(is.na(regCtsemObject$fit[criterion,]))){
      warning("NAs in fit. Only using non-NA results")
    }
    minCriterionValue <- max(which(regCtsemObject$fit[criterion,] == min(regCtsemObject$fit[criterion,], na.rm = TRUE)))
    lambdas <- regCtsemObject$setup$lambdas
    bestLambda <- lambdas[minCriterionValue]
    ret <- append(ret, paste("The best ", criterion, " value was observed for:", sep = ""))
    ret <- append(ret, list("lambda" = bestLambda))
    ret <- append(ret, paste(" --- Parameter estimates: --- "))
    param <- vector("list", length(regCtsemObject$parameters))
    names(param) <- names(regCtsemObject$parameters)
    for(i in 1:length(regCtsemObject$parameters)){
      param[[i]] <- regCtsemObject$parameters[[i]][,minCriterionValue]
    }
    ret <- append(ret, list("ParameterEstimates" = param))
    ret <- append(ret, paste(" Parameters were regularized towards the following values:"))
    ret <- append(ret, list("Targets" = regCtsemObject$setup$targetVector))
    ret <- append(ret, paste(" --- Fit: --- "))
    ret <- append(ret, list("fit" = regCtsemObject$fit[,minCriterionValue]))
    ret <- append(ret, paste0("|",paste0(rep("-", consoleWidth - 10), collapse = ""),"|"))
    return (ret)
  }else{

    ret <- append(ret, paste(" --- Parameter estimates: --- "))
    ret <- append(ret, list("ParameterEstimates" = regCtsemObject$parameters))
    ret <- append(ret, paste(" --- Fit: --- "))
    ret <- append(ret, list("fit" = regCtsemObject$fit))
    ret <- append(ret, paste0("|",paste0(rep("-", consoleWidth - 10), collapse = ""),"|"))
    return (ret)
  }

}





