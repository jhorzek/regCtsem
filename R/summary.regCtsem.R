#' summary.regCtsem
#'
#' @param object Object of type regCtsem
#' @param ... NULL
#' @param criterion select a criterion. Possible are AIC, BIC, cvM2LL
#' @author Jannik Orzek
#' @export
summary.regCtsem <- function(object, ..., criterion = NULL){
  cat("\n")
  consoleWidth <- getOption("width")
  ret <- list()
  ret <- append(ret, paste0("|---", "regCtsem results", paste0(rep("-", consoleWidth - 30), collapse = ""), "|"))

  if(!is.null(criterion)){
    if(!criterion %in% rownames(object$fit)){
      removeEstimatedParameters <- !grepl("estimatedParameters", rownames(object$fit))
      stop(paste("Unknown criterion selected. Possible are: ", paste0(fitLabels[removeEstimatedParameters], collapse = ", ")))
    }
    if(any(is.na(object$fit[criterion,]))){
      warning("NAs in fit. Only using non-NA results")
    }
    minCriterionValue <- max(which(object$fit[criterion,] == min(object$fit[criterion,], na.rm = TRUE)))
    lambdas <- object$setup$lambdas
    bestLambda <- lambdas[minCriterionValue]
    ret <- append(ret, paste("The best ", criterion, " value was observed for:", sep = ""))
    ret <- append(ret, list("lambda" = bestLambda))
    ret <- append(ret, paste(" --- Parameter estimates: --- "))
    ret <- append(ret, list("ParameterEstimates" = object$parameters[,minCriterionValue]))
    ret <- append(ret, paste(" --- Fit: --- "))
    ret <- append(ret, list("fit" = object$fit[,minCriterionValue]))
    ret <- append(ret, paste0("|",paste0(rep("-", consoleWidth - 10), collapse = ""),"|"))
    return (ret)
  }else{

    ret <- append(ret, paste(" --- Parameter estimates: --- "))
    ret <- append(ret, list("ParameterEstimates" = object$parameters))
    ret <- append(ret, paste(" --- Fit: --- "))
    ret <- append(ret, list("fit" = object$fit))
    ret <- append(ret, paste0("|",paste0(rep("-", consoleWidth - 10), collapse = ""),"|"))
    return (ret)
  }
}

#' summary.regCtsemCV
#'
#' @param object Object of type regCtsemCV
#' @param ... NULL
#' @author Jannik Orzek
#' @export
summary.regCtsemCV <- function(object, ...){
  cat("\n")
  consoleWidth <- getOption("width")
  ret <- list()
  ret <- append(ret, paste0("|---", "regCtsem results", paste0(rep("-", consoleWidth - 30), collapse = ""), "|"))

  fitLabels <- c(paste0("fold", 1:object$setup$k), "mean CV fit")

  if(any(is.na(object$fit))){
    warning("NAs in fit. Only using non-NA results")
  }

  minCriterionValue <- max(which(object$fit["mean",] == min(object$fit["mean",], na.rm = TRUE)))

  lambdas <- object$setup$lambdas
  bestLambda <- lambdas[minCriterionValue]
  ret <- append(ret, paste("The best CV fit value was observed for:", sep = ""))
  ret <- append(ret, list("lambda" = bestLambda))
  ret <- append(ret, paste("To compute the final parameter estimates, rerun the model with the full data set and lambda = ", bestLambda,".", sep = ""))
  ret <- append(ret, paste(" --- Fit: --- "))
  ret <- append(ret, list("fit" = object$fit))
  ret <- append(ret, paste0("|",paste0(rep("-", consoleWidth - 10), collapse = ""),"|"))
  return (ret)

}

#' summary.regCtsemMultiSubject
#'
#' @param object Object of type regCtsemMultiSubject
#' @param ... NULL
#' @param criterion select a criterion. Possible are AIC or BIC
#' @author Jannik Orzek
#' @export
summary.regCtsemMultiSubject <- function(object, ..., criterion = NULL){
  cat("\n")
  consoleWidth <- getOption("width")
  ret <- list()
  ret <- append(ret, paste0("|---", "regCtsem results", paste0(rep("-", consoleWidth - 30), collapse = ""), "|"))

  if(!is.null(criterion)){
    if(!criterion %in% rownames(object$fit)){
      removeEstimatedParameters <- !grepl("estimatedParameters", rownames(object$fit))
      stop(paste("Unknown criterion selected. Possible are: ", paste0(fitLabels[removeEstimatedParameters], collapse = ", ")))
    }
    if(any(is.na(object$fit[criterion,]))){
      warning("NAs in fit. Only using non-NA results")
    }
    minCriterionValue <- max(which(object$fit[criterion,] == min(object$fit[criterion,], na.rm = TRUE)))
    lambdas <- object$setup$lambdas
    bestLambda <- lambdas[minCriterionValue]
    ret <- append(ret, paste("The best ", criterion, " value was observed for:", sep = ""))
    ret <- append(ret, list("lambda" = bestLambda))
    ret <- append(ret, paste(" --- Parameter estimates: --- "))
    param <- vector("list", length(object$parameters))
    names(param) <- names(object$parameters)
    for(i in 1:length(object$parameters)){
      param[[i]] <- object$parameters[[i]][,minCriterionValue]
    }
    ret <- append(ret, list("ParameterEstimates" = param))
    ret <- append(ret, paste(" Parameters were regularized towards the following values:"))
    ret <- append(ret, list("Targets" = object$setup$targetVector))
    ret <- append(ret, paste(" --- Fit: --- "))
    ret <- append(ret, list("fit" = object$fit[,minCriterionValue]))
    ret <- append(ret, paste0("|",paste0(rep("-", consoleWidth - 10), collapse = ""),"|"))
    return (ret)
  }else{

    ret <- append(ret, paste(" --- Parameter estimates: --- "))
    ret <- append(ret, list("ParameterEstimates" = object$parameters))
    ret <- append(ret, paste(" --- Fit: --- "))
    ret <- append(ret, list("fit" = object$fit))
    ret <- append(ret, paste0("|",paste0(rep("-", consoleWidth - 10), collapse = ""),"|"))
    return (ret)
  }

}





