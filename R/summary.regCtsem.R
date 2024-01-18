#' print summary of regCtsem
#' @param x object of class regCtsem
#' @param ... not used
#' @return nothing
#' @method print regCtsemSummary
#' @export
print.regCtsemSummary <- function(x, ...){

  if(!is.null(x$criterion)){

    cat("\n")
    cat("#### Results of regCtsem ####\n\n")
    cat(paste0("Criterion for model selection: ", x$criterion, "\n"))
    cat(paste0("Criterion value:               ", round(x$fit[x$criterion],3), "\n"))
    cat(paste0("Lambda value:                  ", round(x$bestLambda,3), "\n\n"))
    cat(paste0("## Estimates ##\n"))

    params <- x$parameterEstimates
    paramsPrint <- round(params, 3)
    paramsPrint[params == 0] <- "."

    print(paramsPrint,
          quote = FALSE,
          right = TRUE)
  }else{

    cat("\n")
    cat("#### Results from regCtsem ####\n\n")

    fits <- x$fit

    cat("## Fit Measures ##\n")

    print(fits,
          quote = FALSE,
          right = TRUE)

    cat(paste0("\n\n## Estimates ##\n"))

    params <- x$parameterEstimates
    paramsPrint <- round(params, 3)
    paramsPrint[params == 0] <- "."

    print(paramsPrint,
          quote = FALSE,
          right = TRUE)

  }
  cat("\n")
}

#' print summary of regCtsemCV
#' @param x object of class regCtsemCV
#' @param ... not used
#' @return nothing
#' @method print regCtsemCVSummary
#' @export
print.regCtsemCVSummary <- function(x, ...){

  cat("\n")

  cat("#### Results of cross-validated regCtsem ####\n\n")

  cat(paste0("Criterion for model selection: cross-validation fit\n"))
  cat(paste0("Lowest mean cv fit:            ", round(x$bestFit, 3), "\n"))
  cat(paste0("Lambda value:                  ", round(x$bestLambda,3), "\n\n"))

  cat(
    paste("To compute the final parameter estimates, rerun the model with the full data set and lambda = ",
          x$bestLambda,".\n", sep = ""))

  cat("\n")
}


#' summary.regCtsem
#'
#' @param object Object of type regCtsem
#' @param ... NULL
#' @param criterion select a criterion. Possible are AIC, BIC, cvM2LL
#' @author Jannik Orzek
#' @export
summary.regCtsem <- function(object, ..., criterion = NULL){

  ret <- list()

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
    ret$criterion = criterion
    ret$bestLambda = bestLambda
    ret$parameterEstimates = object$parameters[,minCriterionValue]
    ret$fit = object$fit[,minCriterionValue]

  }else{

    ret$parameterEstimates <- object$parameters
    ret$fit <- object$fit

  }

  class(ret) <- "regCtsemSummary"

  return(ret)

}

#' summary.regCtsemCV
#'
#' @param object Object of type regCtsemCV
#' @param ... NULL
#' @author Jannik Orzek
#' @export
summary.regCtsemCV <- function(object, ...){

  ret <- list()

  if(any(is.na(object$fit))){
    warning("NAs in fit. Only using non-NA results")
  }

  minCriterionValue <- max(which(object$fit["mean",] == min(object$fit["mean",], na.rm = TRUE)))


  ret$bestFit <- object$fit["mean", minCriterionValue]
  lambdas <- object$setup$lambdas
  ret$bestLambda <- lambdas[minCriterionValue]
  ret$fit <- object$fit

  class(ret) <- "regCtsemCVSummary"

  return(ret)

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
