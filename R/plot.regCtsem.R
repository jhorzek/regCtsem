#' plot.regCtsem
#'
#' @method plot regCtsem
#' @param regCtsemObject fitted regCtsem object
#' @param what what should be plotted? Possbile are: 'drift', 'parameters', 'fit'
#' @param criterion vector with labels of criteria which should be plotted when what = 'fit'
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @param ... additional parameters for plot or matplot (e.g., lty, col, ...)
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
plot.regCtsem <- function(regCtsemObject, what = "drift", criterion = "BIC", xlab = "auto", ylab = "auto",...){
  if(is.character(xlab)){skiptXlabComp <- F}else{skiptXlabComp <- T}
  if(is.character(ylab)){skiptYlabComp <- F}else{skiptYlabComp <- T}
  plotted <- FALSE
  if(tolower(what) == "drift"){
    if(!"parameters" %in% names(regCtsemObject)){stop("Plot of drift values not possible for cross-validation. Use what = 'fit' to plot the fit.")}

    driftIndicator <- grepl("drift", rownames(regCtsemObject$parameters))
    regValues <- regCtsemObject$setup$regValues

    if(any(is.na(regCtsemObject$parameters))){
      warning("NAs in regCtsemObject$parameters. Only plotting non-NA values")
    }

    colsWithNA <- apply(regCtsemObject$parameters,2,function(x) any(is.na(x)))

    matplot(x = regValues[!colsWithNA], y = t(regCtsemObject$parameters[driftIndicator,!colsWithNA]),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","drift values",ylab)),
            xlab = ifelse(skiptXlabComp, xlab,
                          ifelse(xlab == "auto","regValue",xlab)),
            type = "l", ...)
    plotted <- TRUE
  }
  if(tolower(what) == "parameters"){
    if(!"parameters" %in% names(regCtsemObject)){stop("Plot of drift values not possible for cross-validation. Use what = 'fit' to plot the fit.")}

    if(any(is.na(regCtsemObject$parameters))){
      warning("NAs in regCtsemObject$parameters. Only plotting non-NA values")
    }

    colsWithNA <- apply(regCtsemObject$parameters,2,function(x) any(is.na(x)))

    regValues <- regCtsemObject$setup$regValues
    matplot(x = regValues[!colsWithNA], y = t(regCtsemObject$parameters[,!colsWithNA]),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","parameter values",ylab)),
            xlab = ifelse(skiptXlabComp, xlab,
                          ifelse(xlab == "auto","regValue",xlab)),
            type = "l", ...)
    plotted <- TRUE
  }
  if(tolower(what) == "fit"){
    # cross-validation:
    if((!"fitAndParameters" %in% names(regCtsemObject)) & ("fit" %in% names(regCtsemObject))){
      regValues <- regCtsemObject$setup$regValues
      plot(x = regValues, y = regCtsemObject$fit["mean CV fit",],
           ylab = ifelse(skiptYlabComp, ylab,
                         ifelse(ylab == "auto","mean CV fit",ylab)),
           xlab = ifelse(skiptXlabComp, xlab,
                         ifelse(xlab == "auto","regValue",xlab)),
           type = "l", ...)
    }else if("fit" %in% names(regCtsemObject)){
      if(!criterion %in% rownames(fit)){
        stop(paste0(criterion, " not in fit values. Possible are: ", paste0(rownames(fit), collapse = ", "), "."))
      }

      if(any(is.na(regCtsemObject$fit))){
        warning("NAs in regCtsemObject$fit Only plotting non-NA values")
      }

      colsWithNA <- apply(regCtsemObject$fit,2,function(x) any(is.na(x)))

      regValues <- regCtsemObject$setup$regValues

      matplot(x = regValues[!colsWithNA], y = t(matrix(regCtsemObject$fit[criterion,!colsWithNA],
                                                nrow = length(criterion))),
              ylab = ifelse(skiptYlabComp, ylab,
                            ifelse(ylab == "auto","fit values",ylab)),
              xlab = ifelse(skiptXlabComp, xlab,
                            ifelse(xlab == "auto","regValue",xlab)),
              lty = 1:length(criterion),
              type = "l", ...)
    }
    legend(legend = criterion, "bottomright", lty = 1:length(criterion))
    plotted <- TRUE
  }

  if(!plotted){
    stop("Select what = 'drift', 'parameters', or 'fit' to generate a plot.")
  }

}


