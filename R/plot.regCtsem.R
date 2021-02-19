#' plot.regCtsem
#'
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

    drifts <- regCtsemObject$parameters[grepl("drift", rownames(regCtsemObject$parameters)),]
    regIndicators <- regCtsemObject$setup$regIndicators
    drifts_regularized <- drifts[regIndicators,]
    regValues <- regCtsemObject$setup$regValues

    colsWithNA <- apply(regCtsemObject$parameters,2,function(x) any(is.na(x)))

    color <- ifelse(rownames(drifts) %in% regIndicators, yes = "white", "black")

    par(mar=c(4, 3, 5, 2), xpd=TRUE)

    matplot(x = regValues[!colsWithNA], t(drifts[,!colsWithNA]), lty = 2,
            lwd = 2, type = "l",
            col = color, xlab = ifelse(skiptXlabComp, xlab,
                                       ifelse(xlab == "auto",expression(lambda),ylab)),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","drift",ylab)))
    matplot(x = regValues[!colsWithNA], t(drifts_regularized[,!colsWithNA]), lty = 1, lwd = 2, type = "l", col = "#2166AC", add = TRUE)
    tickat <- seq(1, length(regValues[!colsWithNA]), length.out = 10)
    axis(3, at = regValues[tickat],
         labels=apply(drifts == 0,2,sum)[tickat],
         outer= F,
         line=1,col="black",col.ticks="black",col.axis="black")
    mtext("# zeroed parameters",3,line=3,at=mean(regValues),col="black", cex = 1)

    par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
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
    color <- c("#008000", "#008080", "#800080", "#800000", rep("black", 5))
    lty <- c(1,2,3,4, rep(1,5))
    # cross-validation:
    if((!"fitAndParameters" %in% names(regCtsemObject)) & ("fit" %in% names(regCtsemObject))){
      regValues <- regCtsemObject$setup$regValues
      plot(x = regValues, y = regCtsemObject$fit["mean CV fit",],
           ylab = ifelse(skiptYlabComp, ylab,
                         ifelse(ylab == "auto","mean CV fit",ylab)),
           xlab = ifelse(skiptXlabComp, xlab,
                         ifelse(xlab == "auto",expression(lambda),xlab)),
           type = "l", lwd = 2, col = color[2], ...)
    }else if("fit" %in% names(regCtsemObject)){
      if(any(!criterion %in% rownames(regCtsemObject$fit))){
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
                            ifelse(xlab == "auto",expression(lambda),xlab)),
              lty = lty[1:length(criterion)], col = color[1:length(criterion)],
              type = "l", lwd = 2, ...)
    }
    legend(legend = criterion, "bottomright",
           lty = lty[1:length(criterion)],
           col = color[1:length(criterion)])
    plotted <- TRUE
  }

  if(!plotted){
    stop("Select what = 'drift', 'parameters', or 'fit' to generate a plot.")
  }

}

