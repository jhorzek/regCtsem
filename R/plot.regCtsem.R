#' plot.regCtsem
#'
#' @param regCtsemObject fitted regCtsem object
#' @param what what should be plotted? Possbile are: 'drift', 'parameters', 'fit' or labels of specific parameters
#' @param criterion vector with labels of criteria which should be plotted when what = 'fit'
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @param ... additional parameters for plot or matplot (e.g., lty, col, ...)
#' @author Jannik Orzek
#' @export
plot.regCtsem <- function(regCtsemObject, what = "drift", criterion = "BIC", xlab = "auto", ylab = "auto",...){
  if(is.character(xlab)){skiptXlabComp <- F}else{skiptXlabComp <- T}
  if(is.character(ylab)){skiptYlabComp <- F}else{skiptYlabComp <- T}

  ## Plot specific parameter values
  if(all(what %in% rownames(regCtsemObject$parameters))){
    pars <- regCtsemObject$parameters
    parsSelected <- subset(pars, subset = rownames(pars) %in% what)
    regIndicators <- regCtsemObject$setup$regIndicators
    targetVector <- regCtsemObject$setup$targetVector
    pars_regularized <- subset(parsSelected, subset = rownames(parsSelected) %in% regIndicators)

    lambdas <- regCtsemObject$setup$lambdas

    colsWithNA <- apply(pars,2,function(x) any(is.na(x)))

    color <- ifelse(rownames(pars) %in% regIndicators, yes = "white", "black")

    par(mar=c(4, 3, 5, 2), xpd=TRUE)

    matplot(x = lambdas[!colsWithNA], t(subset(parsSelected, select = !colsWithNA)), lty = 2,
            lwd = 2, type = "l",
            col = color, xlab = ifelse(skiptXlabComp, xlab,
                                       ifelse(xlab == "auto",expression(lambda),ylab)),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto", "Values",ylab)))
    matplot(x = lambdas[!colsWithNA], t(subset(pars_regularized, select = !colsWithNA)), lty = 1, lwd = 2, type = "l", col = "#2166AC", add = TRUE)
    tickat <- seq(1, length(lambdas[!colsWithNA]), length.out = 10)
    axis(3, at = lambdas[tickat],
         labels=apply(pars[regIndicators,] == targetVector[regIndicators],2,sum)[tickat],
         outer= F,
         line=1,col="black",col.ticks="black",col.axis="black")
    mtext("# parameters on target",3,line=3,at=mean(lambdas),col="black", cex = 1)

    par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
    return()
  }

  if(tolower(what) == "drift"){
    pars <- regCtsemObject$parameters
    regIndicators <- regCtsemObject$setup$regIndicators
    targetVector <- regCtsemObject$setup$targetVector
    drifts <- subset(pars, subset = grepl("drift", rownames(pars)))
    drifts_regularized <- subset(drifts, rownames(drifts) %in% regIndicators)
    lambdas <- regCtsemObject$setup$lambdas

    colsWithNA <- apply(pars,2,function(x) any(is.na(x)))

    color <- ifelse(rownames(drifts) %in% regIndicators, yes = "white", "black")

    par(mar=c(4, 3, 5, 2), xpd=TRUE)

    matplot(x = lambdas[!colsWithNA], t(subset(drifts, select = !colsWithNA)), lty = 2,
            lwd = 2, type = "l",
            col = color, xlab = ifelse(skiptXlabComp, xlab,
                                       ifelse(xlab == "auto",expression(lambda),ylab)),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","drift",ylab)))
    matplot(x = lambdas[!colsWithNA], t(subset(drifts_regularized, select = !colsWithNA)), lty = 1, lwd = 2, type = "l", col = "#2166AC", add = TRUE)
    tickat <- seq(1, length(lambdas[!colsWithNA]), length.out = 10)
    axis(3, at = lambdas[tickat],
         labels=apply(pars[regIndicators,] == targetVector[regIndicators],2,sum)[tickat],
         outer= F,
         line=1,col="black",col.ticks="black",col.axis="black")
    mtext("# parameters on target",3,line=3,at=mean(lambdas),col="black", cex = 1)

    par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
    return()
  }
  if(tolower(what) == "parameters"){
    if(!"parameters" %in% names(regCtsemObject)){stop("Plot of drift values not possible for cross-validation. Use what = 'fit' to plot the fit.")}

    if(any(is.na(regCtsemObject$parameters))){
      warning("NAs in regCtsemObject$parameters. Only plotting non-NA values")
    }

    colsWithNA <- apply(regCtsemObject$parameters,2,function(x) any(is.na(x)))

    lambdas <- regCtsemObject$setup$lambdas

    matplot(x = lambdas[!colsWithNA], y = t(subset(regCtsemObject$parameters, select = !colsWithNA)),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","parameter values",ylab)),
            xlab = ifelse(skiptXlabComp, xlab,
                          ifelse(xlab == "auto","lambda",xlab)),
            type = "l", ...)
    plotted <- TRUE
  }
  if(tolower(what) == "fit"){
    if(! "fit" %in% names(regCtsemObject)){
      stop("Could not find a field 'fit' in the regCtsem object'")
    }
    if(any(!criterion %in% rownames(regCtsemObject$fit))){
      stop(paste0(criterion, " not in fit values. Possible are: ", paste0(rownames(fit), collapse = ", "), "."))
    }

    color <- c("#008000", "#008080", "#800080", "#800000", rep("black", 5))
    lty <- c(1,2,3,4, rep(1,5))

    if(any(is.na(regCtsemObject$fit))){
      warning("NAs in regCtsemObject$fit Only plotting non-NA values")
    }

    colsWithNA <- apply(regCtsemObject$fit,2,function(x) any(is.na(x)))

    lambdas <- regCtsemObject$setup$lambdas

    matplot(x = lambdas[!colsWithNA], y = t(matrix(regCtsemObject$fit[criterion,!colsWithNA],
                                                   nrow = length(criterion))),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","fit values",ylab)),
            xlab = ifelse(skiptXlabComp, xlab,
                          ifelse(xlab == "auto",expression(lambda),xlab)),
            lty = lty[1:length(criterion)], col = color[1:length(criterion)],
            type = "l", lwd = 2, ...)

    legend(legend = criterion, "bottomright",
           lty = lty[1:length(criterion)],
           col = color[1:length(criterion)])
    return()
  }

  stop("Select what = 'drift', 'parameters', or 'fit' to generate a plot.")

}

#' plot.regCtsemCV
#'
#' @param regCtsemObject fitted regCtsem object
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @param ... additional parameters for plot or matplot (e.g., lty, col, ...)
#' @author Jannik Orzek
#' @export
plot.regCtsemCV <- function(regCtsemObject, xlab = "auto", ylab = "auto",...){
  if(is.character(xlab)){skiptXlabComp <- F}else{skiptXlabComp <- T}
  if(is.character(ylab)){skiptYlabComp <- F}else{skiptYlabComp <- T}

  criterion <- "mean"
  lambdas <- regCtsemObject$setup$lambdas

  plot(x = lambdas, y = regCtsemObject$fit[criterion,],
       ylab = ifelse(skiptYlabComp, ylab,
                     ifelse(ylab == "auto","mean CV fit",ylab)),
       xlab = ifelse(skiptXlabComp, xlab,
                     ifelse(xlab == "auto",expression(lambda),xlab)),
       type = "l", lwd = 2, col = "#008080", ...)

}

#' plot.regCtsemMultiSubject
#'
#' @param regCtsemObject fitted regCtsem object
#' @param what what should be plotted? Possbile are: 'drift', 'parameters', 'fit'
#' @param groups For which groups should the parameters be plotted (numeric vector)?
#' @param criterion vector with labels of criteria which should be plotted when what = 'fit'
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @param ... additional parameters for plot or matplot (e.g., lty, col, ...)
#' @author Jannik Orzek
#' @export
plot.regCtsemMultiSubject <- function(regCtsemObject, what = "drift", groups = NULL, criterion = "BIC", xlab = "auto", ylab = "auto",...){
  if(is.character(xlab)){skiptXlabComp <- F}else{skiptXlabComp <- T}
  if(is.character(ylab)){skiptYlabComp <- F}else{skiptYlabComp <- T}
  if(is.null(groups)){
    groups <- 1:length(regCtsemObject$parameters)
  }

  ## Plot specific parameter values
  # if(any(c(what, paste0(what, "_G1")) %in% rownames(regCtsemObject$parameters[[1]]))){
  #   for(group in groups){
  #     pars <- regCtsemObject$parameters[[group]]
  #     parsSelected <- subset(pars, subset = rownames(pars) %in% what)
  #     targetVector <- regCtsemObject$setup$targetVector
  #     regIndicators <- regCtsemObject$setup$regIndicators
  #     pars_regularized <- subset(pars, subset = rownames(pars) %in% regIndicators)
  #
  #     lambdas <- regCtsemObject$setup$lambdas
  #
  #     colsWithNA <- apply(pars,2,function(x) any(is.na(x)))
  #
  #     color <- ifelse(rownames(pars) %in% regIndicators, yes = "white", "black")
  #
  #     par(mar=c(4, 3, 5, 2), xpd=TRUE)
  #
  #     matplot(x = lambdas[!colsWithNA], t(pars[,!colsWithNA]), lty = 2,
  #             lwd = 2, type = "l",
  #             col = color, xlab = ifelse(skiptXlabComp, xlab,
  #                                        ifelse(xlab == "auto",expression(lambda),ylab)),
  #             ylab = ifelse(skiptYlabComp, ylab,
  #                           ifelse(ylab == "auto", "Values",ylab)))
  #     matplot(x = lambdas[!colsWithNA], t(subset(pars_regularized, select = !colsWithNA)), lty = 1, lwd = 2, type = "l", col = "#2166AC", add = TRUE)
  #     tickat <- seq(1, length(lambdas[!colsWithNA]), length.out = 10)
  #     axis(3, at = lambdas[tickat],
  #          labels=apply(pars[regIndicators[regIndicators %in% rownames(pars)],] == targetVector[regIndicators[regIndicators %in% rownames(pars)]],2,sum)[tickat],
  #          outer= F,
  #          line=1,col="black",col.ticks="black",col.axis="black")
  #     mtext(paste0(paste0("Group ", group), ": # parameters on target"),3,line=3,at=mean(lambdas),col="black", cex = 1)
  #
  #     readline(prompt="Press [Enter] for next plot... ")
  #   }
  #   par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
  #   return()
  # }

  if(tolower(what) == "drift"){
    for(group in groups){
      pars <- regCtsemObject$parameters[[group]]
      drifts <- subset(pars, subset = grepl("drift", rownames(pars)))
      regIndicators <- regCtsemObject$setup$regIndicators
      targetVector <- regCtsemObject$setup$targetVector
      drifts_regularized <- subset(drifts, rownames(drifts) %in% regIndicators)
      lambdas <- regCtsemObject$setup$lambdas

      colsWithNA <- apply(pars,2,function(x) any(is.na(x)))

      color <- ifelse(rownames(drifts) %in% regIndicators, yes = "white", "black")

      par(mar=c(4, 3, 5, 2), xpd=TRUE)

      matplot(x = lambdas[!colsWithNA], t(subset(drifts, select = !colsWithNA)), lty = 2,
              lwd = 2, type = "l",
              col = color, xlab = ifelse(skiptXlabComp, xlab,
                                         ifelse(xlab == "auto",expression(lambda),ylab)),
              ylab = ifelse(skiptYlabComp, ylab,
                            ifelse(ylab == "auto","drift",ylab)))
      matplot(x = lambdas[!colsWithNA], t(subset(drifts_regularized, select = !colsWithNA)), lty = 1, lwd = 2, type = "l", col = "#2166AC", add = TRUE)
      tickat <- seq(1, length(lambdas[!colsWithNA]), length.out = 10)
      axis(3, at = lambdas[tickat],
           labels=apply(pars[regIndicators[regIndicators %in% rownames(pars)],] == targetVector[regIndicators[regIndicators %in% rownames(pars)]],2,sum)[tickat],
           outer= F,
           line=1,col="black",col.ticks="black",col.axis="black")
      mtext(paste0(paste0("Group ", group), ": # parameters on target"),3,line=3,at=mean(lambdas),col="black", cex = 1)

      readline(prompt="Press [Enter] for next plot... ")
    }
    par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
    return()
  }
  if(tolower(what) == "parameters"){
    if(!"parameters" %in% names(regCtsemObject)){stop("Plot of drift values not possible for cross-validation. Use what = 'fit' to plot the fit.")}

    if(any(is.na(regCtsemObject$parameters[[group]]))){
      warning("NAs in regCtsemObject$parameters. Only plotting non-NA values")
    }
    for(group in groups){
      pars <- regCtsemObject$parameters[[group]]
      regIndicators <- regCtsemObject$setup$regIndicators
      targetVector <- regCtsemObject$setup$targetVector

      par(mar=c(4, 3, 5, 2), xpd=TRUE)
      colsWithNA <- apply(regCtsemObject$parameters[[group]],2,function(x) any(is.na(x)))

      lambdas <- regCtsemObject$setup$lambdas

      matplot(x = lambdas[!colsWithNA], y = t(subset(regCtsemObject$parameters[[group]], select = !colsWithNA)),
              ylab = ifelse(skiptYlabComp, ylab,
                            ifelse(ylab == "auto","parameter values",ylab)),
              xlab = ifelse(skiptXlabComp, xlab,
                            ifelse(xlab == "auto","lambda",xlab)),
              type = "l", ...)
      tickat <- seq(1, length(lambdas[!colsWithNA]), length.out = 10)
      axis(3, at = lambdas[tickat],
           labels=apply(pars[regIndicators[regIndicators %in% rownames(pars)],] == targetVector[regIndicators[regIndicators %in% rownames(pars)]],2,sum)[tickat],
           outer= F,
           line=1,col="black",col.ticks="black",col.axis="black")
      mtext(paste0(paste0("Group ", group), ": # parameters on target"),3,line=3,at=mean(lambdas),col="black", cex = 1)

      readline(prompt="Press [Enter] for next plot... ")
    }
    par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
    return()
  }
  if(tolower(what) == "fit"){
    if(! "fit" %in% names(regCtsemObject)){
      stop("Could not find a field 'fit' in the regCtsem object'")
    }
    if(any(!criterion %in% rownames(regCtsemObject$fit))){
      stop(paste0(criterion, " not in fit values. Possible are: ", paste0(rownames(fit), collapse = ", "), "."))
    }

    color <- c("#008000", "#008080", "#800080", "#800000", rep("black", 5))
    lty <- c(1,2,3,4, rep(1,5))

    if(any(is.na(regCtsemObject$fit))){
      warning("NAs in regCtsemObject$fit Only plotting non-NA values")
    }

    colsWithNA <- apply(regCtsemObject$fit,2,function(x) any(is.na(x)))

    lambdas <- regCtsemObject$setup$lambdas

    matplot(x = lambdas[!colsWithNA], y = t(matrix(regCtsemObject$fit[criterion,!colsWithNA],
                                                   nrow = length(criterion))),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","fit values",ylab)),
            xlab = ifelse(skiptXlabComp, xlab,
                          ifelse(xlab == "auto",expression(lambda),xlab)),
            lty = lty[1:length(criterion)], col = color[1:length(criterion)],
            type = "l", lwd = 2, ...)

    legend(legend = criterion, "bottomright",
           lty = lty[1:length(criterion)],
           col = color[1:length(criterion)])
    return()
  }

  stop("Select what = 'drift', 'parameters', or 'fit' to generate a plot.")

}
