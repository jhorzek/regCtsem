#' plot.regCtsem
#'
#' @param x fitted regCtsem object
#' @param y NULL
#' @param ... additional parameters for plot or matplot (e.g., lty, col, ...)
#' @param what what should be plotted? Possbile are: 'drift', 'parameters', 'fit' or labels of specific parameters
#' @param criterion vector with labels of criteria which should be plotted when what = 'fit'
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @author Jannik Orzek
#' @export
plot.regCtsem <- function(x, y = NULL, ..., what = "drift", criterion = "BIC", xlab = "auto", ylab = "auto"){
  if(is.character(xlab)){skiptXlabComp <- F}else{skiptXlabComp <- T}
  if(is.character(ylab)){skiptYlabComp <- F}else{skiptYlabComp <- T}

  if(what == "drift"){
    what <- c(x$setup$ctsemObject$mxobj$DRIFT$labels)
  }

  ## Plot specific parameter values
  if(all(what %in% rownames(x$parameters))){
    pars <- x$parameters
    parsSelected <- subset(pars, subset = rownames(pars) %in% what)
    regIndicators <- x$setup$regIndicators
    targetVector <- x$setup$targetVector
    pars_regularized <- subset(parsSelected, subset = rownames(parsSelected) %in% regIndicators)

    lambdas <- x$setup$lambdas

    colsWithNA <- apply(pars,2,function(x) any(is.na(x)))

    color <- ifelse(rownames(pars) %in% regIndicators, yes = "white", "black")

    graphics::par(mar=c(4, 3, 5, 2), xpd=TRUE)

    graphics::matplot(x = lambdas[!colsWithNA], t(subset(parsSelected, select = !colsWithNA)), lty = 2,
            lwd = 2, type = "l",
            col = color, xlab = ifelse(skiptXlabComp, xlab,
                                       ifelse(xlab == "auto",expression(lambda),ylab)),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto", "Values",ylab)))
    graphics::matplot(x = lambdas[!colsWithNA], t(subset(pars_regularized, select = !colsWithNA)), lty = 1, lwd = 2, type = "l", col = "#2166AC", add = TRUE)
    tickat <- seq(1, length(lambdas[!colsWithNA]), length.out = 10)
    graphics::axis(3, at = lambdas[tickat],
         labels=apply(pars[regIndicators,] == targetVector[regIndicators],2,sum)[tickat],
         outer= F,
         line=1,col="black",col.ticks="black",col.axis="black")
    graphics::mtext("# parameters on target",3,line=3,at=mean(lambdas),col="black", cex = 1)

    graphics::par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
    return(invisible())
  }

  if(tolower(what) == "drift"){
    pars <- x$parameters
    regIndicators <- x$setup$regIndicators
    targetVector <- x$setup$targetVector
    drifts <- subset(pars, subset = grepl("drift", rownames(pars)))
    drifts_regularized <- subset(drifts, rownames(drifts) %in% regIndicators)
    lambdas <- x$setup$lambdas

    colsWithNA <- apply(pars,2,function(x) any(is.na(x)))

    color <- ifelse(rownames(drifts) %in% regIndicators, yes = "white", "black")

    graphics::par(mar=c(4, 3, 5, 2), xpd=TRUE)

    graphics::matplot(x = lambdas[!colsWithNA], t(subset(drifts, select = !colsWithNA)), lty = 2,
            lwd = 2, type = "l",
            col = color, xlab = ifelse(skiptXlabComp, xlab,
                                       ifelse(xlab == "auto",expression(lambda),ylab)),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","drift",ylab)))
    graphics::matplot(x = lambdas[!colsWithNA], t(subset(drifts_regularized, select = !colsWithNA)), lty = 1, lwd = 2, type = "l", col = "#2166AC", add = TRUE)
    tickat <- seq(1, length(lambdas[!colsWithNA]), length.out = 10)
    graphics::axis(3, at = lambdas[tickat],
         labels=apply(pars[regIndicators,] == targetVector[regIndicators],2,sum)[tickat],
         outer= F,
         line=1,col="black",col.ticks="black",col.axis="black")
    graphics::mtext("# parameters on target",3,line=3,at=mean(lambdas),col="black", cex = 1)

    graphics::par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
    return(invisible())
  }
  if(tolower(what) == "parameters"){
    if(!"parameters" %in% names(x)){stop("Plot of drift values not possible for cross-validation. Use what = 'fit' to plot the fit.")}

    if(any(is.na(x$parameters))){
      warning("NAs in x$parameters. Only plotting non-NA values")
    }

    colsWithNA <- apply(x$parameters,2,function(x) any(is.na(x)))

    lambdas <- x$setup$lambdas

    graphics::matplot(x = lambdas[!colsWithNA], y = t(subset(x$parameters, select = !colsWithNA)),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","parameter values",ylab)),
            xlab = ifelse(skiptXlabComp, xlab,
                          ifelse(xlab == "auto","lambda",xlab)),
            type = "l", ...)
    plotted <- TRUE
  }
  if(tolower(what) == "fit"){
    if(! "fit" %in% names(x)){
      stop("Could not find a field 'fit' in the regCtsem object'")
    }
    if(any(!criterion %in% rownames(x$fit))){
      stop(paste0(criterion, " not in fit values. Possible are: ", paste0(rownames(fit), collapse = ", "), "."))
    }

    color <- c("#008000", "#008080", "#800080", "#800000", rep("black", 5))
    lty <- c(1,2,3,4, rep(1,5))

    if(any(is.na(x$fit))){
      warning("NAs in x$fit Only plotting non-NA values")
    }

    colsWithNA <- apply(x$fit,2,function(x) any(is.na(x)))

    lambdas <- x$setup$lambdas

    graphics::matplot(x = lambdas[!colsWithNA], y = t(matrix(x$fit[criterion,!colsWithNA],
                                                   nrow = length(criterion))),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","fit values",ylab)),
            xlab = ifelse(skiptXlabComp, xlab,
                          ifelse(xlab == "auto",expression(lambda),xlab)),
            lty = lty[1:length(criterion)], col = color[1:length(criterion)],
            type = "l", lwd = 2, ...)

    graphics::legend(legend = criterion, "bottomright",
           lty = lty[1:length(criterion)],
           col = color[1:length(criterion)])
    return(invisible())
  }

  stop("Select what = 'drift', 'parameters', or 'fit' to generate a plot.")

}

#' plot.regCtsemCV
#'
#' @param x fitted regCtsem object
#' @param y NULL
#' @param ... additional parameters for plot or matplot (e.g., lty, col, ...)
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @author Jannik Orzek
#' @export
plot.regCtsemCV <- function(x, y = NULL, ..., xlab = "auto", ylab = "auto"){
  if(is.character(xlab)){skiptXlabComp <- F}else{skiptXlabComp <- T}
  if(is.character(ylab)){skiptYlabComp <- F}else{skiptYlabComp <- T}

  criterion <- "mean"
  lambdas <- x$setup$lambdas

  plot(x = lambdas, y = x$fit[criterion,],
       ylab = ifelse(skiptYlabComp, ylab,
                     ifelse(ylab == "auto","mean CV fit",ylab)),
       xlab = ifelse(skiptXlabComp, xlab,
                     ifelse(xlab == "auto",expression(lambda),xlab)),
       type = "l", lwd = 2, col = "#008080", ...)

}

#' plot.regCtsemMultiSubject
#'
#' @param x fitted regCtsem object
#' @param y NULL
#' @param ... additional parameters for plot or matplot (e.g., lty, col, ...)
#' @param what what should be plotted? Possbile are: 'drift', 'parameters', 'fit'
#' @param groups For which groups should the parameters be plotted (numeric vector)?
#' @param criterion vector with labels of criteria which should be plotted when what = 'fit'
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @author Jannik Orzek
#' @export
plot.regCtsemMultiSubject <- function(x, y = NULL, ..., what = "drift", groups = NULL, criterion = "BIC", xlab = "auto", ylab = "auto"){
  if(is.character(xlab)){skiptXlabComp <- F}else{skiptXlabComp <- T}
  if(is.character(ylab)){skiptYlabComp <- F}else{skiptYlabComp <- T}
  if(is.null(groups)){
    groups <- 1:length(x$parameters)
  }

  ## Plot specific parameter values
  # if(any(c(what, paste0(what, "_G1")) %in% rownames(x$parameters[[1]]))){
  #   for(group in groups){
  #     pars <- x$parameters[[group]]
  #     parsSelected <- subset(pars, subset = rownames(pars) %in% what)
  #     targetVector <- x$setup$targetVector
  #     regIndicators <- x$setup$regIndicators
  #     pars_regularized <- subset(pars, subset = rownames(pars) %in% regIndicators)
  #
  #     lambdas <- x$setup$lambdas
  #
  #     colsWithNA <- apply(pars,2,function(x) any(is.na(x)))
  #
  #     color <- ifelse(rownames(pars) %in% regIndicators, yes = "white", "black")
  #
  #     graphics::par(mar=c(4, 3, 5, 2), xpd=TRUE)
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
  #   graphics::par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
  #   return()
  # }

  if(tolower(what) == "drift"){
    for(group in groups){
      pars <- x$parameters[[group]]
      drifts <- subset(pars, subset = grepl("drift", rownames(pars)))
      regIndicators <- x$setup$regIndicators
      targetVector <- x$setup$targetVector
      drifts_regularized <- subset(drifts, rownames(drifts) %in% regIndicators)
      lambdas <- x$setup$lambdas

      colsWithNA <- apply(pars,2,function(x) any(is.na(x)))

      color <- ifelse(rownames(drifts) %in% regIndicators, yes = "white", "black")

      graphics::par(mar=c(4, 3, 5, 2), xpd=TRUE)

      graphics::matplot(x = lambdas[!colsWithNA], t(subset(drifts, select = !colsWithNA)), lty = 2,
              lwd = 2, type = "l",
              col = color, xlab = ifelse(skiptXlabComp, xlab,
                                         ifelse(xlab == "auto",expression(lambda),ylab)),
              ylab = ifelse(skiptYlabComp, ylab,
                            ifelse(ylab == "auto","drift",ylab)))
      graphics::matplot(x = lambdas[!colsWithNA], t(subset(drifts_regularized, select = !colsWithNA)), lty = 1, lwd = 2, type = "l", col = "#2166AC", add = TRUE)
      tickat <- seq(1, length(lambdas[!colsWithNA]), length.out = 10)
      graphics::axis(3, at = lambdas[tickat],
           labels=apply(pars[regIndicators[regIndicators %in% rownames(pars)],] == targetVector[regIndicators[regIndicators %in% rownames(pars)]],2,sum)[tickat],
           outer= F,
           line=1,col="black",col.ticks="black",col.axis="black")
      graphics::mtext(paste0(paste0("Group ", group), ": # parameters on target"),3,line=3,at=mean(lambdas),col="black", cex = 1)

      readline(prompt="Press [Enter] for next plot... ")
    }
    graphics::par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
    return(invisible())
  }
  if(tolower(what) == "parameters"){
    if(!"parameters" %in% names(x)){stop("Plot of drift values not possible for cross-validation. Use what = 'fit' to plot the fit.")}

    if(any(is.na(x$parameters[[group]]))){
      warning("NAs in x$parameters. Only plotting non-NA values")
    }
    for(group in groups){
      pars <- x$parameters[[group]]
      regIndicators <- x$setup$regIndicators
      targetVector <- x$setup$targetVector

      graphics::par(mar=c(4, 3, 5, 2), xpd=TRUE)
      colsWithNA <- apply(x$parameters[[group]],2,function(x) any(is.na(x)))

      lambdas <- x$setup$lambdas

      graphics::matplot(x = lambdas[!colsWithNA], y = t(subset(x$parameters[[group]], select = !colsWithNA)),
              ylab = ifelse(skiptYlabComp, ylab,
                            ifelse(ylab == "auto","parameter values",ylab)),
              xlab = ifelse(skiptXlabComp, xlab,
                            ifelse(xlab == "auto","lambda",xlab)),
              type = "l", ...)
      tickat <- seq(1, length(lambdas[!colsWithNA]), length.out = 10)
      graphics::axis(3, at = lambdas[tickat],
           labels=apply(pars[regIndicators[regIndicators %in% rownames(pars)],] == targetVector[regIndicators[regIndicators %in% rownames(pars)]],2,sum)[tickat],
           outer= F,
           line=1,col="black",col.ticks="black",col.axis="black")
      graphics::mtext(paste0(paste0("Group ", group), ": # parameters on target"),3,line=3,at=mean(lambdas),col="black", cex = 1)

      readline(prompt="Press [Enter] for next plot... ")
    }
    graphics::par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
    return(invisible())
  }
  if(tolower(what) == "fit"){
    if(! "fit" %in% names(x)){
      stop("Could not find a field 'fit' in the regCtsem object'")
    }
    if(any(!criterion %in% rownames(x$fit))){
      stop(paste0(criterion, " not in fit values. Possible are: ", paste0(rownames(fit), collapse = ", "), "."))
    }

    color <- c("#008000", "#008080", "#800080", "#800000", rep("black", 5))
    lty <- c(1,2,3,4, rep(1,5))

    if(any(is.na(x$fit))){
      warning("NAs in x$fit Only plotting non-NA values")
    }

    colsWithNA <- apply(x$fit,2,function(x) any(is.na(x)))

    lambdas <- x$setup$lambdas

    graphics::matplot(x = lambdas[!colsWithNA], y = t(matrix(x$fit[criterion,!colsWithNA],
                                                   nrow = length(criterion))),
            ylab = ifelse(skiptYlabComp, ylab,
                          ifelse(ylab == "auto","fit values",ylab)),
            xlab = ifelse(skiptXlabComp, xlab,
                          ifelse(xlab == "auto",expression(lambda),xlab)),
            lty = lty[1:length(criterion)], col = color[1:length(criterion)],
            type = "l", lwd = 2, ...)

    graphics::legend(legend = criterion, "bottomright",
           lty = lty[1:length(criterion)],
           col = color[1:length(criterion)])
    return(invisible())
  }

  stop("Select what = 'drift', 'parameters', or 'fit' to generate a plot.")

}
