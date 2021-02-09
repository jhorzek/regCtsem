#' generateDRIFTPlot
#'
#' generates plot of drift values
#'
#' @param model model from regCtsem
#' @param ylab label for y-axis. auto will set to lambda
#' @param xlab label for x-axis. auto will set to drift
#' @param skiptYlabComp boolean: is ylab a string
#' @param skiptXlabComp boolean: is xlab a string
#' @export
generateDRIFTPlot <- function(model, ylab = "auto", xlab = "auto", skiptYlabComp, skiptXlabComp){
  drifts <- model$parameters[grepl("drift", rownames(model$parameters)),]
  regIndicators <- model$setup$regIndicators
  drifts_regularized <- drifts[regIndicators,]
  regValues <- model$setup$regValues

  colsWithNA <- apply(model$parameters,2,function(x) any(is.na(x)))

  color <- ifelse(rownames(drifts) %in% regIndicators, yes = "white", "black")

  par(mar=c(4, 3, 5, 2), xpd=TRUE)

  matplot(x = regValues[!colsWithNA], t(drifts[!colsWithNA]), lty = 2,
          lwd = 2, type = "l",
          col = color, xlab = ifelse(skiptXlabComp, xlab,
                                     ifelse(xlab == "auto",expression(lambda),ylab)),
          ylab = ifelse(skiptYlabComp, ylab,
                        ifelse(ylab == "auto","drift",ylab)))
  matplot(x = regValues[!colsWithNA], t(drifts_regularized[!colsWithNA]), lty = 1, lwd = 2, type = "l", col = "#2166AC", add = TRUE)
  tickat <- seq(1, length(regValues[!colsWithNA]), length.out = 10)
  axis(3, at = regValues[tickat],
       labels=apply(drifts == 0,2,sum)[tickat],
       outer= F,
       line=1,col="black",col.ticks="black",col.axis="black")
  mtext("# zeroed parameters",3,line=3,at=mean(regValues),col="black", cex = 1)

  par(mar=c(5, 4, 4, 2) + 0.1, xpd=TRUE)
}

