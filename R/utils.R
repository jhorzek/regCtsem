#' checkNonConvexity
#'
#' checks the non-convexity of a provided model within the provided parameter bounds numerically.
#'
#' This is a very simple implementation of the concept described in Tamura, K., & Gallagher, M. (2019). Quantitative measure of nonconvexity for black-box continuous functions. Information Sciences, 476, 64–82. https://doi.org/10.1016/j.ins.2018.10.009
#' The basic idea is as follows: The function generates nSample points within the provided parameter bounds. For each combination of these sample points the following condition is evaluated:
#' f(x1) + f(x2) >= 2*f((x1 + x2)/2)  (midpoint convexity, see Equation 2 in Tamura, K., & Gallagher, M., 2019)
#' This condition should evaluate to TRUE for all possible combinations of x1 and x2. If it evaluates to FALSE, this implies that the model fitting function is not convex.
#' WARNING: Even if all evaluations indicate that midpoint convexity is given, this is no proof that the function is convex! However, if a counter-example is found (i.e., at least one evaluation results in FALSE), this implies
#' that the function might be nonconvex. Keep in mind that this is only a numerical check and susceptible to numerical errors
#' @param model model of type ctsemFit, Rcpp_cpptsemRAMmodel, Rcpp_cpptsemKalmanModel, or MxRAMModel
#' @param lowerBound vector of the same length as the number of parameters in the model specifying the lower bound for each parameter
#' @param upperBound vector of the same length as the number of parameters in the model specifying the upper bound for each parameter
#' @param nSample number of sample points between the lower and upper bound
#' @return returns the generated samplePoints, a logical vector indicating if the midpoint convexity was given for the comparison of these sample points (isConvex) and the difference f(x1) + f(x2) - 2*f((x1 + x2)/2) for each combination of sample points. This difference should be positive for all comparisons!
#' @return
checkNonConvexity <- function(model, lowerBound = NULL, upperBound = NULL, nSample = 50){
  if(any(class(model) == "ctsemFit")){
    model <- model$mxobj
    pars <- omxGetParameters(model)
    getFit <- function(mxObject, pars){
      mxObject <- omxSetParameters(mxObject, labels = names(pars), values = pars)
      mxObject <- mxRun(mxObject, useOptimizer = FALSE, silent = TRUE)
      return(mxObject$fitfunction$result[[1]])
    }
  }else if(any(class(model) == "Rcpp_cpptsemRAMmodel")){
    pars <- model$getParameterValues()
    getFit <- function(Rcpp_cpptsemRAMmodel, pars){
      Rcpp_cpptsemRAMmodel$setParameterValues(pars, names(pars))
      o <- try(Rcpp_cpptsemRAMmodel$computeRAM())
      if(any(class(o) == "try-error")){return(NA)}
      o <- try(Rcpp_cpptsemRAMmodel$fitRAM())
      if(any(class(o) == "try-error")){return(NA)}
      return(Rcpp_cpptsemRAMmodel$m2LL)
    }
  }else if(any(class(model) == "Rcpp_cpptsemKalmanModel")){
    pars <- model$getParameterValues()
    getFit <- function(Rcpp_cpptsemKalmanModel, pars){
      Rcpp_cpptsemKalmanModel$setParameterValues(pars, names(pars))
      o <- try(Rcpp_cpptsemKalmanModel$computeAndFitKalman())
      if(any(class(o) == "try-error")){return(NA)}
      return(Rcpp_cpptsemKalmanModel$m2LL)
    }
  }else if(any(class(model) == "MxRAMModel")){
    pars <- omxGetParameters(model)
    getFit <- function(mxObject, pars){
      mxObject <- omxSetParameters(mxObject, labels = names(pars), values = pars)
      mxObject <- mxRun(mxObject, useOptimizer = FALSE, silent = TRUE)
      return(mxObject$fitfunction$result[[1]])
    }
  }else{stop("model of unknown class")}

  if(is.null(lowerBound)){lowerBound <- pars - 3*abs(pars)}
  if(is.null(upperBound)){upperBound <- pars + 3*abs(pars)}

  # sample points between lower and upper bounds
  bounds <- cbind(lowerBound, upperBound)
  samplePoints <- apply(bounds, 1,
                        function(x) runif(n = nSample,
                                          min = x[1],
                                          max = x[2]))


  message("Computing fit at sample points: \n")
  fitAtSamplePoints <- c()
  pb <- txtProgressBar(min = 0, max = nSample, initial = 0, char = "=",
                       width = NA, title, label, style = 3, file = "")
  it <- 1
  for(samplePoint in 1:nSample){
    it <- it+1
    setTxtProgressBar(pb, it)
    fitAtSamplePoint <- try(getFit(model, pars = samplePoints[samplePoint,]), silent = T)
    if(any(class(fitAtSamplePoint) == "try-error") || is.na(fitAtSamplePoint)){
      observedNA <- TRUE
      fitAtSamplePoints[samplePoint] <- NA
      next}
    fitAtSamplePoints[samplePoint] <- fitAtSamplePoint
  }
  cat("\n")

  message("Evaluating midpoint convexity for all pairs of sample points: \n")
  isConvex <- c()
  elemNames <- c()
  midpointDifference <- c()
  observedNA <- FALSE
  maxChecks <- (sum(!is.na(fitAtSamplePoints))-1)*(sum(!is.na(fitAtSamplePoints)))/2
  pb <- txtProgressBar(min = 0, max = maxChecks, initial = 0, char = "=",
                       width = NA, title, label, style = 3, file = "")
  it <- 1
  for(samplePoint1 in 1:nSample){
    f_1 <- fitAtSamplePoints[samplePoint1]
    if(samplePoint1 == nSample){break}
    if(is.na(f_1)){next}

    for(samplePoint2 in (samplePoint1+1):nSample){
      it <- it+1
      setTxtProgressBar(pb, it)

      f_2 <- fitAtSamplePoints[samplePoint2]
      if(is.na(f_2)){next}
      f_3 <- try(getFit(model, pars = (samplePoints[samplePoint1,] + samplePoints[samplePoint2,])/2), silent = T)
      if(any(class(f_3) == "try-error")){
        observedNA <- TRUE
        next}
      if(anyNA(c(f_1,f_2,f_3))){
        observedNA <- TRUE
        next}
      # check midpoint convexity
      isConvex <- c(isConvex, f_1 + f_2 >= 2*f_3)
      elemNames <- c(elemNames, paste0(samplePoint1, " vs ", samplePoint2))
      midpointDifference <- c(midpointDifference, f_1 + f_2 - 2*f_3)
    }
  }
  cat("\n")
  if(observedNA){warning("Some sample points resulted in errors / NA. Only the results for non-NA evaluations are reported.")}
  names(isConvex) <- elemNames
  names(midpointDifference) <- elemNames
  rownames(samplePoints) <- 1:nSample

  message(paste0("Of in total ", length(isConvex), " non-NA checks ", sum(!isConvex), " resulted in non-convex evaluations at the mid point."))
  return(list("samplePoints" = samplePoints, "isConvex" = isConvex, "midpointDifference" = midpointDifference))
}

#' shinyfy
#'
#' creates a shiny application from a regCtsemObject
#' @param regCtsemObject regCtsemObject
#' @param runShiny boolean: if TRUE, the app is started directly. If false, the code to generate the app is returned
#' @param figHeight height of the figure
#' @export
shinyfy <- function(regCtsemObject, runShiny = TRUE, figHeight = 600){
  if(regCtsemObject$setup$autoCV){stop("Shinyfy currently not supported for models with automatic cross-validation.")}
  filename <- tempfile(pattern = "modelClone", tmpdir = tempdir(), fileext = ".RData")
  save(regCtsemObject, file = filename)
  lambdas <- regCtsemObject$setup$lambdas
  criteria <- c("m2LL","AIC", "BIC", "cvM2LL")
  criteria <- criteria[criteria %in% rownames(regCtsemObject$fit)]

  selCriterium <- paste0("list(", paste0(paste0("'", criteria , "' = '" , criteria, "'"), collapse = ", "), ")")
  selLambda <- paste0("list(", paste0(paste0("'", lambdas , "' = " , seq_len(length(lambdas))), collapse = ", "), ")")

  shinyCode <- paste0(
    'library(shiny)
library(ggplot2)
library(gridExtra)
library(regCtsem)
library(qgraph)
load("', filename, '")

latentNames <- regCtsemObject$setup$ctsemObject$ctmodelobj$latentNames
lambdas <- regCtsemObject$setup$lambdas
fit <- regCtsemObject$fit
parameters <- regCtsemObject$parameterEstimatesRaw
parameterNames <- rownames(parameters)
DRIFTlabels <- regCtsemObject$setup$mxObject$DRIFT$labels
DRIFTvalues <- regCtsemObject$setup$mxObject$DRIFT$values
DIFFUSIONbaseLabels <- regCtsemObject$setup$mxObject$DIFFUSIONbase$labels
DIFFUSIONbaseValues <- regCtsemObject$setup$mxObject$DIFFUSIONbase$values

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("regCtsem model inspection"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            selectInput("criterion", h3("criterion"), ', selCriterium,'
                        ),
            checkboxInput("discrete", "discrete", value = FALSE),
            numericInput("dT",
                         h3("time interval"),
                         value = 1),
            checkboxInput("manualLambda", "manually select lambda", value = FALSE),
            selectInput("lambdaVal", h3("lambda"),', selLambda,', selected = 1)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("qplot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$qplot <- renderPlot({
        makeEdges <- function(mat){
            from <- 1:ncol(mat)
            to <- 1:nrow(mat)
            weights <- as.vector(mat)
            edges <- matrix(nrow = length(weights), ncol = 3)
            colnames(edges) <- c("from", "to", "weight")
            edges[,"weight"] <- weights
            edges[,"from"] <- rep(from, each = length(to))
            edges[,"to"] <- rep(to, length(from))
            return(edges)
        }

        manualLambda <- input$manualLambda
        dT <- as.numeric(input$dT)
        criterion <- input$criterion
        lambdaVal <- as.integer(input$lambdaVal)
        discrete <- input$discrete

        ## check if manualLambda was set
        if(manualLambda){
            lambda <- lambdas[lambdaVal]
        }else{
            # else use criterion
            finalPars <- getFinalParameters(regCtsemObject, criterion = criterion, raw = TRUE)
            lambda <- finalPars$lambda
        }

        # set values
        for(DRIFTlabel in unique(c(DRIFTlabels[!is.na(DRIFTlabels)]))){
            DRIFTvalues[DRIFTlabels == DRIFTlabel & !is.na(DRIFTlabels)] <- parameters[DRIFTlabel,lambdas == lambda]
        }
        DRIFTHASH <- regCtsem::computeDRIFTHASH(DRIFTvalues)

        for(DIFFUSIONbaseLabel in unique(c(DIFFUSIONbaseLabels[!is.na(DIFFUSIONbaseLabels)]))){
            DIFFUSIONbaseValues[DIFFUSIONbaseLabels == DIFFUSIONbaseLabel & !is.na(DIFFUSIONbaseLabels)] <- parameters[DIFFUSIONbaseLabel,lambdas == lambda]
        }
        DIFFUSION <- regCtsem::getVarianceFromVarianceBase2(varianceBaseValues = DIFFUSIONbaseValues)

        if(discrete){

        discreteDRIFT <- expm(DRIFTvalues*dT)
        discreteDiffusion <- matrix(solve(DRIFTHASH) %*% (expm(DRIFTHASH*dT) - diag(nrow(DRIFTHASH)))%*%c(DIFFUSION),
                                    nrow = nrow(DIFFUSION), ncol = ncol(DIFFUSION))

        layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), heights = c(2,1))
        edges <- makeEdges(mat = discreteDRIFT)
        edge.color <- ifelse(edges[,"weight"]>0, "#008080", "#800000")
        qgraph(title = "discrete autoregressive and cross-lagged parameters",
               edges,
               maximum = .6,
               fade=FALSE,
               layout="circular",
               labels = latentNames,
               lty=ifelse(edges[,"weight"]>0,1,5),
               edge.labels=T, font = 2,
               edge.color=edge.color)

        edges <- makeEdges(mat = discreteDiffusion)
        edge.color <- ifelse(edges[,"weight"]>0, "#008080", "#800000")
        qgraph(title = "discrete diffusion",
               edges,
               maximum = .6,
               fade=FALSE,
               bidirectional = TRUE,
               layout="circular",
               labels = latentNames,
               lty=ifelse(edges[,"weight"]>0,1,5),
               edge.labels=T, font = 2,
               edge.color=edge.color)
        regCtsem::plot.regCtsem(regCtsemObject = regCtsemObject, what = "fit", criterion = criterion)
        }else{
        edges <- makeEdges(mat = DRIFTvalues)

        edge.color <- ifelse(edges[,"weight"]>0, "#008080", "#800000")
        layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), heights = c(2,1))
        qgraph(title = "DRIFT",
               edges,
               maximum = .6,
               fade=FALSE,
               layout="circular",
               labels = latentNames,
               lty=ifelse(edges[,"weight"]>0,1,5),
               edge.labels=T, font = 2,
               edge.color=edge.color)

        edges <- makeEdges(mat = DIFFUSION)
        edge.color <- ifelse(edges[,"weight"]>0, "#008080", "#800000")
        qgraph(title = "DIFFUSION",
               edges,
               maximum = .6,
               fade=FALSE,
               bidirectional = TRUE,
               layout="circular",
               labels = latentNames,
               lty=ifelse(edges[,"weight"]>0,1,5),
               edge.labels=T, font = 2,
               edge.color=edge.color)
        regCtsem::plot.regCtsem(regCtsemObject = regCtsemObject, what = "fit", criterion = criterion)

        }
    }, height = ', figHeight, ')
  }

  # Run the application
shinyApp(ui = ui, server = server)
  ')

if(runShiny){
  filename <- tempfile(pattern = "regCtsem_Shiny", tmpdir = tempdir(), fileext = ".R")
  fileConn <- file(filename)
  writeLines(shinyCode, fileConn)
  close(fileConn)

  shiny::runApp(filename)

}else{
  return(shinyCode)
}
}
