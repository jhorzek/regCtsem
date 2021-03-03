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
