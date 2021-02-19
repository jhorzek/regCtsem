extractCtsemMatrices <- function(mxObject, nlatent, nmanifest){
  ctMatrices <- list()
  ## latent
  # T0MEANS
  cat("Identified ")
  if(!is.null(mxObject$T0MEANS)){
    cat("T0MEANS, ")
    T0MEANS <- list("values" = NULL, "names" = NULL)
    T0MEANS$values <- deepCopyNumericMatrix(mxObject$T0MEANS$values)
    T0MEANS$names <- deepCopyStringMatrix(mxObject$T0MEANS$labels)
    ctMatrices[["T0MEANS"]] <- T0MEANS
  }

  # T0VARbase
  if(!is.null(mxObject$T0VARbase)){
    cat("T0VARbase, ")
    T0VARbase <- list("values" = NULL, "names" = NULL)
    T0VARbase$values <- deepCopyNumericMatrix(mxObject$T0VARbase$values)
    T0VARbase$names <- deepCopyStringMatrix(mxObject$T0VARbase$labels)
    ctMatrices[["T0VARbase"]] <- T0VARbase
  }


  # DRIFT
  if(!is.null(mxObject$DRIFT)){
    cat("DRIFT, ")
    DRIFT <- list("values" = NULL, "names" = NULL)
    DRIFT$values <- deepCopyNumericMatrix(mxObject$DRIFT$values)
    DRIFT$names <- deepCopyStringMatrix(mxObject$DRIFT$labels)
    ctMatrices[["DRIFT"]] <- DRIFT
  }


  # DIFFUSIONbase
  if(!is.null(mxObject$DIFFUSIONbase)){
    cat("DIFFUSIONbase, ")
    DIFFUSIONbase <- list("values" = NULL, "names" = NULL)
    DIFFUSIONbase$values <- deepCopyNumericMatrix(mxObject$DIFFUSIONbase$values)
    DIFFUSIONbase$names <- deepCopyStringMatrix(mxObject$DIFFUSIONbase$labels)
    ctMatrices[["DIFFUSIONbase"]] <- DIFFUSIONbase
  }

  # TRAITVARbase
  if(!is.null(mxObject$TRAITVARbase$values)){
    cat("TRAITVARbase, ")
    TRAITVARbase <- list("values" = NULL, "names" = NULL)
    TRAITVARbase$values <- deepCopyNumericMatrix(mxObject$TRAITVARbase$values)
    TRAITVARbase$names <- deepCopyStringMatrix(mxObject$TRAITVARbase$labels)
    ctMatrices[["TRAITVARbase"]] <- TRAITVARbase
  }

  # T0TRAITEFFECT
  if(!is.null(mxObject$T0TRAITEFFECT)){
    cat("T0TRAITEFFECT, ")
    T0TRAITEFFECT <- list("values" = NULL, "names" = NULL)
    T0TRAITEFFECT$values <- deepCopyNumericMatrix(mxObject$T0TRAITEFFECT$values)
    T0TRAITEFFECT$names <- deepCopyStringMatrix(mxObject$T0TRAITEFFECT$labels)
    ctMatrices[["T0TRAITEFFECT"]] <- T0TRAITEFFECT
  }

  # CINT
  if(!is.null(mxObject$CINT)){
    cat("CINT, ")
    CINT <- list("values" = NULL, "names" = NULL)
    CINT$values <- deepCopyNumericMatrix(mxObject$CINT$values)
    CINT$names <- deepCopyStringMatrix(mxObject$CINT$labels)
    ctMatrices[["CINT"]] <- CINT
  }

  ## manifest

  # MANIFESTMEANS
  if(!is.null(mxObject$MANIFESTMEANS)){
    cat("MANIFESTMEANS, ")
    MANIFESTMEANS <- list("values" = NULL, "names" = NULL)
    MANIFESTMEANS$values <- deepCopyNumericMatrix(mxObject$MANIFESTMEANS$values)
    MANIFESTMEANS$names <- deepCopyStringMatrix(mxObject$MANIFESTMEANS$labels)
    ctMatrices[["MANIFESTMEANS"]] <- MANIFESTMEANS
  }

  # LAMBDA
  if(!is.null(mxObject$LAMBDA)){
    cat("LAMBDA, ")
  LAMBDA <- list("values" = NULL, "names" = NULL)
  LAMBDA$values <- deepCopyNumericMatrix(mxObject$LAMBDA$values)
  LAMBDA$names <- deepCopyStringMatrix(mxObject$LAMBDA$labels)
  ctMatrices[["LAMBDA"]] <- LAMBDA
  }

  # MANIFESTVAR
  if(!is.null(mxObject$MANIFESTVARbase$values)){
    cat("MANIFESTVARbase ")
    MANIFESTVARbase <- list("values" = NULL, "names" = NULL)
    MANIFESTVARbase$values <- deepCopyNumericMatrix(mxObject$MANIFESTVARbase$values)
    MANIFESTVARbase$names <- deepCopyStringMatrix(mxObject$MANIFESTVARbase$labels)
    ctMatrices[["MANIFESTVARbase"]] <- MANIFESTVARbase
  }

  cat("\n")
  return(ctMatrices)

}

