#' exact_initialHessian
#'
#' computes an initial positive-definite Hessian matrix
#'
#' @param mxObject Fitted object of class MxObject
#' @param approximationType which Hessian should be used? Currently available are "ident" for an identity matrix and "OpenMx" for the Hessian from OpenMx Gradient Descent (recommended)
#' @param estimatedHessian estimated Hessian from OpenMx
#' @author Jannik Orzek
#' @import ctsemOMX
#' @export
exact_initialHessian <- function(mxObject,
                                 approximationType,
                                 estimatedHessian){

  numberParameters <- nrow(estimatedHessian)

  # identity matrix
  if(approximationType == "ident"){
    return(diag(numberParameters))
  }

  if(approximationType == "OpenMx"){
    #positiveDefiniteHessianModel <- suppressWarnings(try(mxOption(model = mxObject, key = "Calculate Hessian", value = "No")))
    #positiveDefiniteHessianModel <- suppressWarnings(try(mxRun(positiveDefiniteHessianModel, useOptimizer = TRUE, silent = TRUE)))

    hessian <- mxObject$output$hessian

    if(is.null(hessian)){
      # if the hessian is not returned in the mxObject:
      hessian <- regCtsem::exact_initialHessian(mxObject = mxObject,
                                                 approximationType = "ident",
                                                 estimatedHessian = estimatedHessian)
    }

    return(hessian)
  }

  # approximation from Nocedal, p. 200
  if(approximationType == "nocedal"){
    initialHessian <- (t(y_1)%*%s_1)/(t(y_1)%*%y_1)%*%diag(numberParameters)
    return(initialHessian)
  }

  stop("Undefined approximation type")
}



