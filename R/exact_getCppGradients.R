#' exact_getCppGradients
#'
#' exact_getCppGradients will try to compute gradients with decreasing precision starting from the default in OpenMx. Sometimes the gradients will result in NA for a very specific setting of the precisions. Then it can help to slightly alter the precision
#'
#' @param cppmodel model of type cpptsem
#' @param objective ml or Kalman
#' @export
exact_getCppGradients <- function(cppmodel, objective){
  # will try different precisions for the gradients
  defaultPrecision <- OpenMx::imxAutoOptionValue("Gradient step size")
  defaultPrecision2 <- (1.1 * 10^(-16))^(1/3)
  precisions <- rev(seq(defaultPrecision, defaultPrecision2, length.out = 5))
  for(precision in precisions){
    if(tolower(objective) == "ml"){
      invisible(capture.output(gradients <- try(cppmodel$approxRAMGradients(precision), silent = T), type = "message"))
    }else{
      invisible(capture.output(gradients <- try(cppmodel$approxKalmanGradients(precision), silent = T), type = "message"))
    }
    if(!(any(class(gradients) == "try-error")) &
       !anyNA(gradients)){
      break
    }
  }
  return(gradients)
}




