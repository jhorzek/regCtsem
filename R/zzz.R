#' regCtsem
#'
#' This package is based on ctsem and can regularize Drift parameters in continuous time structural equation models
#'
#' @docType package
#' @author Jannik Orzek <orzek@mpib-berlin.mpg.de>
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @useDynLib regCtsem
#' @name regCtsem

Rcpp::loadModule("cpptsemKalmanModel_cpp", TRUE)
Rcpp::loadModule("cpptsemRAMmodel_cpp", TRUE)

