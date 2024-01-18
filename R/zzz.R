#' regCtsem
#'
#' regCtsem uses objects from ctsemOMX and implements lasso, adaptive lasso and ridge regularization.
#'
#' @docType package
#' @author Jannik Orzek <orzek@mpib-berlin.mpg.de>
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @importFrom methods is
#' @useDynLib regCtsem
#' @name regCtsem

Rcpp::loadModule("cpptsemKalmanModel_cpp", TRUE)
Rcpp::loadModule("cpptsemRAMmodel_cpp", TRUE)

