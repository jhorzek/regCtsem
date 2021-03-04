#include <RcppArmadillo.h>
#include "getMatrixExponential.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the matrix exponential
using namespace Rcpp;
// [[Rcpp::export]]

arma::mat getMatrixExponential(arma::mat mat){
  return(arma::expmat(mat));
}
