#include <RcppArmadillo.h>
#include "getVarianceFromVarianceBase.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the covariance matrix given the covariancebase
using namespace Rcpp;
// [[Rcpp::export]]

arma::mat getVarianceFromVarianceBase(const arma::mat& varianceBaseValues){
  arma::mat varianceCholValues = arma::diagmat((arma::exp(arma::diagvec( varianceBaseValues )))) +
    varianceBaseValues - arma::diagmat((arma::diagvec( varianceBaseValues )));
  arma::mat varianceValues = varianceCholValues * arma::trans(varianceCholValues);
  return(varianceValues);
}
