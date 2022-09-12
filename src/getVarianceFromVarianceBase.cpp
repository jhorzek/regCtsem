#include <RcppArmadillo.h>
#include "getVarianceFromVarianceBase.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' getVarianceFromVarianceBase
//'
//' Computes the covariance matrix given the covariancebase
//' The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
//' @param varianceBaseValues matrix with base values
//' @keywords internal
// [[Rcpp::export]]
arma::mat getVarianceFromVarianceBase(const arma::mat& varianceBaseValues){
  arma::mat varianceCholValues = arma::diagmat((arma::exp(arma::diagvec( varianceBaseValues )))) +
    varianceBaseValues - arma::diagmat((arma::diagvec( varianceBaseValues )));
  arma::mat varianceValues = varianceCholValues * arma::trans(varianceCholValues);
  return(varianceValues);
}
