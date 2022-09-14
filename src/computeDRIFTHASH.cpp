#include <RcppArmadillo.h>
#include "computeDiscreteParameters.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeDRIFTHASH
//'
//' Computes the computeDRIFTHASH given the drift. DRIFTHASH = DRIFT otimes I + I otimes DRIFT.
//'
//' The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
//'
//' @param DRIFTValues matrix with drift values
//' @returns matrix with drift hash
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat computeDRIFTHASH(const arma::mat& DRIFTValues) {

  arma::mat DRIFTHASH =  kron(DRIFTValues, eye( size(DRIFTValues) )) +  kron(eye( size(DRIFTValues) ),DRIFTValues);

  return(DRIFTHASH);
}
