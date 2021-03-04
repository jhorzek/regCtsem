#include <RcppArmadillo.h>
#include "computeDiscreteParameters.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' Computes the computeDRIFTHASH given the drift. DRIFTHASH = DRIFT otimes I + I otimes DRIFT.
//'
//' @export
// [[Rcpp::export]]
arma::mat computeDRIFTHASH(arma::mat DRIFTValues) {

  arma::mat DRIFTHASH =  kron(DRIFTValues, eye( size(DRIFTValues) )) +  kron(eye( size(DRIFTValues) ),DRIFTValues);

  return(DRIFTHASH);
}
