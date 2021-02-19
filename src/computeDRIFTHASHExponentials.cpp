#include <RcppArmadillo.h>
#include "computeDiscreteParameters.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the DRIFTHASHExponential for the discrete time diffusion given the drifthash and DRIFTHASHExponential (list with discrete trait names, dTs and results)
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List computeDRIFTHASHExponentials(arma::mat DRIFTHASH, Rcpp::List DRIFTHASHExponentialUnique) {
  Rcpp::StringVector DRIFTHASHExponentialUniqueLabels = DRIFTHASHExponentialUnique["labels"];
  arma::vec dT = DRIFTHASHExponentialUnique["dT"];
  for(int i = 0; i < DRIFTHASHExponentialUniqueLabels.size(); i++){

    Rcpp::String DRIFTHASHExponentialUniqueLabel = DRIFTHASHExponentialUniqueLabels(i);
    arma::mat currentDRIFTHASHExponential = DRIFTHASHExponentialUnique[DRIFTHASHExponentialUniqueLabel];
    currentDRIFTHASHExponential = arma::expmat( DRIFTHASH * dT(i));

    // update DRIFTHASHExponentialUnique
    DRIFTHASHExponentialUnique[DRIFTHASHExponentialUniqueLabel] = currentDRIFTHASHExponential;
  }

  return(DRIFTHASHExponentialUnique);
}
