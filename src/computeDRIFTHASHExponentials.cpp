#include <RcppArmadillo.h>
#include "computeDiscreteParameters.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeDRIFTHASHExponentials
//'
//' Computes the DRIFTHASHExponential for the discrete time diffusion given the drifthash and DRIFTHASHExponential (list with discrete trait names, dTs and results)
//' The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
//' @param DRIFTHASH matrix with drift hash
//' @param DRIFTHASHExponentialUnique list with discrete trait names, dTs and results
//' @returns list with drift hash exponentials
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List computeDRIFTHASHExponentials(const arma::mat& DRIFTHASH, Rcpp::List DRIFTHASHExponentialUnique) {
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
