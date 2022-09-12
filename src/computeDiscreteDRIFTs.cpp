#include <RcppArmadillo.h>
#include "computeDiscreteParameters.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeDiscreteDRIFTs
//'
//' Computes the discrete time drift matrices given the drift and discreteDRIFTUnique (list with discrete drift names, dTs and results)
//' The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
//' @param DRIFTValues matrix with drift values
//' @param discreteDRIFTUnique list with discrete drift names, dTs and results
//' @returns list
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List computeDiscreteDRIFTs(const arma::mat& DRIFTValues, Rcpp::List discreteDRIFTUnique) {

  Rcpp::StringVector discreteDRIFTUniqueLabels = discreteDRIFTUnique["labels"];
  arma::vec dT = discreteDRIFTUnique["dT"];

  for(int i = 0; i < discreteDRIFTUniqueLabels.size(); i++){

    Rcpp::String discreteDRIFTUniqueLabel = discreteDRIFTUniqueLabels(i);
    arma::mat currentDiscreteDrift = discreteDRIFTUnique[discreteDRIFTUniqueLabel];

    currentDiscreteDrift = arma::expmat( DRIFTValues * dT(i));

    // update discreteDrift
    discreteDRIFTUnique[discreteDRIFTUniqueLabel] = currentDiscreteDrift;
  }

  return(discreteDRIFTUnique);
}


