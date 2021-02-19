#include <RcppArmadillo.h>
#include "computeDiscreteParameters.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the discrete time drift matrices given the drift and discreteDRIFTUnique (list with discrete drift names, dTs and results)
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List computeDiscreteDRIFTs(arma::mat DRIFTValues, Rcpp::List discreteDRIFTUnique) {

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


