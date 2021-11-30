#include <RcppArmadillo.h>
#include "computeDiscreteParameters.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the discrete time trait matrices given the discreteDRIFTUnique (list with discrete drift names, dTs and results) and discreteTRAITUnique (list with discrete trait names, dTs and results)
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List computeDiscreteTRAITs(const Rcpp::List& discreteDRIFTUnique, Rcpp::List discreteTRAITUnique) {

  Rcpp::StringVector discreteTRAITUniqueLabels = discreteTRAITUnique["labels"];
  arma::vec dTTRAIT = discreteTRAITUnique["dT"];

  Rcpp::StringVector discreteDRIFTUniqueLabels = discreteDRIFTUnique["labels"];
  arma::vec dTDRIFT = discreteDRIFTUnique["dT"];

  for(int i = 0; i < discreteTRAITUniqueLabels.size(); i++){
    // check if discreteDRIFT and discreteTRAIT habe same dT
    if(!(dTTRAIT(i) == dTDRIFT(i))){
      Rcpp::Rcout << "Error! Non-matching time intervals in computeDiscreteTRAITs" << std::endl;
    }

    Rcpp::String discreteTRAITUniqueLabel = discreteTRAITUniqueLabels(i);
    arma::mat currentDiscreteTRAITValues = discreteTRAITUnique[discreteTRAITUniqueLabel];

    Rcpp::String discreteDRIFTUniqueLabel = discreteDRIFTUniqueLabels(i);
    arma::mat discreteDriftValues = discreteDRIFTUnique[discreteDRIFTUniqueLabel];

    arma::mat diagMat = arma::eye( arma::size(discreteDriftValues) );

    currentDiscreteTRAITValues = diagMat - discreteDriftValues;

    // update discreteDrift
    discreteTRAITUnique[discreteTRAITUniqueLabel] = currentDiscreteTRAITValues;
  }

  return(discreteTRAITUnique);
}


