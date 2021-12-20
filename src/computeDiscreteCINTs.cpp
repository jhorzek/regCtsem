#include <RcppArmadillo.h>
#include "computeDiscreteParameters.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the discrete time latent intercept values (list with discreteCINT-names, dTs and results) given the inverse of the drift (DRIFTInverseValues), the discreteDRIFTUnique (list with discreteDRIFT-names, dTs and results) and CINTValues
// The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List computeDiscreteCINTs(Rcpp::List discreteCINTUnique, const arma::mat& DRIFTInverseValues, const Rcpp::List& discreteDRIFTUnique, const arma::colvec& CINTValues) {

  Rcpp::StringVector discreteCINTUniqueLabels = discreteCINTUnique["labels"];
  arma::vec dtCINT = discreteCINTUnique["dT"];
  Rcpp::StringVector discreteDRIFTUniqueLabels = discreteDRIFTUnique["labels"];
  arma::vec dtDRIFT = discreteDRIFTUnique["dT"];

  for(int i = 0; i < discreteCINTUniqueLabels.size(); i++){
    if(!(dtCINT(i) == dtDRIFT(i))){
      Rcpp::Rcout << "Error in computeDiscreteCINTs: different time intervals!" << std::endl;
    }
    Rcpp::String discreteCINTUniqueLabel = discreteCINTUniqueLabels(i); // get specific discreteCINT
    arma::colvec currentDiscreteCINT = discreteCINTUnique[discreteCINTUniqueLabel]; // the results will be saved here

    Rcpp::String discreteDRIFTUniqueLabel = discreteDRIFTUniqueLabels(i);
    arma::mat currentDiscreteDRIFT = discreteDRIFTUnique[discreteDRIFTUniqueLabel];

    currentDiscreteCINT = DRIFTInverseValues * (currentDiscreteDRIFT - arma::eye(arma::size(currentDiscreteDRIFT)))*CINTValues;

    // update DiscreteCINTUnique
    discreteCINTUnique[discreteCINTUniqueLabel] = currentDiscreteCINT;
  }

  return(discreteCINTUnique);
}


