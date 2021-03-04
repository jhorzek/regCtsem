#include <RcppArmadillo.h>
#include "computeDiscreteParameters.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the discrete time diffusion matrices given the DRIFTHASHInverse, DIFFUSION,
// DRIFTHASHExponentialUnique (list with expm(drifthash)-names, dTs and results), and discreteDIFFUSIONUnique (list with diffusion-names, dTs and results)
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List computeDiscreteDIFFUSIONs(arma::mat DRIFTHASHInverse, arma::mat DIFFUSION,
                                     Rcpp::List DRIFTHASHExponentialUnique, Rcpp::List discreteDIFFUSIONUnique) {

  Rcpp::StringVector discreteDIFFUSIONUniqueLabels = discreteDIFFUSIONUnique["labels"];
  arma::vec dTDIFFUSION = discreteDIFFUSIONUnique["dT"];
  Rcpp::StringVector DRIFTHASHExponentialUniqueLabels = DRIFTHASHExponentialUnique["labels"];
  arma::vec dTDRIFTHASH = DRIFTHASHExponentialUnique["dT"];

  for(int i = 0; i < discreteDIFFUSIONUniqueLabels.size(); i++){
    if(!(dTDIFFUSION(i) == dTDRIFTHASH(i))){
      Rcpp::Rcout << "Error in computeDiscreteDIFFUSIONs: different time intervals!" << std::endl;
    }
    Rcpp::String discreteDIFFUSIONUniqueLabel = discreteDIFFUSIONUniqueLabels(i);
    arma::mat currentDiscreteDIFFUSION = discreteDIFFUSIONUnique[discreteDIFFUSIONUniqueLabel];
    unsigned int targetRows = currentDiscreteDIFFUSION.n_rows;
    unsigned int targetCols = currentDiscreteDIFFUSION.n_cols;
    Rcpp::String DRIFTHASHExponentialUniqueLabel = DRIFTHASHExponentialUniqueLabels(i);
    arma::mat currentDRIFTHASHExponential = DRIFTHASHExponentialUnique[DRIFTHASHExponentialUniqueLabel];

    currentDiscreteDIFFUSION = DRIFTHASHInverse * (currentDRIFTHASHExponential - arma::eye(arma::size(currentDRIFTHASHExponential)))*arma::vectorise(DIFFUSION);

    currentDiscreteDIFFUSION = arma::reshape( currentDiscreteDIFFUSION, targetRows, targetCols );
    // update DiscreteDIFFUSION
    discreteDIFFUSIONUnique[discreteDIFFUSIONUniqueLabel] = currentDiscreteDIFFUSION;
  }

  return(discreteDIFFUSIONUnique);
}


