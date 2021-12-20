#include <RcppArmadillo.h>
#include "computeDiscreteParameters.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the discrete time diffusion matrices given the DRIFTHASHInverse, DIFFUSION,
// DRIFTHASHExponentialUnique (list with expm(drifthash)-names, dTs and results), and discreteDIFFUSIONUnique (list with diffusion-names, dTs and results)
// The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
using namespace Rcpp;
// [[Rcpp::export]]
Rcpp::List computeDiscreteDIFFUSIONs(const arma::mat& DRIFTHASHInverse, const arma::mat& DIFFUSION,
                                     const Rcpp::List& DRIFTHASHExponentialUnique, Rcpp::List discreteDIFFUSIONUnique) {

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


