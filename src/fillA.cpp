#include <RcppArmadillo.h>
#include "fillRAMMatrices.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Fills the A matrix (directed effects) given the A matrix, the discreteDRIFTUnique (list with discrete drift names, dTs and results),
// discreteTRAITUnique (list with discrete trait names, dTs and results), LAMBDA (loadings), and AParameterIndicators (tells fillA where to put the discreteDRIFT, discreteTRAIT, and LAMBDA)
// The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05

using namespace Rcpp;
// [[Rcpp::export]]

arma::mat fillA(arma::mat A,
                const bool hasDiscreteDRIFTUnique,
                const Rcpp::List& discreteDRIFTUnique,
                const bool hasDiscreteTRAITUnique,
                const Rcpp::List& discreteTRAITUnique,
                const arma::mat& LAMBDA,
                const Rcpp::DataFrame& AParameterIndicators){
  // extract labels of elements in A
  Rcpp::StringVector  AElementLabels = AParameterIndicators["label"];
  arma::uvec  AElementRowStart = AParameterIndicators["row_start"];
  arma::uvec  AElementRowEnd = AParameterIndicators["row_end"];
  arma::uvec  AElementColStart = AParameterIndicators["col_start"];
  arma::uvec  AElementColEnd = AParameterIndicators["col_end"];

  if(hasDiscreteDRIFTUnique){
    Rcpp::StringVector discreteDRIFTUniqueLabels = discreteDRIFTUnique["labels"];
    // setting discrete drifts:
    for(int i = 0; i < discreteDRIFTUniqueLabels.size(); i++){
      Rcpp::String discreteDRIFTUniqueLabel = discreteDRIFTUniqueLabels(i);
      arma::mat discreteDRIFTUniqueValues = discreteDRIFTUnique[discreteDRIFTUniqueLabel];

      // check for this discrete drift in AElementLabels:
      for(int j = 0; j < AElementLabels.size(); j++){

        if(discreteDRIFTUniqueLabels(i) == AElementLabels(j)){
          // set elements in A to the corresponding values
          A(arma::span(AElementRowStart(j),AElementRowEnd(j)), arma::span(AElementColStart(j), AElementColEnd(j)) ) = discreteDRIFTUniqueValues;
        }
      }
    }

  }
  if(hasDiscreteTRAITUnique){
    Rcpp::StringVector discreteTRAITUniqueLabels = discreteTRAITUnique["labels"];

    // setting discrete traits:
    for(int i = 0; i < discreteTRAITUniqueLabels.size(); i++){
      Rcpp::String discreteTRAITUniqueLabel = discreteTRAITUniqueLabels(i);
      arma::mat discreteTRAITUniqueValues = discreteTRAITUnique[discreteTRAITUniqueLabel];

      // check for this discrete trait in AElementLabels:
      for(int j = 0; j < AElementLabels.size(); j++){

        if(discreteTRAITUniqueLabels(i) == AElementLabels(j)){
          // set elements in A to the corresponding values
          A(arma::span(AElementRowStart(j),AElementRowEnd(j)), arma::span(AElementColStart(j), AElementColEnd(j)) ) = discreteTRAITUniqueValues;
        }
      }
    }
  }

  // setting Lambda:

  // check for LAMBDA in AElementLabels:
  Rcpp::String LAMBDALabel = "LAMBDA";
  for(int j = 0; j < AElementLabels.size(); j++){
    if(LAMBDALabel == AElementLabels(j)){
      // set elements in A to the corresponding values
      A(arma::span(AElementRowStart(j),AElementRowEnd(j)), arma::span(AElementColStart(j), AElementColEnd(j)) ) = LAMBDA;
    }
  }

  return(A);
}
