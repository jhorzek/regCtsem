#include <RcppArmadillo.h>
#include "fillRAMMatrices.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Fills the M vector (intercepts) given the M vector, the MANIFESTMEANS,
// discreteCINTUnique (list with discrete continuous time intercepts, dTs and results), and cppMParameterIndicators (tells fillM where to put the MANIFESTMEANS and discreteCINT)

using namespace Rcpp;
// [[Rcpp::export]]

arma::colvec fillM(arma::colvec M, arma::colvec T0MEANS, arma::colvec MANIFESTMEANS,
                   bool hasDiscreteCINTUnique, Rcpp::List discreteCINTUnique,
                   Rcpp::DataFrame cppMParameterIndicators){
  // extract labels of elements in A
  Rcpp::StringVector  MElementLabels = cppMParameterIndicators["label"];
  arma::uvec  MElementRowStart = cppMParameterIndicators["row_start"];
  arma::uvec  MElementRowEnd = cppMParameterIndicators["row_end"];
  arma::uvec  MElementColStart = cppMParameterIndicators["col_start"];
  arma::uvec  MElementColEnd = cppMParameterIndicators["col_end"];

  if(hasDiscreteCINTUnique){
    Rcpp::StringVector discreteCINTUniqueLabels = discreteCINTUnique["labels"];
    // setting discrete cints:
    for(int i = 0; i < discreteCINTUniqueLabels.size(); i++){
      Rcpp::String discreteCINTUniqueLabel = discreteCINTUniqueLabels(i);
      arma::colvec discreteCINTUniqueValues = Rcpp::as<arma::colvec> (discreteCINTUnique[discreteCINTUniqueLabel]);

      // check for this discrete drift in MElementLabels:
      for(int j = 0; j < MElementLabels.size(); j++){

        if(discreteCINTUniqueLabels(i) == MElementLabels(j)){
          // set elements in A to the corresponding values
          M(arma::span(MElementRowStart(j),MElementRowEnd(j))) = discreteCINTUniqueValues;
        }
      }
    }

  }

  // setting T0MEANS:
  // check for T0MEANS in MElementLabels:
  Rcpp::String T0MEANSLabel = "T0MEANS";
  for(int j = 0; j < MElementLabels.size(); j++){
    if(T0MEANSLabel == MElementLabels(j)){
      // set elements in A to the corresponding values
      M(arma::span(MElementRowStart(j),MElementRowEnd(j))) = T0MEANS;
    }
  }


  // setting MANIFESTMEANS:
  // check for MANIFESTMEANS in MElementLabels:
  Rcpp::String MANIFESTMEANSLabel = "MANIFESTMEANS";
  for(int j = 0; j < MElementLabels.size(); j++){
    if(MANIFESTMEANSLabel == MElementLabels(j)){
      // set elements in A to the corresponding values
      M(arma::span(MElementRowStart(j),MElementRowEnd(j))) = MANIFESTMEANS;
    }
  }
  return(M);
}
