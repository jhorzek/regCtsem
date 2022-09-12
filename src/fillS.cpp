#include <RcppArmadillo.h>
#include "fillRAMMatrices.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' fillS
//'
//' Fills the S matrix (covariances) given the S matrix, the discreteDIFFUSIONUnique (list with discrete diffusion names, dTs and results),
//' MANIFESTVAR (manifest variance), and cppSParameterIndicators (tells fillS where to put the discreteDIFFUSIONUnique and MANIFESTVAR)
//' The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
//' @param S matrix with undirected effects
//' @param T0VAR matrix with initial covariances
//' @param TRAITVAR trait variance
//' @param hasDiscreteTRAITUnique boolean
//' @param hasDiscreteDIFFUSIONUnique boolean
//' @param discreteDIFFUSIONUnique list with discrete diffusion names, dTs and results
//' @param MANIFESTVAR manifest variance
//' @param SParameterIndicators tells fillS where to put the discreteDIFFUSIONUnique and MANIFESTVAR
//'
// [[Rcpp::export]]

arma::mat  fillS(arma::mat S,
                 const arma::mat& T0VAR,
                 const arma::mat& MANIFESTVAR,
                 const bool hasDiscreteTRAITUnique,
                 const arma::mat& TRAITVAR,
                 const bool hasDiscreteDIFFUSIONUnique,
                 const Rcpp::List& discreteDIFFUSIONUnique,
                 const Rcpp::DataFrame& SParameterIndicators){
  // extract labels of elements in S
  Rcpp::StringVector  SElementLabels = SParameterIndicators["label"];
  arma::uvec  SElementRowStart = SParameterIndicators["row_start"];
  arma::uvec  SElementRowEnd = SParameterIndicators["row_end"];
  arma::uvec  SElementColStart = SParameterIndicators["col_start"];
  arma::uvec  SElementColEnd = SParameterIndicators["col_end"];

  // setting T0VAR
  Rcpp::String T0VARLabel = "T0VAR";
  for(int j = 0; j < SElementLabels.size(); j++){

    if(T0VARLabel == SElementLabels(j)){
      // set elements in A to the corresponding values
      S(arma::span(SElementRowStart(j),SElementRowEnd(j)), arma::span(SElementColStart(j), SElementColEnd(j)) ) = T0VAR;
    }
  }

  if(hasDiscreteDIFFUSIONUnique){
    // setting discrete diffusion:
    Rcpp::StringVector discreteDIFFUSIONUniqueLabels = discreteDIFFUSIONUnique["labels"];
    for(int i = 0; i < discreteDIFFUSIONUniqueLabels.size(); i++){
      Rcpp::String discreteDIFFUSIONUniqueLabel = discreteDIFFUSIONUniqueLabels(i);
      arma::mat discreteDIFFUSIONUniqueValues = discreteDIFFUSIONUnique[discreteDIFFUSIONUniqueLabel];

      // check for this discrete drift in SElementLabels:
      for(int j = 0; j < SElementLabels.size(); j++){
        Rcpp::String SElementLabel = SElementLabels(j);

        if(discreteDIFFUSIONUniqueLabel == SElementLabel){
          // set elements in S to the corresponding values
          S(arma::span(SElementRowStart(j),SElementRowEnd(j)), arma::span(SElementColStart(j), SElementColEnd(j)) ) = discreteDIFFUSIONUniqueValues;
        }
      }
    }
  }

  if(hasDiscreteTRAITUnique){
    // setting TRAITVAR
    Rcpp::String TRAITVARLabel = "TRAITVAR";
    for(int j = 0; j < SElementLabels.size(); j++){
      if(TRAITVARLabel == SElementLabels(j)){
        // set elements in A to the corresponding values
        S(arma::span(SElementRowStart(j),SElementRowEnd(j)), arma::span(SElementColStart(j), SElementColEnd(j)) ) = TRAITVAR;
      }
    }
  }

  // setting MANIFESTVAR:

  // check for MANIFESTVAR in SElementLabels:
  Rcpp::String MANIFESTVARLabel = "MANIFESTVAR";
  for(int j = 0; j < SElementLabels.size(); j++){
    if(MANIFESTVARLabel == SElementLabels(j)){
      // set elements in A to the corresponding values
      S(arma::span(SElementRowStart(j),SElementRowEnd(j)), arma::span(SElementColStart(j), SElementColEnd(j)) ) = MANIFESTVAR;
    }
  }
  return(S);
}
