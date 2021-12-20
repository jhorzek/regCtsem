#include <RcppArmadillo.h>
#include "computeRAMExpectations.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the expected means in a RAM model give filter matrix F, directed effects matrix A, and means vector M.
// The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05

// [[Rcpp::export]]

arma::mat computeRAMExpectedMeans(const arma::mat& F, const arma::mat& A, const arma::colvec& M){

  arma::mat I = arma::eye(arma::size(A));
  arma::mat expectedMeans = F*arma::inv(I-A)*M;

  return(expectedMeans);
}
