#include <RcppArmadillo.h>
#include "computeRAMExpectations.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the expected means in a RAM model give filter matrix F, directed effects matrix A, and means vector M.
// [[Rcpp::export]]

arma::mat computeRAMExpectedMeans(arma::mat F, arma::mat A, arma::colvec M){

  arma::mat I = arma::eye(arma::size(A));
  arma::mat expectedMeans = F*arma::inv(I-A)*M;

  return(expectedMeans);
}
