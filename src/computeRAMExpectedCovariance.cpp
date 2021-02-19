#include <RcppArmadillo.h>
#include "computeRAMExpectations.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the expected covariance in a RAM model give filter matrix F, directed effects matrix A, and covariance matrix S.
// [[Rcpp::export]]

arma::mat computeRAMExpectedCovariance(arma::mat F, arma::mat A, arma::mat S){

  arma::mat I = arma::eye(arma::size(A));
  arma::mat expectedCovariance = F*arma::inv(I-A)*S*arma::trans(arma::inv(I-A))*arma::trans(F);

  return(expectedCovariance);
}
