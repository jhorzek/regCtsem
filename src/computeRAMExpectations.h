#ifndef computeRAMExpectations_H
#define computeRAMExpectations_H

#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat computeRAMExpectedCovariance(const arma::mat& F, const arma::mat& A, const arma::mat& S);
arma::mat computeRAMExpectedMeans(const arma::mat& F, const arma::mat& A, const arma::colvec& M);

#endif
