#ifndef computeRAMExpectations_H
#define computeRAMExpectations_H

#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat computeRAMExpectedCovariance(arma::mat F, arma::mat A, arma::mat S);
arma::mat computeRAMExpectedMeans(arma::mat F, arma::mat A, arma::colvec M);

#endif
