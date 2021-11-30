#ifndef computeRAMM2LL_H
#define computeRAMM2LL_H

#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

// computes the -2 log likelihood for RAM specified models
double computeRAMM2LL(const Rcpp::List& RAMdata, const arma::colvec& expectedMeans, const arma::mat& expectedCovariance);

// for N = 1
double computeIndividualM2LL(const int nObservedVariables, const arma::colvec& rawData, const arma::colvec& expectedMeans, const arma::mat& expectedCovariance);

// for N > 1
double computeGroupM2LL(const int sampleSize, const int nObservedVariables, const arma::colvec& observedMeans, const arma::mat& observedCov,
                        const arma::colvec& expectedMeans, const arma::mat& expectedCovariance);

#endif
