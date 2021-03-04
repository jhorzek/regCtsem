#ifndef computeRAMM2LL_H
#define computeRAMM2LL_H

#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

// computes the -2 log likelihood for RAM specified models
double computeRAMM2LL(Rcpp::List RAMdata, arma::colvec expectedMeans, arma::mat expectedCovariance);

// for N = 1
double computeIndividualM2LL(int nObservedVariables, arma::colvec rawData,  arma::colvec expectedMeans, arma::mat expectedCovariance);

// for N > 1
double computeGroupM2LL(int sampleSize, int nObservedVariables, arma::colvec observedMeans, arma::mat observedCov,
                        arma::colvec expectedMeans, arma::mat expectedCovariance);

#endif
