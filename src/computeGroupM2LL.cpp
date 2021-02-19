#include <RcppArmadillo.h>
#include "computeRAMM2LL.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the -2log likelihood in a RAM model for a  group of people with identical missing structure given the sample size, number of non-missing variables (nObservedVariables),
// observedMeans within this sample, observedCov within this sample, filtered expectedMeans, and filtered expectedCovariance.
using namespace Rcpp;
// [[Rcpp::export]]
double computeGroupM2LL(int sampleSize, int nObservedVariables, arma::colvec observedMeans, arma::mat observedCov,
                        arma::colvec expectedMeans, arma::mat expectedCovariance){
  double m2LL;
  double Nklog2pi = sampleSize*nObservedVariables*std::log(2*M_PI);
  double NlogDetExpCov = sampleSize*std::log(arma::det(expectedCovariance));
  double NtrSSigma = sampleSize*arma::trace(observedCov * arma::inv(expectedCovariance));
  arma::mat Ndist = sampleSize*arma::trans(observedMeans - expectedMeans)*arma::inv(expectedCovariance)*(observedMeans - expectedMeans);
  m2LL = Nklog2pi +
    NlogDetExpCov +
    NtrSSigma+
    Ndist(0,0); // note: dist is a 1x1 matrix; extraction is necessary for the data type to be compatible
  return(m2LL);
}
