#include <RcppArmadillo.h>
#include "computeRAMM2LL.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeIndividualM2LL
//'
//' Computes the -2log likelihood in a RAM model for a  single person given the number of non-missing variables (nObservedVariables), x (rawData), filtered expectedMeans, and filtered expectedCovariance.
//' @param nObservedVariables number of non-missing variables
//' @param rawData x
//' @param expectedMeans filtered expectedMeans
//' @param expectedCovariance filtered expectedCovariance
//' @keywords internal
// [[Rcpp::export]]
double computeIndividualM2LL(const int nObservedVariables, const arma::colvec& rawData, const arma::colvec& expectedMeans, const arma::mat& expectedCovariance){
  double m2LL;
  double klog2pi = nObservedVariables*std::log(2*M_PI);
  double logDetExpCov = std::log(arma::det(expectedCovariance));
  arma::mat dist = arma::trans(rawData - expectedMeans)*arma::inv(expectedCovariance)*(rawData - expectedMeans);
  m2LL = klog2pi +
    logDetExpCov +
    dist(0,0); // note: dist is a 1x1 matrix; extraction is necessary for the data type to be compatible
  return(m2LL);
}
