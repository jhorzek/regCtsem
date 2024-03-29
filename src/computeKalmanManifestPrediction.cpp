#include <RcppArmadillo.h>
#include "computeKalmanFunctions.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeKalmanManifestPrediction
//'
//' Computes the Prediction for the manifest variables given the predicted latent states when using the Kalman filter
//' The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
//' @param LAMBDAValues matrix with loadings
//' @param predictedStates vector with predicted latent states
//' @param MANIFESTMEANSValues vector with means of manifest variables
//' @returns predicted manifest values
//' @keywords internal
// [[Rcpp::export]]
arma::colvec computeKalmanManifestPrediction(const arma::mat& LAMBDAValues,
                                             const arma::colvec& predictedStates,
                                             const arma::colvec& MANIFESTMEANSValues){
  arma::colvec predictedMANIFEST = LAMBDAValues*predictedStates + MANIFESTMEANSValues;

  return(predictedMANIFEST);
}
