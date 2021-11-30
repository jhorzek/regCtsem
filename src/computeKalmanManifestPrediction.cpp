#include <RcppArmadillo.h>
#include "computeKalmanFunctions.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the Prediction for the manifest variables given the predicted latent states when using the Kalman filter
using namespace Rcpp;
// [[Rcpp::export]]
arma::colvec computeKalmanManifestPrediction(const arma::mat& LAMBDAValues,
                                             const arma::colvec& predictedStates,
                                             const arma::colvec& MANIFESTMEANSValues){
  arma::colvec predictedMANIFEST = LAMBDAValues*predictedStates + MANIFESTMEANSValues;

  return(predictedMANIFEST);
}
