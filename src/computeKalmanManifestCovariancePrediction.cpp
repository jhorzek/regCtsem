#include <RcppArmadillo.h>
#include "computeKalmanFunctions.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the predicted manifest residual covariance matrix given the latent residual covariance when using the Kalman filter
using namespace Rcpp;
// [[Rcpp::export]]

arma::mat computeKalmanManifestCovariancePrediction(const arma::mat& LAMBDA, const arma::mat& predictedLatentCovariances, const arma::mat& currentMANIFESTVAR){
  arma::mat predictedManifestCovariances = LAMBDA * predictedLatentCovariances * arma::trans(LAMBDA) + currentMANIFESTVAR;

  return(predictedManifestCovariances);
}
