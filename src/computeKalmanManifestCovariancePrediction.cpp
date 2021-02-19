#include <RcppArmadillo.h>
#include "computeKalmanFunctions.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the predicted manifest residual covariance matrix given the latent residual covariance when using the Kalman filter
using namespace Rcpp;
// [[Rcpp::export]]

arma::mat computeKalmanManifestCovariancePrediction(arma::mat LAMBDA, arma::mat predictedLatentCovariances, arma::mat currentMANIFESTVAR){
  arma::mat predictedManifestCovariances = LAMBDA * predictedLatentCovariances * arma::trans(LAMBDA) + currentMANIFESTVAR;

  return(predictedManifestCovariances);
}
