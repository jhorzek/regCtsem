#include <RcppArmadillo.h>
#include "computeKalmanFunctions.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the prediction step for the residual covariance matrix of the latent states when using the Kalman filter
using namespace Rcpp;
// [[Rcpp::export]]

arma::mat computeKalmanLatentCovariancePrediction(arma::mat discreteDRIFTValues,
                                                  arma::mat previousCovariances,
                                                  arma::mat discreteDIFFUSIONValues){
  arma::mat predictedLatentCovariances = discreteDRIFTValues * previousCovariances * arma::trans(discreteDRIFTValues) + discreteDIFFUSIONValues;

  return(predictedLatentCovariances);
}
