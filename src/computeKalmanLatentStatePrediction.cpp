#include <RcppArmadillo.h>
#include "computeKalmanFunctions.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the prediction step for the latent states when using the Kalman filter
using namespace Rcpp;
// [[Rcpp::export]]
arma::colvec computeKalmanLatentStatePrediction(arma::mat discreteDRIFTValues,
                                                arma::colvec previousStates,
                                                arma::colvec discreteCINTValues){
  arma::colvec predictedStates = discreteDRIFTValues*previousStates + discreteCINTValues;

  return(predictedStates);
}
