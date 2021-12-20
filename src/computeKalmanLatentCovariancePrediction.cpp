#include <RcppArmadillo.h>
#include "computeKalmanFunctions.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the prediction step for the residual covariance matrix of the latent states when using the Kalman filter
// The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05

using namespace Rcpp;
// [[Rcpp::export]]

arma::mat computeKalmanLatentCovariancePrediction(const arma::mat& discreteDRIFTValues,
                                                  const arma::mat& previousCovariances,
                                                  const arma::mat& discreteDIFFUSIONValues){
  arma::mat predictedLatentCovariances = discreteDRIFTValues * previousCovariances * arma::trans(discreteDRIFTValues) + discreteDIFFUSIONValues;

  return(predictedLatentCovariances);
}
