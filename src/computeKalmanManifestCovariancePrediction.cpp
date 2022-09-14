#include <RcppArmadillo.h>
#include "computeKalmanFunctions.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeKalmanManifestCovariancePrediction
//'
//' Computes the predicted manifest residual covariance matrix given the latent residual covariance when using the Kalman filter
//' The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
//' @param LAMBDA loadings matrix
//' @param predictedLatentCovariances predicted latent covariances
//' @param currentMANIFESTVAR current manifest covariance
//' @keywords internal
// [[Rcpp::export]]
arma::mat computeKalmanManifestCovariancePrediction(const arma::mat& LAMBDA,
                                                    const arma::mat& predictedLatentCovariances,
                                                    const arma::mat& currentMANIFESTVAR){
  arma::mat predictedManifestCovariances = LAMBDA * predictedLatentCovariances * arma::trans(LAMBDA) + currentMANIFESTVAR;

  return(predictedManifestCovariances);
}
