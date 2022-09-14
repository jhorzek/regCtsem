#include <RcppArmadillo.h>
#include "computeKalmanFunctions.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeKalmanLatentStatePrediction
//' Computes the prediction step for the latent states when using the Kalman filter
//' The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
//' @param discreteDRIFTValues discrete drift values
//' @param previousStates states of previous time point
//' @param discreteCINTValues discrete continuous-time-intercept values
//' @returns vector with predicted states
//' @keywords internal
// [[Rcpp::export]]
arma::colvec computeKalmanLatentStatePrediction(const arma::mat& discreteDRIFTValues,
                                                const arma::colvec& previousStates,
                                                const arma::colvec& discreteCINTValues){
  arma::colvec predictedStates = discreteDRIFTValues*previousStates + discreteCINTValues;

  return(predictedStates);
}
