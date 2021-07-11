#ifndef COMPUTEKALMANFUNCTIONS_H
#define COMPUTEKALMANFUNCTIONS_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::colvec kalmanFit(const int sampleSize,
                       const int Tpoints,
                       const int nlatent,
                       int nmanifest,
                       const arma::mat kalmanData,
                       arma::mat &latentScores,
                       arma::mat &predictedManifestValues,
                       const Rcpp::List& discreteTimeParameterNames,
                       const arma::colvec& T0MEANSValues,
                       const arma::mat& T0VARValues,
                       const Rcpp::List& discreteDRIFTUnique,
                       const Rcpp::List& discreteCINTUnique,
                       const Rcpp::List& discreteTRAITUnique,
                       const Rcpp::List& discreteDIFFUSIONUnique,
                       const arma::mat& LAMBDAValues,
                       const arma::colvec& MANIFESTMEANSValues,
                       const arma::mat& MANIFESTVARValues);

arma::colvec computeKalmanLatentStatePrediction(arma::mat discreteDRIFTValues,
                                                arma::colvec previousStates,
                                                arma::colvec discreteCINTValues);

arma::mat computeKalmanLatentCovariancePrediction(arma::mat discreteDRIFTValues,
                                                  arma::mat previousCovariances,
                                                  arma::mat discreteDIFFUSIONValues);

arma::colvec computeKalmanManifestPrediction(arma::mat LAMBDAValues,
                                             arma::colvec predictedStates,
                                             arma::colvec MANIFESTMEANSValues);

arma::mat computeKalmanManifestCovariancePrediction(arma::mat LAMBDA,
                                                       arma::mat predictedLatentCovariances,
                                                       arma::mat currentMANIFESTVAR);

#endif
