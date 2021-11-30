#ifndef COMPUTEKALMANFUNCTIONS_H
#define COMPUTEKALMANFUNCTIONS_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

arma::colvec kalmanFit(
    bool update,
    const int sampleSize,
                       const int Tpoints,
                       const int nlatent,
                       const int nmanifest,
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

arma::colvec computeKalmanLatentStatePrediction(const arma::mat& discreteDRIFTValues,
                                                const arma::colvec& previousStates,
                                                const arma::colvec& discreteCINTValues);

arma::mat computeKalmanLatentCovariancePrediction(const arma::mat& discreteDRIFTValues,
                                                  const arma::mat& previousCovariances,
                                                  const arma::mat& discreteDIFFUSIONValues);

arma::colvec computeKalmanManifestPrediction(const arma::mat& LAMBDAValues,
                                             const arma::colvec& predictedStates,
                                             const arma::colvec& MANIFESTMEANSValues);

arma::mat computeKalmanManifestCovariancePrediction(const arma::mat& LAMBDA,
                                                    const arma::mat& predictedLatentCovariances,
                                                    const arma::mat& currentMANIFESTVAR);

#endif
