#ifndef COMPUTEKALMANFUNCTIONS_H
#define COMPUTEKALMANFUNCTIONS_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

double kalmanFit(int sampleSize,
                            int Tpoints,
                            int nlatent,
                            int nmanifest,
                            arma::mat kalmanData,
                            arma::mat &latentScores,
                            arma::mat &predictedManifestValues,
                            Rcpp::List discreteTimeParameterNames,
                            arma::colvec T0MEANSValues,
                            arma::mat T0VARValues,
                            Rcpp::List discreteDRIFTUnique,
                            Rcpp::List discreteCINTUnique,
                            Rcpp::List discreteTRAITUnique,
                            Rcpp::List discreteDIFFUSIONUnique,
                            arma::mat LAMBDAValues,
                            arma::colvec MANIFESTMEANSValues,
                            arma::mat MANIFESTVARValues);

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
