#include <RcppArmadillo.h>
#include "computeKalmanFunctions.h"
#include "computeRAMM2LL.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Performs the prediciton and updating step for the Kalman Filter
using namespace Rcpp;
// [[Rcpp::export]]
arma::colvec kalmanFit(int sampleSize,
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
                 arma::mat MANIFESTVARValues) {

  double m2LL = 0; // -2 log likelihood
  double currentM2LL = 0;
  arma::colvec indM2LL(sampleSize, arma::fill::zeros); // individual -2 log likelihood
  arma::uvec nonmissing, missing;
  arma::mat LAMBDANotNA, MANIFESTVARNotNA, KalmanGain, predictedLatentVariance, predictedLatentVariance_tMinus1, predictedMANIFESTVariance, filteredCovariance,
  currentDiscreteDRIFTValues, currentDiscreteCINTValues, currentDiscreteDIFFUSIONValues;
  arma::colvec observed, observedNotNA, residual, currentLatentScores, predictedManifest, latentScore_tMinus1, filteredData, filteredMeans;

  // The following String vectors contain the name for each discrete time parameter
  // (e.g.  "discreteDRIFT_1" "discreteDRIFT_1" "discreteDRIFT_2" "discreteDRIFT_2")
  Rcpp::StringMatrix discreteDRIFTNames = discreteTimeParameterNames["discreteDRIFTNames"];
  Rcpp::StringMatrix discreteDIFFUSIONNames = discreteTimeParameterNames["discreteDIFFUSIONNames"];
  Rcpp::StringMatrix discreteCINTNames = discreteTimeParameterNames["discreteCINTNames"];
  Rcpp::String discreteDRIFTName, discreteDIFFUSIONName, discreteCINTName;


  for(int person = 0; person < sampleSize; person++){
    // extract data as colvec
    observed = arma::trans(kalmanData.submat(person, 0,
                                             person, nmanifest-1));
    // Step 1: Prediction
    // set T0MEANS
    currentLatentScores = T0MEANSValues;
    // set T0VAR
    predictedLatentVariance = T0VARValues;
    // Manifest Variables
    predictedManifest = computeKalmanManifestPrediction(LAMBDAValues,
                                                        currentLatentScores,
                                                        MANIFESTMEANSValues);

    // find non-missing
    nonmissing = arma::find_finite(observed);
    missing = arma::find_nonfinite(observed);
    int nObservedVariables = nonmissing.size();

    observedNotNA = observed(nonmissing);
    LAMBDANotNA = LAMBDAValues.rows(nonmissing);
    MANIFESTVARNotNA = MANIFESTVARValues.submat(nonmissing, nonmissing);

    // Step 2: Updating

    // update predictions:
    if(nObservedVariables > 0){
      // Manifest Residual Covariances
      predictedMANIFESTVariance = computeKalmanManifestCovariancePrediction(LAMBDANotNA,
                                                                            predictedLatentVariance,
                                                                            MANIFESTVARNotNA);

      // Kalman Gain
      KalmanGain = predictedLatentVariance * arma::trans(LAMBDANotNA) * arma::inv(predictedMANIFESTVariance);

      // compute residuals: observed-predicted
      residual = observed - predictedManifest;
      // update predicted latent scores and latent variance
      currentLatentScores = currentLatentScores + KalmanGain*(residual(nonmissing));
      predictedLatentVariance = predictedLatentVariance - KalmanGain * LAMBDANotNA * predictedLatentVariance;

      // compute Likelihood --nonmissing
      filteredData = observed(nonmissing);
      filteredMeans = predictedManifest(nonmissing);
      currentM2LL = computeIndividualM2LL(nObservedVariables, filteredData, filteredMeans, predictedMANIFESTVariance);
      indM2LL(person) += currentM2LL;
      m2LL += currentM2LL;
    }
    // save latent scores
    latentScores.submat(person, 0,
                        person, nlatent-1) =
                          arma::trans(currentLatentScores);
    predictedManifestValues.submat(person, 0,
                                   person, nmanifest-1) = arma::trans(predictedManifest);

    // Iterate over all following time points
    for(int Tpoint = 1; Tpoint < Tpoints; Tpoint++){
      //Step 0: prepare

      predictedLatentVariance_tMinus1 = predictedLatentVariance;
      latentScore_tMinus1 = currentLatentScores;

      // extract discrete time parameters for the given time interval
      discreteDRIFTName = discreteDRIFTNames(person, Tpoint-1);
      discreteDIFFUSIONName = discreteDIFFUSIONNames(person, Tpoint-1);
      discreteCINTName = discreteCINTNames(person, Tpoint-1);

      currentDiscreteDRIFTValues = Rcpp::as<arma::mat>(discreteDRIFTUnique[discreteDRIFTName]);
      currentDiscreteDIFFUSIONValues = Rcpp::as<arma::mat>(discreteDIFFUSIONUnique[discreteDIFFUSIONName]);
      currentDiscreteCINTValues = Rcpp::as<arma::mat>(discreteCINTUnique[discreteCINTName]);

      // Extract observed data for person as colvec
      observed = arma::trans(kalmanData.submat(person, Tpoint*nmanifest,
                                               person, (Tpoint+1) * nmanifest-1));

      // find non-missing data
      nonmissing = arma::find_finite(observed);
      missing = arma::find_nonfinite(observed);
      nObservedVariables = nonmissing.size();

      // Step 1: Prediction

      // States:
      currentLatentScores = computeKalmanLatentStatePrediction(currentDiscreteDRIFTValues,
                                                               latentScore_tMinus1,
                                                               currentDiscreteCINTValues);

      // Latent Residual Covariances
      predictedLatentVariance = computeKalmanLatentCovariancePrediction(currentDiscreteDRIFTValues,
                                                                        predictedLatentVariance_tMinus1,
                                                                        currentDiscreteDIFFUSIONValues);

      // Manifest Variables
      predictedManifest = computeKalmanManifestPrediction(LAMBDAValues,
                                                          currentLatentScores,
                                                          MANIFESTMEANSValues);

      // Step 2: Updating

      observedNotNA = observed(nonmissing);
      LAMBDANotNA = LAMBDAValues.rows(nonmissing);
      MANIFESTVARNotNA = MANIFESTVARValues.submat(nonmissing, nonmissing);

      // update predictions:
      if(nObservedVariables > 0){
        // Manifest Residual Covariances
        predictedMANIFESTVariance = computeKalmanManifestCovariancePrediction(LAMBDANotNA,
                                                                              predictedLatentVariance,
                                                                              MANIFESTVARNotNA);

        // Kalman Gain
        KalmanGain = predictedLatentVariance * arma::trans(LAMBDANotNA) * arma::inv(predictedMANIFESTVariance);

        // compute residuals
        residual = observed - predictedManifest;
        // update scores and latent variance
        currentLatentScores = currentLatentScores + KalmanGain*(residual(nonmissing));
        predictedLatentVariance = predictedLatentVariance - KalmanGain * LAMBDANotNA * predictedLatentVariance;

        // compute Likelihood --nonmissing
        filteredData = observed(nonmissing);
        filteredMeans = predictedManifest(nonmissing);
        currentM2LL = computeIndividualM2LL(nObservedVariables, filteredData, filteredMeans, predictedMANIFESTVariance);
        indM2LL(person) += currentM2LL;
        m2LL += currentM2LL;
      }

      // save predicted scores
      latentScores.submat(person, Tpoint*nlatent,
                          person, (Tpoint+1) * nlatent-1) = arma::trans(
                            currentLatentScores);

      // save predicted manifests
      predictedManifestValues.submat(person, Tpoint*nmanifest,
                                     person, (Tpoint+1) * nmanifest-1) = arma::trans(predictedManifest);

    }

  }


  return(indM2LL);
}


