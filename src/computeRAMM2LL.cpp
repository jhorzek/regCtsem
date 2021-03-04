#include <RcppArmadillo.h>
#include "computeRAMM2LL.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the discrete time diffusion matrices given the DRIFTHASHInverse, DIFFUSION,
// DRIFTHASHExponentialUnique (list with expm(drifthash)-names, dTs and results), and discreteDIFFUSIONUnique (list with diffusion-names, dTs and results)
using namespace Rcpp;
// [[Rcpp::export]]
double computeRAMM2LL(Rcpp::List RAMdata, arma::colvec expectedMeans, arma::mat expectedCovariance){
  double m2LL = 0, currentM2LL;
  int sampleSize;
  int nObservedVariables;
  Rcpp::List currentMissingnessPattern;
  arma::colvec rawData;
  arma::colvec observedMeans;
  arma::mat  observedCov;
  arma::colvec  currentExpectedMeans;
  arma::mat  currentExpectedCovariance;
  arma::ucolvec  cppExpectedSelector;
  Rcpp::String sampleSizeName = "sampleSize", nObservedVariablesName = "nObservedVariables",
    rawDataName = "rawData", observedMeansName = "observedMeans",
    observedCovName = "observedCov", currentExpectedMeansName = "expectedMeans",
    currentExpectedCovarianceName = "expectedCovariance", cppExpectedSelectorName  = "cppExpectedSelector", currentMissingnessPatternName, m2LLName = "m2LL";

  Rcpp::StringVector missingnessPatternNames = RAMdata["names"];

  for(int i = 0; i < missingnessPatternNames.size(); i++){
    currentMissingnessPatternName = missingnessPatternNames(i);
    currentMissingnessPattern = RAMdata[currentMissingnessPatternName];
    sampleSize = currentMissingnessPattern[sampleSizeName];
    nObservedVariables = currentMissingnessPattern[nObservedVariablesName];
    cppExpectedSelector = Rcpp::as<arma::ucolvec>(currentMissingnessPattern[cppExpectedSelectorName]);

    // Expected Means
    currentExpectedMeans = Rcpp::as<arma::colvec>(currentMissingnessPattern[currentExpectedMeansName]);
    // subset expected means
    currentExpectedMeans = expectedMeans.elem( cppExpectedSelector );

    // Expected Covariance
    currentExpectedCovariance = Rcpp::as<arma::mat>(currentMissingnessPattern[currentExpectedCovarianceName]);
    // subset expected covariance
    currentExpectedCovariance = expectedCovariance.submat( cppExpectedSelector , cppExpectedSelector );

    currentM2LL = currentMissingnessPattern( m2LLName );

    if(sampleSize > 1){
      observedMeans = Rcpp::as<arma::colvec>(currentMissingnessPattern[observedMeansName]);
      observedCov = Rcpp::as<arma::mat>(currentMissingnessPattern[observedCovName]);
      currentM2LL = computeGroupM2LL(sampleSize, nObservedVariables, observedMeans, observedCov,
                                     currentExpectedMeans, currentExpectedCovariance);
      m2LL = m2LL + currentM2LL;
    }else{
      rawData = Rcpp::as<arma::colvec>(currentMissingnessPattern[rawDataName]);
      currentM2LL = computeIndividualM2LL(nObservedVariables, rawData,  currentExpectedMeans, currentExpectedCovariance);
      m2LL = m2LL + currentM2LL;
    }
  }

  return(m2LL);
}
