#ifndef CPPTSEMKALMANMODEL_H
#define CPPTSEMKALMANMODEL_H

#include <RcppArmadillo.h>
#include "cpptsemmodel.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

// cpptsemKalmanModel class

class cpptsemKalmanModel: public cpptsemmodel{
  // inherits from cpptsemmodel
public:
  int sampleSize, Tpoints, nlatent, nmanifest;
  Rcpp::List discreteTimeParameterNames;
  arma::mat kalmanData, latentScores, predictedManifestValues;

  // constrcutor
  cpptsemKalmanModel(std::string mName,
                     Rcpp::List mCtMatrixList,
                     Rcpp::DataFrame mParameterTable,
                     bool mStationaryT0VAR,
                     bool mStationaryT0MEANS);

  // setter
  void setKalmanData(Rcpp::List mKalmanData);
  void setDiscreteTimeParameterNames(Rcpp::List mDiscreteTimeParameterNames);
  void setKalmanMatrices(Rcpp::List mKalmanMatrices);

  // computation
  void computeAndFitKalman();
  Rcpp::NumericVector approxKalmanGradients(double epsilon = .000001);

};

#endif
