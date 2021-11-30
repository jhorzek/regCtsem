#ifndef CPPTSEMKALMANMODEL_H
#define CPPTSEMKALMANMODEL_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

// cpptsemKalmanModel class

class cpptsemKalmanModel{
  // inherits from cpptsemmodel
public:
  int sampleSize, Tpoints, nlatent, nmanifest;
  Rcpp::List discreteTimeParameterNames;
  arma::mat kalmanData, latentScores, predictedManifestValues;
  std::string name;
  Rcpp::List ctMatrixList;
  Rcpp::List discreteDRIFTUnique;
  Rcpp::List discreteTRAITUnique;
  Rcpp::List DRIFTHASHExponentialUnique;
  Rcpp::List discreteDIFFUSIONUnique;
  Rcpp::List discreteCINTUnique;
  Rcpp::DataFrame parameterTable;
  arma::mat DRIFTValues;
  arma::mat DIFFUSIONValues;
  arma::mat DIFFUSIONbaseValues;
  arma::mat T0VARValues;
  arma::colvec T0MEANSValues, indM2LL, persons, group;
  arma::mat TRAITVARValues;
  arma::mat MANIFESTVARValues;
  arma::mat LAMBDAValues;
  arma::colvec MANIFESTMEANSValues;
  arma::mat asymptoticDIFFUSION;
  bool update = true,
    hasDiscreteDRIFTUnique = false,
    hasDiscreteTRAITUnique = false,
    hasDRIFTHASHExponentialUnique = false,
    hasDiscreteDIFFUSIONUnique = false,
    hasDiscreteCINTUnique = false,
    computeWasCalled = false,
    stationaryT0VAR = false,
    stationaryT0MEANS = false,
    hasDefinitionVariables = false;

  double m2LL;

  // constructor
  cpptsemKalmanModel(std::string mName,
                     Rcpp::List mCtMatrixList,
                     Rcpp::DataFrame mParameterTable,
                     bool mStationaryT0VAR,
                     bool mStationaryT0MEANS);

  // setter
  void setUpdate(bool mUpdate);
  void setParameterValues(const Rcpp::NumericVector& mParameters, Rcpp::StringVector& parameterLabels);
  void setKalmanMatrixValues(const int selectedGroup);
  void setDiscreteDRIFTUnique(Rcpp::List mDiscreteDRIFTUnique);
  void setDiscreteTRAITUnique(Rcpp::List mDiscreteTRAITUnique);
  void setDRIFTHASHExponentialUnique(Rcpp::List mDRIFTHASHExponentialUnique);
  void setDiscreteDIFFUSIONUnique(Rcpp::List mDiscreteDIFFUSIONUnique);
  void setDiscreteCINTUnique(Rcpp::List mDiscreteCINTUnique);

  void setKalmanData(const Rcpp::List mKalmanData, const bool mhasDefinitionVariables);
  void setKalmanGroupings(arma::colvec mPersons, arma::colvec mGroup);
  void setDiscreteTimeParameterNames(const Rcpp::List mDiscreteTimeParameterNames);
  void setKalmanMatrices(const Rcpp::List& mKalmanMatrices);

  // computation
  void computeAndFitKalman();
  Rcpp::NumericVector approxKalmanGradients(const double epsilon = .000001);

  Rcpp::NumericVector approxKalmanGradient(const double epsilon, const Rcpp::String parName);

  // getter
  Rcpp::NumericVector getParameterValues();

};

#endif
