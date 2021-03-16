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
  arma::mat T0VARValues;
  arma::colvec T0MEANSValues, indM2LL;
  arma::mat TRAITVARValues;
  arma::mat MANIFESTVARValues;
  arma::mat LAMBDAValues;
  arma::colvec MANIFESTMEANSValues;
  bool hasDiscreteDRIFTUnique = false,
    hasDiscreteTRAITUnique = false,
    hasDRIFTHASHExponentialUnique = false,
    hasDiscreteDIFFUSIONUnique = false,
    hasDiscreteCINTUnique = false,
    computeWasCalled = false,
    stationaryT0VAR = false,
    stationaryT0MEANS = false;

  double m2LL;

  // constructor
  cpptsemKalmanModel(std::string mName,
                     Rcpp::List mCtMatrixList,
                     Rcpp::DataFrame mParameterTable,
                     bool mStationaryT0VAR,
                     bool mStationaryT0MEANS);

  // setter
  void setParameterValues(Rcpp::NumericVector mParameters, Rcpp::StringVector parameterLabels);
  void setDiscreteDRIFTUnique(Rcpp::List mDiscreteDRIFTUnique);
  void setDiscreteTRAITUnique(Rcpp::List mDiscreteTRAITUnique);
  void setDRIFTHASHExponentialUnique(Rcpp::List mDRIFTHASHExponentialUnique);
  void setDiscreteDIFFUSIONUnique(Rcpp::List mDiscreteDIFFUSIONUnique);
  void setDiscreteCINTUnique(Rcpp::List mDiscreteCINTUnique);

  void setKalmanData(Rcpp::List mKalmanData);
  void setDiscreteTimeParameterNames(Rcpp::List mDiscreteTimeParameterNames);
  void setKalmanMatrices(Rcpp::List mKalmanMatrices);

  // computation
  void computeAndFitKalman();
  Rcpp::NumericVector approxKalmanGradients(double epsilon = .000001);

  Rcpp::NumericVector approxKalmanGradient(double epsilon, Rcpp::String parName);

  double computePenalty(Rcpp::NumericVector pars,
                        Rcpp::StringVector regIndicators,
                        Rcpp::NumericVector lambda);

  Rcpp::NumericVector computeSubgradients(Rcpp::NumericVector pars,
                                          Rcpp::NumericVector gradients,
                                          Rcpp::StringVector regIndicators,
                                          Rcpp::NumericVector lambdas);

  Rcpp::List GIST(
      Rcpp::NumericVector pars,
      Rcpp::StringVector regIndicators,
      Rcpp::NumericVector lambda,
      double eta, double sig,
      Rcpp::NumericVector gradient_epsilons,
      double initialStepsize, double stepsizeMin = 0, double stepsizeMax = 999999999,
      std::string GISTLinesearchCriterion = "monotone",
      int GISTNonMonotoneNBack = 5,
      int maxIter_out = 100, int maxIter_in = 100,
      double break_outer = .00000001,
      std::string break_crit = "parameterChange", int verbose = 0);

  // getter
  Rcpp::NumericVector getParameterValues();

};

#endif
