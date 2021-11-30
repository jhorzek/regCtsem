#ifndef CPPTSEMRAMMODEL_H
#define CPPTSEMRAMMODEL_H

#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

// cpptsemRAMmodel class

class cpptsemRAMmodel{
  // inherits from cpptsemmodel
  public:

  Rcpp::List RAMdata;
  arma::mat A;
  arma::mat S;
  arma::colvec M;
  arma::mat F;
  Rcpp::DataFrame cppAParameterIndicators;
  Rcpp::DataFrame cppSParameterIndicators;
  Rcpp::DataFrame cppMParameterIndicators;
  arma::mat expectedCovariance;
  arma::mat expectedMeans;

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
  arma::colvec T0MEANSValues;
  arma::mat TRAITVARValues;
  arma::mat MANIFESTVARValues;
  arma::mat LAMBDAValues;
  arma::colvec MANIFESTMEANSValues;
  arma::mat asymptoticDIFFUSION;
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
  cpptsemRAMmodel(std::string mName,
                  Rcpp::List mCtMatrixList,
                  Rcpp::DataFrame mParameterTable,
                  bool mStationaryT0VAR,
                  bool mStationaryT0MEANS);

  // setter
  void setParameterValues(const Rcpp::NumericVector& mParameters, Rcpp::StringVector& parameterLabels);
  void setDiscreteDRIFTUnique(Rcpp::List mDiscreteDRIFTUnique);
  void setDiscreteTRAITUnique(Rcpp::List mDiscreteTRAITUnique);
  void setDRIFTHASHExponentialUnique(Rcpp::List mDRIFTHASHExponentialUnique);
  void setDiscreteDIFFUSIONUnique(Rcpp::List mDiscreteDIFFUSIONUnique);
  void setDiscreteCINTUnique(Rcpp::List mDiscreteCINTUnique);

  void setRAMMatrices(arma::mat mA,
                      arma::mat mS,
                      arma::mat mM,
                      arma::mat mF,
                      Rcpp::DataFrame mCppAParameterIndicators,
                      Rcpp::DataFrame mCppSParameterIndicators,
                      Rcpp::DataFrame mCppMParameterIndicators);

  void setRAMData(Rcpp::List mRAMdata);

  // computation
  void computeRAM();
  void fitRAM();
  Rcpp::NumericVector approxRAMGradients(const double epsilon);
  Rcpp::NumericVector approxRAMGradient(const double epsilon, const Rcpp::String parName);

  // getter
  Rcpp::NumericVector getParameterValues();

};

#endif
