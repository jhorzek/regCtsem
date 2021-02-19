#ifndef CPPTSEMMODEL_H
#define CPPTSEMMODEL_H

#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

// cpptsemmodel class

class cpptsemmodel{
public:
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
  arma::colvec T0MEANSValues;
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
  cpptsemmodel(std::string mName, Rcpp::List mCtMatrixList, Rcpp::DataFrame mParameterTable,
               bool mStationaryT0VAR, bool mStationaryT0MEANS);

  // setter
  //void setData(arma::mat mDataset, arma::vec mDT);
  void setParameterValues(Rcpp::NumericVector mParameters, Rcpp::StringVector parameterLabels);
  void setDiscreteDRIFTUnique(Rcpp::List mDiscreteDRIFTUnique);
  void setDiscreteTRAITUnique(Rcpp::List mDiscreteTRAITUnique);
  void setDRIFTHASHExponentialUnique(Rcpp::List mDRIFTHASHExponentialUnique);
  void setDiscreteDIFFUSIONUnique(Rcpp::List mDiscreteDIFFUSIONUnique);
  void setDiscreteCINTUnique(Rcpp::List mDiscreteCINTUnique);

  // getter
  Rcpp::NumericVector getParameterValues();

};

#endif
