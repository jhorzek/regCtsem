#ifndef CPPTSEMRAMMODEL_H
#define CPPTSEMRAMMODEL_H

#include <RcppArmadillo.h>
#include "cpptsemmodel.h"

// [[Rcpp :: depends ( RcppArmadillo )]]

// cpptsemRAMmodel class

class cpptsemRAMmodel: public cpptsemmodel{
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

  // constrcutor
  cpptsemRAMmodel(std::string mName,
                  Rcpp::List mCtMatrixList,
                  Rcpp::DataFrame mParameterTable,
                  bool mStationaryT0VAR,
                  bool mStationaryT0MEANS);

  // setter
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
  Rcpp::NumericVector approxRAMGradients(double epsilon);

};

#endif
