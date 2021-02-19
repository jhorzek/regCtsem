#include "cpptsemmodel.h"
#include "computeDiscreteParameters.h"
#include "getVarianceFromVarianceBase.h"
#include "fillRAMMatrices.h"
#include "computeRAMExpectations.h"
#include "computeRAMM2LL.h"
#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

using namespace Rcpp;

// cpptsemmodel class

// Constructor:

cpptsemmodel::cpptsemmodel(std::string mName, Rcpp::List mCtMatrixList, Rcpp::DataFrame mParameterTable,
                           bool mStationaryT0VAR, bool mStationaryT0MEANS){
  // Define model name and the ct matrices
  name = mName;
  ctMatrixList = mCtMatrixList;

  // Adjust the row and col indicators to start at 0
  Rcpp::NumericVector mParameterTableRow = mParameterTable["row"];
  mParameterTableRow = mParameterTableRow-1;
  Rcpp::NumericVector mParameterTableCol = mParameterTable["col"];
  mParameterTableCol = mParameterTableCol-1;

  parameterTable = mParameterTable;

  stationaryT0VAR = mStationaryT0VAR;
  stationaryT0MEANS = mStationaryT0MEANS;

}

void cpptsemmodel::setDiscreteDRIFTUnique(Rcpp::List mDiscreteDRIFTUnique){
  // list with labels of the unique discrete drifts, corresponding time intervals and results
  discreteDRIFTUnique = mDiscreteDRIFTUnique;
  hasDiscreteDRIFTUnique = true;
}
void cpptsemmodel::setDiscreteTRAITUnique(Rcpp::List mDiscreteTRAITUnique){
  // list with labels of the unique discrete traits, corresponding time intervals and results
  discreteTRAITUnique = mDiscreteTRAITUnique;
  hasDiscreteTRAITUnique = true;
}

void cpptsemmodel::setDRIFTHASHExponentialUnique(Rcpp::List mDRIFTHASHExponentialUnique){
  // list with labels of the unique discrete traits, corresponding time intervals and results
  DRIFTHASHExponentialUnique = mDRIFTHASHExponentialUnique;
  hasDRIFTHASHExponentialUnique = true;
}

void cpptsemmodel::setDiscreteDIFFUSIONUnique(Rcpp::List mDiscreteDIFFUSIONUnique){
  discreteDIFFUSIONUnique = mDiscreteDIFFUSIONUnique;
  hasDiscreteDIFFUSIONUnique = true;
}

void cpptsemmodel::setDiscreteCINTUnique(Rcpp::List mDiscreteCINTUnique){
  discreteCINTUnique = mDiscreteCINTUnique;
  hasDiscreteCINTUnique = true;
}

void cpptsemmodel::setParameterValues(Rcpp::NumericVector mParameters, Rcpp::StringVector parameterLabels){
  // change status of RAM computation and RAM fitting
  computeWasCalled = false;
  // set the parameters in (1) the parameterTable (2) the ctMatrices
  bool wasSet; // checks if all parameters were changed
  Rcpp::StringVector expectedParameterLabels = parameterTable["label"];
  Rcpp::StringVector matrixLabels = parameterTable["matrix"]; // vector with labels of the ct matrices
  Rcpp::NumericVector currentParameters = parameterTable["value"]; // values before change
  Rcpp::NumericVector row =  parameterTable["row"];
  Rcpp::NumericVector col =  parameterTable["col"];

  for(int i = 0; i < expectedParameterLabels.size(); i++){
    wasSet = false;
    for(int j = 0; j < parameterLabels.size(); j++){
      if(parameterLabels(j) == expectedParameterLabels(i)){
        wasSet = true;
        currentParameters(i) = mParameters(j);
      }
    }

    if(!wasSet){
      Rcpp::Rcout << "Error while setting parameters: The parameter " << expectedParameterLabels(i) << " was not found!" << std::endl;
    }
  }

  for(int i = 0; i < expectedParameterLabels.size(); i++){
    Rcpp::String matrixLabel = matrixLabels(i); // matrix in which parameter i is situated
    Rcpp::List currentMatrixList = ctMatrixList[matrixLabel]; // values and labels of this matrix
    Rcpp::NumericMatrix currentMatrixValues = currentMatrixList["values"]; // values of this matrix
    currentMatrixValues(row(i), col(i)) = currentParameters(i); // set the values to the new parameter values
  }

}

// returns the parameter values of a cpptsemmodel
Rcpp::NumericVector cpptsemmodel::getParameterValues() {

  Rcpp::StringVector parameterLabels = parameterTable["label"];
  Rcpp::StringVector uniqueParameterLabels = unique(parameterLabels);
  Rcpp::NumericVector parameterValuesRep = parameterTable["value"]; // values with repeated elements

  Rcpp::NumericVector paramterValues (uniqueParameterLabels.length(),-9999.99);

  for(int i = 0; i < uniqueParameterLabels.size(); i++){

    for(int j = 0; j < parameterLabels.length(); j++)
      if(uniqueParameterLabels(i) == parameterLabels(j)){
        paramterValues(i) = parameterValuesRep(j);
        break;
      }
  }
  paramterValues.names() = uniqueParameterLabels;

  return(Rcpp::clone(paramterValues));
}



RCPP_MODULE(cpptsemmodule){
  Rcpp::class_<cpptsemmodel>( "cpptsemmodel" )

  .constructor<std::string, Rcpp::List, Rcpp::DataFrame, bool, bool>("Creates a new model. Expects a model name, ctMatrixList, parameterTable, boolean indicating stationaryT0VAR, boolean indicating T0MEANS")
  .field_readonly( "name", &cpptsemmodel::name, "Name of the model")
  .field_readonly( "ctMatrixList", &cpptsemmodel::ctMatrixList, "List of ct matrices")
  .field_readonly( "parameterTable", &cpptsemmodel::parameterTable, "Data frame of model parameters")
  .field_readonly( "discreteDRIFTUnique", &cpptsemmodel::discreteDRIFTUnique, "discreteDRIFTUnique")
  .field_readonly( "discreteTRAITUnique", &cpptsemmodel::discreteTRAITUnique, "discreteTRAITUnique")
  .field_readonly( "DRIFTHASHExponentialUnique", &cpptsemmodel::DRIFTHASHExponentialUnique, "DRIFTHASHExponentialUnique")
  .field_readonly( "discreteDIFFUSIONUnique", &cpptsemmodel::discreteDIFFUSIONUnique, "discreteDIFFUSIONUnique")
  .field_readonly( "DRIFTValues", &cpptsemmodel::DRIFTValues, "DRIFTValues")
  .field_readonly( "DIFFUSIONValues", &cpptsemmodel::DIFFUSIONValues, "DIFFUSIONValues")
  .field_readonly( "T0VARValues", &cpptsemmodel::T0VARValues, "T0VARValues")
  .field_readonly( "T0MEANSValues", &cpptsemmodel::T0MEANSValues, "T0MEANSValues")
  .field_readonly( "TRAITVARValues", &cpptsemmodel::TRAITVARValues, "TRAITVARValues")
  .field_readonly( "MANIFESTMEANSValues", &cpptsemmodel::MANIFESTMEANSValues, "MANIFESTMEANSValues")
  .field_readonly( "MANIFESTVARValues", &cpptsemmodel::MANIFESTVARValues, "MANIFESTVARValues")
  .field_readonly( "LAMBDAValues", &cpptsemmodel::LAMBDAValues, "LAMBDAValues")
  .field_readonly( "m2LL", &cpptsemmodel::m2LL, "-2 log likelihood")

  // methods
  .method( "setParameterValues", &cpptsemmodel::setParameterValues, "Set the parameters. Expects a vector with parametervalues and a stringvector with labels")
  .method( "setDiscreteDRIFTUnique", &cpptsemmodel::setDiscreteDRIFTUnique, "Set the number of discrete time drifts, corresponding dts and values")
  .method( "setDiscreteTRAITUnique", &cpptsemmodel::setDiscreteTRAITUnique, "Set the number of discrete time traits, corresponding dts and values")
  .method( "setDRIFTHASHExponentialUnique", &cpptsemmodel::setDRIFTHASHExponentialUnique, "Set the number of drift hash exponentials, corresponding dts and values")
  .method( "setDiscreteDIFFUSIONUnique", &cpptsemmodel::setDiscreteDIFFUSIONUnique, "Set the number of discrete diffusions, corresponding dts and values")
  .method( "setDiscreteCINTUnique", &cpptsemmodel::setDiscreteCINTUnique, "Set the number of discrete continuous time intercepts, corresponding dts and values")
  .method( "getParameterValues", &cpptsemmodel::getParameterValues, "Get current parameter values")
  ;
}


