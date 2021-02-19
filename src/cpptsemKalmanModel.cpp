#include "cpptsemKalmanModel.h"
#include "computeDiscreteParameters.h"
#include "getVarianceFromVarianceBase.h"
#include "computeKalmanFunctions.h"
#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

using namespace Rcpp;

// Constructor
cpptsemKalmanModel::cpptsemKalmanModel(std::string mName,
                                       Rcpp::List mCtMatrixList,
                                       Rcpp::DataFrame mParameterTable,
                                       bool mStationaryT0VAR,
                                       bool mStationaryT0MEANS){
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

// Basic functions

void cpptsemKalmanModel::setDiscreteDRIFTUnique(Rcpp::List mDiscreteDRIFTUnique){
  // list with labels of the unique discrete drifts, corresponding time intervals and results
  discreteDRIFTUnique = mDiscreteDRIFTUnique;
  hasDiscreteDRIFTUnique = true;
}
void cpptsemKalmanModel::setDiscreteTRAITUnique(Rcpp::List mDiscreteTRAITUnique){
  // list with labels of the unique discrete traits, corresponding time intervals and results
  discreteTRAITUnique = mDiscreteTRAITUnique;
  hasDiscreteTRAITUnique = true;
}

void cpptsemKalmanModel::setDRIFTHASHExponentialUnique(Rcpp::List mDRIFTHASHExponentialUnique){
  // list with labels of the unique discrete traits, corresponding time intervals and results
  DRIFTHASHExponentialUnique = mDRIFTHASHExponentialUnique;
  hasDRIFTHASHExponentialUnique = true;
}

void cpptsemKalmanModel::setDiscreteDIFFUSIONUnique(Rcpp::List mDiscreteDIFFUSIONUnique){
  discreteDIFFUSIONUnique = mDiscreteDIFFUSIONUnique;
  hasDiscreteDIFFUSIONUnique = true;
}

void cpptsemKalmanModel::setDiscreteCINTUnique(Rcpp::List mDiscreteCINTUnique){
  discreteCINTUnique = mDiscreteCINTUnique;
  hasDiscreteCINTUnique = true;
}

void cpptsemKalmanModel::setParameterValues(Rcpp::NumericVector mParameters, Rcpp::StringVector parameterLabels){
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

// returns the parameter values of a cpptsemKalmanModel
Rcpp::NumericVector cpptsemKalmanModel::getParameterValues() {

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

// Kalman specific functions

void cpptsemKalmanModel::setKalmanMatrices(Rcpp::List mKalmanMatrices){

  latentScores = Rcpp::as<arma::mat>(mKalmanMatrices["latentScores"]);
  predictedManifestValues = Rcpp::as<arma::mat>(mKalmanMatrices["predictedManifestValues"]);

}

void cpptsemKalmanModel::setKalmanData(Rcpp::List mKalmanData){
  kalmanData = Rcpp::as<arma::mat>(mKalmanData["dataset"]);
  sampleSize = mKalmanData["sampleSize"];
  Tpoints = mKalmanData["Tpoints"];
  nlatent = mKalmanData["nlatent"];
  nmanifest = mKalmanData["nmanifest"];
}

void cpptsemKalmanModel::setDiscreteTimeParameterNames(Rcpp::List mDiscreteTimeParameterNames){
  discreteTimeParameterNames = mDiscreteTimeParameterNames;
}

void cpptsemKalmanModel::computeAndFitKalman(){
  computeWasCalled = true;
  arma::mat DRIFTInverseValues;
  arma::mat DRIFTHASH;
  arma::mat DRIFTHASHInverse;

  if( ctMatrixList.containsElementNamed("DRIFT") ){
    // get Drift
    Rcpp::List DRIFTList = ctMatrixList["DRIFT"];
    DRIFTValues = Rcpp::as<arma::mat> (DRIFTList["values"]);

    // get drift inverse
    DRIFTInverseValues = arma::inv(DRIFTValues);

    // get DRIFTHASH
    DRIFTHASH = computeDRIFTHASH(DRIFTValues);
    DRIFTHASHInverse = arma::inv(DRIFTHASH);

  }else{
    Rcpp::Rcout << "Could not find a DRIFT!" << std::endl;
  }

  if( ctMatrixList.containsElementNamed("DIFFUSIONbase") ){
    // get DIFFUSION
    Rcpp::List DIFFUSIONbaseList = ctMatrixList["DIFFUSIONbase"];
    arma::mat DIFFUSIONbaseValues = DIFFUSIONbaseList["values"];
    DIFFUSIONValues = getVarianceFromVarianceBase(DIFFUSIONbaseValues);
  }else{
    Rcpp::Rcout << "Could not find a DIFFUSION!" << std::endl;
  }

  if( ctMatrixList.containsElementNamed("T0VARbase") ){
    // get T0VAR
    Rcpp::List T0VARbaseList = ctMatrixList["T0VARbase"];
    arma::mat T0VARbaseValues = T0VARbaseList["values"];
    T0VARValues = getVarianceFromVarianceBase(T0VARbaseValues);
  }

  if(stationaryT0VAR){
    // compute asymptotic diffusion
    arma::uword nrows = DIFFUSIONValues.n_rows;
    arma::uword ncols = DIFFUSIONValues.n_cols;
    arma::mat asymptoticDIFFUSION = -1*arma::inv(DRIFTHASH) * arma::vectorise(DIFFUSIONValues);
    T0VARValues = asymptoticDIFFUSION;
    T0VARValues.reshape(nrows, ncols);
  }

  if( ctMatrixList.containsElementNamed("TRAITVARbase") ){
    // get TRAITVAR
    Rcpp::List TRAITVARbaseList = ctMatrixList["TRAITVARbase"];
    arma::mat TRAITVARbaseValues = TRAITVARbaseList["values"];
    TRAITVARValues = getVarianceFromVarianceBase(TRAITVARbaseValues);
  }

  if( ctMatrixList.containsElementNamed("MANIFESTVARbase") ){
    // MANIFESTVAR
    Rcpp::List MANIFESTVARbaseList = ctMatrixList["MANIFESTVARbase"];
    arma::mat MANIFESTVARbaseValues = MANIFESTVARbaseList["values"];
    MANIFESTVARValues = getVarianceFromVarianceBase(MANIFESTVARbaseValues);
  }

  // Manifest Means and T0MEANS
  Rcpp::List MANIFESTMEANSList = ctMatrixList["MANIFESTMEANS"];
  MANIFESTMEANSValues = Rcpp::as<arma::colvec> (MANIFESTMEANSList["values"]);

  if(stationaryT0MEANS){
    Rcpp::List CINTList = ctMatrixList["CINT"];
    arma::colvec CINTValues = Rcpp::as<arma::colvec> (CINTList["values"]);
    T0MEANSValues = -1*DRIFTInverseValues * CINTValues;
  }else{
    Rcpp::List T0MEANSList = ctMatrixList["T0MEANS"];
    T0MEANSValues = Rcpp::as<arma::colvec> (T0MEANSList["values"]);
  }

  Rcpp::List LAMBDAList = ctMatrixList["LAMBDA"];
  LAMBDAValues = Rcpp::as<arma::mat> (LAMBDAList["values"]);

  // compute discrete time drift
  if(hasDiscreteDRIFTUnique){
    discreteDRIFTUnique = computeDiscreteDRIFTs(DRIFTValues, discreteDRIFTUnique);
  }else{
    Rcpp::Rcout << "discreteDRIFTUnique seems to be missing!" << std::endl;
  }

  // compute discrete time traits
  if(hasDiscreteTRAITUnique){
    discreteTRAITUnique = computeDiscreteTRAITs(discreteDRIFTUnique, discreteTRAITUnique);
  }

  // compute discrete DRIFTHASHs
  if(hasDRIFTHASHExponentialUnique){
    DRIFTHASHExponentialUnique = computeDRIFTHASHExponentials(DRIFTHASH, DRIFTHASHExponentialUnique);
  }else{
    Rcpp::Rcout << "DRIFTHASHExponentialUnique seems to be missing!" << std::endl;
  }

  // compute discreteDIFFUSION
  if(hasDiscreteDIFFUSIONUnique){
    if(!hasDiscreteDIFFUSIONUnique){
      Rcpp::Rcout << "Error: DRIFTHASHExponentialUnique missing!" << std::endl;
    }
    discreteDIFFUSIONUnique = computeDiscreteDIFFUSIONs(DRIFTHASHInverse, DIFFUSIONValues,
                                                        DRIFTHASHExponentialUnique, discreteDIFFUSIONUnique);
  }

  // compute discreteCINT
  if(hasDiscreteCINTUnique){
    if(!hasDiscreteDRIFTUnique){
      Rcpp::Rcout << "Error: hasDiscreteDRIFTUnique missing!" << std::endl;
    }
    Rcpp::List CINTList = ctMatrixList["CINT"];
    arma::colvec CINTValues = Rcpp::as<arma::colvec>(CINTList["values"]);

    discreteCINTUnique = computeDiscreteCINTs(discreteCINTUnique, DRIFTInverseValues, discreteDRIFTUnique, CINTValues);
  }


  // Kalman filter: Prediction and Updating

  m2LL = kalmanFit(sampleSize,
                   Tpoints,
                   nlatent,
                   nmanifest,
                   kalmanData,
                   latentScores,
                   predictedManifestValues,
                   discreteTimeParameterNames,
                   T0MEANSValues,
                   T0VARValues,
                   discreteDRIFTUnique,
                   discreteCINTUnique,
                   discreteTRAITUnique,
                   discreteDIFFUSIONUnique,
                   LAMBDAValues,
                   MANIFESTMEANSValues,
                   MANIFESTVARValues);

}

// compute gradient vector. Returns the central gradient approximation
Rcpp::NumericVector cpptsemKalmanModel::approxKalmanGradients(double epsilon){
  Rcpp::NumericVector parameterValues = getParameterValues();
  Rcpp::StringVector parameterNames = parameterValues.names();
  Rcpp::NumericMatrix likelihoods( parameterValues.length() , 2 );
  Rcpp::NumericVector gradients( parameterValues.length());


  rownames(likelihoods) = parameterNames;
  Rcpp::StringVector colNames = {"stepBackward", "stepForward"};
  colnames(likelihoods) = colNames;
  gradients.names() = parameterNames;

  // avoid changing the parameterValues by reference:
  Rcpp::NumericVector currentParameterValues = Rcpp::clone( parameterValues );
  for(int parameter = 0; parameter < parameterNames.length(); parameter++){
    // forward step
    currentParameterValues(parameter) = currentParameterValues(parameter) + epsilon;
    setParameterValues(currentParameterValues, parameterNames);
    computeAndFitKalman();
    likelihoods(parameter,0) = m2LL;

    // reset parameter
    currentParameterValues(parameter) = currentParameterValues(parameter) - epsilon;
    setParameterValues(currentParameterValues, parameterNames);

    // backward step
    currentParameterValues(parameter) = currentParameterValues(parameter) - epsilon;
    setParameterValues(currentParameterValues, parameterNames);
    computeAndFitKalman();
    likelihoods(parameter,1) = m2LL;

    // reset parameter
    currentParameterValues(parameter) = currentParameterValues(parameter) + epsilon;
    setParameterValues(currentParameterValues, parameterNames);
    gradients(parameter) = (likelihoods(parameter,0) - likelihoods(parameter,1))/(2*epsilon);
  }

  return(Rcpp::clone(gradients));
}

RCPP_EXPOSED_CLASS(cpptsemKalmanModel)
RCPP_MODULE(cpptsemKalmanModel_cpp){
  using namespace Rcpp;
  Rcpp::class_<cpptsemKalmanModel>( "cpptsemKalmanModel" )
  .constructor<std::string, Rcpp::List, Rcpp::DataFrame, bool, bool>("Creates a new model. Expects a model name, ctMatrixList, parameterTable")
  .field_readonly( "name", &cpptsemKalmanModel::name, "Name of the model")
  .field_readonly( "ctMatrixList", &cpptsemKalmanModel::ctMatrixList, "List of ct matrices")
  .field_readonly( "parameterTable", &cpptsemKalmanModel::parameterTable, "Data frame of model parameters")
  .field_readonly( "discreteDRIFTUnique", &cpptsemKalmanModel::discreteDRIFTUnique, "discreteDRIFTUnique")
  .field_readonly( "discreteTRAITUnique", &cpptsemKalmanModel::discreteTRAITUnique, "discreteTRAITUnique")
  .field_readonly( "discreteCINTUnique", &cpptsemKalmanModel::discreteCINTUnique, "discreteCINTUnique")
  .field_readonly( "DRIFTHASHExponentialUnique", &cpptsemKalmanModel::DRIFTHASHExponentialUnique, "DRIFTHASHExponentialUnique")
  .field_readonly( "discreteDIFFUSIONUnique", &cpptsemKalmanModel::discreteDIFFUSIONUnique, "discreteDIFFUSIONUnique")
  .field_readonly( "DRIFTValues", &cpptsemKalmanModel::DRIFTValues, "DRIFTValues")
  .field_readonly( "DIFFUSIONValues", &cpptsemKalmanModel::DIFFUSIONValues, "DIFFUSIONValues")
  .field_readonly( "T0VARValues", &cpptsemKalmanModel::T0VARValues, "T0VARValues")
  .field_readonly( "T0MEANSValues", &cpptsemKalmanModel::T0MEANSValues, "T0MEANSValues")
  .field_readonly( "TRAITVARValues", &cpptsemKalmanModel::TRAITVARValues, "TRAITVARValues")
  .field_readonly( "MANIFESTMEANSValues", &cpptsemKalmanModel::MANIFESTMEANSValues, "MANIFESTMEANSValues")
  .field_readonly( "MANIFESTVARValues", &cpptsemKalmanModel::MANIFESTVARValues, "MANIFESTVARValues")
  .field_readonly( "LAMBDAValues", &cpptsemKalmanModel::LAMBDAValues, "LAMBDAValues")
  .field_readonly( "m2LL", &cpptsemKalmanModel::m2LL, "-2 log likelihood")
  .field_readonly( "kalmanData", &cpptsemKalmanModel::kalmanData, "Data for Kalman model")
  .field_readonly( "latentScores", &cpptsemKalmanModel::latentScores, "latentScores")
  .field_readonly( "predictedManifestValues", &cpptsemKalmanModel::predictedManifestValues, "predictedManifestValues")


  // methods
  .method( "setParameterValues", &cpptsemKalmanModel::setParameterValues, "Set the parameters. Expects a vector with parametervalues and a stringvector with labels")
  .method( "setDiscreteDRIFTUnique", &cpptsemKalmanModel::setDiscreteDRIFTUnique, "Set the number of discrete time drifts, corresponding dts and values")
  .method( "setDiscreteTRAITUnique", &cpptsemKalmanModel::setDiscreteTRAITUnique, "Set the number of discrete time traits, corresponding dts and values")
  .method( "setDRIFTHASHExponentialUnique", &cpptsemKalmanModel::setDRIFTHASHExponentialUnique, "Set the number of drift hash exponentials, corresponding dts and values")
  .method( "setDiscreteDIFFUSIONUnique", &cpptsemKalmanModel::setDiscreteDIFFUSIONUnique, "Set the number of discrete diffusions, corresponding dts and values")
  .method( "setDiscreteCINTUnique", &cpptsemKalmanModel::setDiscreteCINTUnique, "Set the number of discrete continuous time intercepts, corresponding dts and values")
  .method( "getParameterValues", &cpptsemKalmanModel::getParameterValues, "Get current parameter values")
  .method( "setKalmanMatrices", &cpptsemKalmanModel::setKalmanMatrices, "Set up the Kalman matrices. Requires numeric matrices for A, S, M, F and DataFrames with cpp compatible row and column indicators for A, S, M which specify where in the matrices the discrete time parameters go.")
  .method( "setKalmanData", &cpptsemKalmanModel::setKalmanData, "Set up the dataset for Kalman models. Expects a list with sampleSize, nObservedVariables, rawData, observedMeans, observedCov, expectedSelector, expectedMeans, expectedCovariance, m2LL for each missingness pattern")
  .method( "setDiscreteTimeParameterNames", &cpptsemKalmanModel::setDiscreteTimeParameterNames, "Set up the names of the discrete time parameters of the Kalman model")
  .method( "computeAndFitKalman", &cpptsemKalmanModel::computeAndFitKalman, "Computes the Kalman matrices")
  .method( "approxKalmanGradients", &cpptsemKalmanModel::approxKalmanGradients, "Returns a central approximation of the gradients (Jacobian).")
  ;
}


