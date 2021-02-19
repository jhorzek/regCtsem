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
                                       bool mStationaryT0MEANS): cpptsemmodel(mName, mCtMatrixList, mParameterTable,
                                       mStationaryT0VAR, mStationaryT0MEANS){}

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

RCPP_MODULE(cpptsemKalmanmodule){
  Rcpp::class_<cpptsemmodel>( "cpptsemmodel" )

  .field_readonly( "name", &cpptsemmodel::name, "Name of the model")
  .field_readonly( "ctMatrixList", &cpptsemmodel::ctMatrixList, "List of ct matrices")
  .field_readonly( "parameterTable", &cpptsemmodel::parameterTable, "Data frame of model parameters")
  .field_readonly( "discreteDRIFTUnique", &cpptsemmodel::discreteDRIFTUnique, "discreteDRIFTUnique")
  .field_readonly( "discreteTRAITUnique", &cpptsemmodel::discreteTRAITUnique, "discreteTRAITUnique")
  .field_readonly( "discreteCINTUnique", &cpptsemmodel::discreteCINTUnique, "discreteCINTUnique")
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

  // additional fields:
  Rcpp::class_<cpptsemKalmanModel>( "cpptsemKalmanModel" )
    .derives<cpptsemmodel>("cpptsemmodel")
    .constructor<std::string, Rcpp::List, Rcpp::DataFrame, bool, bool>("Creates a new model. Expects a model name, ctMatrixList, parameterTable")
    .field_readonly( "kalmanData", &cpptsemKalmanModel::kalmanData, "Data for Kalman model")
    .field_readonly( "latentScores", &cpptsemKalmanModel::latentScores, "latentScores")
    .field_readonly( "predictedManifestValues", &cpptsemKalmanModel::predictedManifestValues, "predictedManifestValues")

  .method( "setKalmanMatrices", &cpptsemKalmanModel::setKalmanMatrices, "Set up the Kalman matrices. Requires numeric matrices for A, S, M, F and DataFrames with cpp compatible row and column indicators for A, S, M which specify where in the matrices the discrete time parameters go.")
  .method( "setKalmanData", &cpptsemKalmanModel::setKalmanData, "Set up the dataset for Kalman models. Expects a list with sampleSize, nObservedVariables, rawData, observedMeans, observedCov, expectedSelector, expectedMeans, expectedCovariance, m2LL for each missingness pattern")
  .method( "setDiscreteTimeParameterNames", &cpptsemKalmanModel::setDiscreteTimeParameterNames, "Set up the names of the discrete time parameters of the Kalman model")
  .method( "computeAndFitKalman", &cpptsemKalmanModel::computeAndFitKalman, "Computes the Kalman matrices")
  .method( "approxKalmanGradients", &cpptsemKalmanModel::approxKalmanGradients, "Returns a central approximation of the gradients (Jacobian).")
  ;
}


