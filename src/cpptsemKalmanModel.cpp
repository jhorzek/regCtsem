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
void cpptsemKalmanModel::setUpdate(bool mUpdate){
  // Should the predictions be updated in the Kalman filter?
  update = mUpdate;
}

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

void cpptsemKalmanModel::setKalmanGroupings(arma::colvec mPersons, arma::colvec mGroup){
  persons = mPersons;
  group = mGroup;
}

void cpptsemKalmanModel::setParameterValues(Rcpp::NumericVector mParameters, Rcpp::StringVector parameterLabels){
  Rcpp::StringVector personExpectedParameterLabels, personMatrixLabels;
  Rcpp::NumericVector personCurrentParameters, personRow, personCol;
  Rcpp::LogicalVector groupSelector, personChanged;

  // change status of Kalman computation and Kalman fitting
  computeWasCalled = false;
  // set the parameters in (1) the parameterTable (2) the ctMatrices
  bool wasSet; // checks if all parameters were changed
  Rcpp::StringVector expectedParameterLabels = parameterTable["label"];
  Rcpp::StringVector matrixLabels = parameterTable["matrix"]; // vector with labels of the ct matrices
  Rcpp::NumericVector currentParameters = parameterTable["value"]; // values before change
  Rcpp::NumericVector row =  parameterTable["row"];
  Rcpp::NumericVector col =  parameterTable["col"];
  Rcpp::NumericVector groupID =  parameterTable["groupID"];
  Rcpp::NumericVector uniqueIDs = unique(groupID);
  Rcpp::LogicalVector changed = parameterTable["changed"];

  for(int group = 0; group < uniqueIDs.length(); group++){
    // select the rows of the specified group
    groupSelector = (groupID==uniqueIDs[group]);

    personExpectedParameterLabels = expectedParameterLabels[groupSelector];
    personMatrixLabels = matrixLabels[groupSelector];
    personCurrentParameters = currentParameters[groupSelector];
    personRow = row[groupSelector];
    personCol = col[groupSelector];
    personChanged = changed[groupSelector];

    for(int i = 0; i < personExpectedParameterLabels.size(); i++){
      wasSet = false;
      for(int j = 0; j < parameterLabels.size(); j++){
        if(parameterLabels(j) == personExpectedParameterLabels(i)){
          wasSet = true;
          if(!(personCurrentParameters(i) == mParameters(j))){
            personChanged.fill(true);
          }
          personCurrentParameters(i) = mParameters(j);
        }
      }

      if(!wasSet){
        Rcpp::stop("Error while setting parameters: The parameter" +  personExpectedParameterLabels(i) + " was not found!");
      }
    }

    currentParameters[groupSelector] = personCurrentParameters;
    changed[groupSelector] = personChanged;
  }

  // set parameters to values of last group
  for(int i = 0; i < personExpectedParameterLabels.size(); i++){
    Rcpp::String matrixLabel = personMatrixLabels(i); // matrix in which parameter i is situated
    Rcpp::List currentMatrixList = ctMatrixList[matrixLabel]; // values and labels of this matrix
    Rcpp::NumericMatrix currentMatrixValues = currentMatrixList["values"]; // values of this matrix
    currentMatrixValues(personRow(i), personCol(i)) = personCurrentParameters(i); // set the values to the new parameter values
  }

}

void cpptsemKalmanModel::setKalmanMatrixValues(int selectedGroup){
  // sets the values of the drift, diffusion, etc matrices to the values of a specific group.
  Rcpp::StringVector personExpectedParameterLabels, personMatrixLabels;
  Rcpp::NumericVector personCurrentParameters, personRow, personCol;
  Rcpp::LogicalVector groupSelector;

  // set the parameters in the ctMatrices
  Rcpp::StringVector expectedParameterLabels = parameterTable["label"];
  Rcpp::StringVector matrixLabels = parameterTable["matrix"]; // vector with labels of the ct matrices
  Rcpp::NumericVector currentParameters = parameterTable["value"]; // values before change
  Rcpp::NumericVector row =  parameterTable["row"];
  Rcpp::NumericVector col =  parameterTable["col"];
  Rcpp::NumericVector groupID =  parameterTable["groupID"];

  // select the rows of the specified group
  groupSelector = (groupID==selectedGroup);

  personExpectedParameterLabels = expectedParameterLabels[groupSelector];
  personMatrixLabels = matrixLabels[groupSelector];
  personCurrentParameters = currentParameters[groupSelector];
  personRow = row[groupSelector];
  personCol = col[groupSelector];

  // set parameters to values of last group
  for(int i = 0; i < personExpectedParameterLabels.size(); i++){
    Rcpp::String matrixLabel = personMatrixLabels(i); // matrix in which parameter i is situated
    Rcpp::List currentMatrixList = ctMatrixList[matrixLabel]; // values and labels of this matrix
    Rcpp::NumericMatrix currentMatrixValues = currentMatrixList["values"]; // values of this matrix
    currentMatrixValues(personRow(i), personCol(i)) = personCurrentParameters(i); // set the values to the new parameter values
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

void cpptsemKalmanModel::setKalmanData(Rcpp::List mKalmanData, bool mhasDefinitionVariables){
  kalmanData = Rcpp::as<arma::mat>(mKalmanData["dataset"]);
  sampleSize = mKalmanData["sampleSize"];
  Tpoints = mKalmanData["Tpoints"];
  nlatent = mKalmanData["nlatent"];
  nmanifest = mKalmanData["nmanifest"];
  hasDefinitionVariables = mhasDefinitionVariables;
  arma::colvec tempIndM2ll(sampleSize);
  tempIndM2ll.fill(arma::datum::nan);
  indM2LL = tempIndM2ll;
}

void cpptsemKalmanModel::setDiscreteTimeParameterNames(Rcpp::List mDiscreteTimeParameterNames){
  discreteTimeParameterNames = mDiscreteTimeParameterNames;
}

void cpptsemKalmanModel::computeAndFitKalman(){
  computeWasCalled = true;
  arma::mat DRIFTInverseValues;
  arma::mat DRIFTHASH;
  arma::mat DRIFTHASHInverse;
  arma::uvec currentPersonRows;
  int currentSampleSize;
  arma::mat currentDataSet, currentLatentScores, currentPredictedManifest, currentIndM2LL;

  // iterate over all unique sets of parameters
  Rcpp::LogicalVector changed = parameterTable["changed"];
  Rcpp::LogicalVector currentParameters, currentChange;
  Rcpp::NumericVector groupID =  parameterTable["groupID"];
  Rcpp::NumericVector uniqueIDs = unique(groupID);

  for(int idIndex = 0; idIndex < uniqueIDs.length(); idIndex++){
    // check if the parameters of the current id-group changed. If not, we can reuse the previous -2log Likelihood
    currentParameters = groupID == uniqueIDs[idIndex];
    currentChange = changed[currentParameters];
    if(Rcpp::is_false(Rcpp::any(currentChange))){
      // if none of the parameters changed: skip to next group
      continue;
    }
    // set parameters to values of persons with current id
    setKalmanMatrixValues(uniqueIDs[idIndex]);
    currentPersonRows = arma::find(group == uniqueIDs[idIndex]);
    currentSampleSize = currentPersonRows.n_elem;
    currentDataSet = kalmanData.rows(currentPersonRows);

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
      Rcpp::stop("Could not find a DRIFT!");
    }

    if( ctMatrixList.containsElementNamed("DIFFUSIONbase") ){
      // get DIFFUSION
      Rcpp::List DIFFUSIONbaseList = ctMatrixList["DIFFUSIONbase"];
      DIFFUSIONbaseValues = Rcpp::as<arma::mat>(DIFFUSIONbaseList["values"]);
      DIFFUSIONValues = getVarianceFromVarianceBase(DIFFUSIONbaseValues);
    }else{
      Rcpp::stop("Could not find a DIFFUSION!");
    }

    if( ctMatrixList.containsElementNamed("T0VARbase") ){
      // get T0VAR
      Rcpp::List T0VARbaseList = ctMatrixList["T0VARbase"];
      arma::mat T0VARbaseValues = T0VARbaseList["values"];
      T0VARValues = getVarianceFromVarianceBase(T0VARbaseValues);
    }

    // compute asymptotic diffusion
    asymptoticDIFFUSION = -1*arma::inv(DRIFTHASH) * arma::vectorise(DIFFUSIONValues);
    asymptoticDIFFUSION.reshape(DIFFUSIONValues.n_rows, DIFFUSIONValues.n_cols);
    if(stationaryT0VAR){
      T0VARValues = asymptoticDIFFUSION;
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
      Rcpp::stop("discreteDRIFTUnique seems to be missing!");
    }

    // compute discrete time traits
    if(hasDiscreteTRAITUnique){
      discreteTRAITUnique = computeDiscreteTRAITs(discreteDRIFTUnique, discreteTRAITUnique);
    }

    // compute discrete DRIFTHASHs
    if(hasDRIFTHASHExponentialUnique){
      DRIFTHASHExponentialUnique = computeDRIFTHASHExponentials(DRIFTHASH, DRIFTHASHExponentialUnique);
    }else{
      Rcpp::stop("DRIFTHASHExponentialUnique seems to be missing!");
    }

    // compute discreteDIFFUSION
    if(hasDiscreteDIFFUSIONUnique){
      if(!hasDiscreteDIFFUSIONUnique){
        Rcpp::stop("Error: DRIFTHASHExponentialUnique missing!");
      }
      discreteDIFFUSIONUnique = computeDiscreteDIFFUSIONs(DRIFTHASHInverse, DIFFUSIONValues,
                                                          DRIFTHASHExponentialUnique, discreteDIFFUSIONUnique);
    }

    // compute discreteCINT
    if(hasDiscreteCINTUnique){
      if(!hasDiscreteDRIFTUnique){
        Rcpp::stop("Error: hasDiscreteDRIFTUnique missing!");
      }
      Rcpp::List CINTList = ctMatrixList["CINT"];
      arma::colvec CINTValues = Rcpp::as<arma::colvec>(CINTList["values"]);

      discreteCINTUnique = computeDiscreteCINTs(discreteCINTUnique, DRIFTInverseValues, discreteDRIFTUnique, CINTValues);
    }

    currentLatentScores.resize(currentSampleSize, latentScores.n_cols);
    currentLatentScores.fill(arma::datum::nan);
    currentPredictedManifest.resize(currentSampleSize, predictedManifestValues.n_cols);
    currentPredictedManifest.fill(arma::datum::nan);

    // Kalman filter: Prediction and Updating
    currentIndM2LL = kalmanFit(update,
                               currentSampleSize,
                               Tpoints,
                               nlatent,
                               nmanifest,
                               currentDataSet,
                               currentLatentScores,
                               currentPredictedManifest,
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
    indM2LL.rows(currentPersonRows) = currentIndM2LL;
    latentScores.rows(currentPersonRows) = currentLatentScores;
    predictedManifestValues.rows(currentPersonRows) = currentPredictedManifest;
  }
  m2LL = sum(indM2LL);

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

// compute single gradient value. Returns the central gradient approximation
Rcpp::NumericVector cpptsemKalmanModel::approxKalmanGradient(double epsilon, Rcpp::String parName){
  Rcpp::NumericVector parameterValues = getParameterValues();
  Rcpp::StringVector parameterNames = parameterValues.names();
  Rcpp::NumericMatrix likelihoods( 1 , 2 );
  Rcpp::NumericVector gradient(1);

  Rcpp::StringVector colNames = {"stepBackward", "stepForward"};
  colnames(likelihoods) = colNames;
  gradient.names() = parName;

  // avoid changing the parameterValues by reference:
  Rcpp::NumericVector currentParameterValues = Rcpp::clone( parameterValues );

  // forward step
  currentParameterValues(parName) = currentParameterValues(parName) + epsilon;
  setParameterValues(currentParameterValues, parameterNames);
  computeAndFitKalman();
  likelihoods(0,0) = m2LL;

  // reset parameter
  currentParameterValues(parName) = currentParameterValues(parName) - epsilon;
  setParameterValues(currentParameterValues, parameterNames);

  // backward step
  currentParameterValues(parName) = currentParameterValues(parName) - epsilon;
  setParameterValues(currentParameterValues, parameterNames);
  computeAndFitKalman();
  likelihoods(0,1) = m2LL;

  // reset parameter
  currentParameterValues(parName) = currentParameterValues(parName) + epsilon;
  setParameterValues(currentParameterValues, parameterNames);
  gradient(0) = (likelihoods(0,0) - likelihoods(0,1))/(2*epsilon);


  return(Rcpp::clone(gradient));
}

Rcpp::List cpptsemKalmanModel::GIST(
    Rcpp::NumericVector pars,
    Rcpp::StringVector regIndicators,
    Rcpp::NumericVector lambda,
    double eta, double sig,
    Rcpp::NumericVector gradient_epsilons,
    double initialStepsize, double stepsizeMin, double stepsizeMax,
    std::string GISTLinesearchCriterion,
    int GISTNonMonotoneNBack,
    int maxIter_out, int maxIter_in,
    double break_outer, std::string break_crit, int verbose){

  if(!(pars.length() == lambda.length())){Rcpp::stop("Lambda has to be a vector of the same length as pars.");}

  int k_out = 0, k = 0;
  bool convergence = false, breakCriterion = false,
    breakOuter = false, smallStep = false, largeStep = false, tempComparison = false;
  double gradient_epsilon, m2LL_k, m2LL_kp1, regM2LL_k, regM2LL_kp1, lambda_i;
  arma::mat stepsize(1,1, arma::fill::zeros);
  Rcpp::String parameterName;
  Rcpp::NumericVector regM2LLs;

  // initialize parameters
  Rcpp::StringVector parameterNames = pars.names();
  setParameterValues(pars, parameterNames);
  Rcpp::NumericVector gradients_km1, gradients_k, gradients_kp1, subgradients,
  parameters_km1 = Rcpp::clone(pars), parameters_k = Rcpp::clone(pars), parameters_kp1 = Rcpp::clone(pars),
  x_kRcpp, y_kRcpp, u_kRcpp;
  arma::colvec x_k, y_k;
  arma::rowvec u_k;

  parameters_km1.fill(arma::datum::nan);
  parameters_kp1.fill(arma::datum::nan);

  // get initial fit and gradients
  computeAndFitKalman();

  if (indM2LL.has_nan()){
    Rcpp::stop("Infeasible initial values: Could not compute gradients.");
  }

  m2LL_k = m2LL;
  regM2LL_k = m2LL + computePenalty(parameters_k,
                                    regIndicators,
                                    lambda);
  try{
    for(int i = 0; i < gradient_epsilons.length(); i++){
      gradient_epsilon = gradient_epsilons[i];
      gradients_k = approxKalmanGradients(gradient_epsilon);
      tempComparison = Rcpp::any(Rcpp::is_na(gradients_k));
      if (!tempComparison){
        break;
      }
      if(i == gradient_epsilons.length()-1){
        Rcpp::stop("Infeasible initial values: Could not compute gradients.");
      }
    }
  }catch(...){
    Rcpp::stop("Infeasible initial values: Could not compute gradients.");
  }

  while(k_out < maxIter_out){
    Rcpp::checkUserInterrupt();
    k_out += 1;

    // set step size
    if(k_out == 1){
      stepsize = initialStepsize;
    }else{
      x_k = parameters_k - parameters_km1;
      y_k = gradients_k - gradients_km1;
      x_kRcpp = x_k;
      y_kRcpp = y_k;
      x_kRcpp.names() = parameterNames;
      y_kRcpp.names() = parameterNames;

      stepsize = (arma::trans(x_k)*y_k)/(arma::trans(x_k)*x_k);

      if(stepsize.has_nan() || stepsize.has_inf()){
        stepsize = initialStepsize;
      }

      smallStep = stepsize(0,0) < stepsizeMin;
      largeStep = stepsize(0,0) > stepsizeMax;
      if(smallStep || largeStep){
        stepsize = initialStepsize;
      }
    }// end step size setting

    // inner iterations
    while(k < maxIter_in){
      //Rcpp::Rcout << "stepsize: " << stepsize(0,0) << std::endl;
      //Rcpp::Rcout << "parameters_k: " << parameters_k << std::endl;
      //Rcpp::Rcout << "gradients_k: " << gradients_k << std::endl;
      u_k = parameters_k - gradients_k/(stepsize(0,0));
      u_kRcpp = u_k;
      u_kRcpp.names() = parameterNames;

      parameters_kp1.fill(arma::datum::nan);

      for(int i = 0; i < parameters_kp1.length(); i++){
        parameterName = parameterNames[i];
        lambda_i = lambda[parameterName];

        bool isRegularized = false;
        for(int j = 0; j < regIndicators.length(); j++){
          isRegularized = (parameterName == regIndicators[j]);
          if(isRegularized){break;}
        }

        if(isRegularized){
          // update parameter i with lasso
          parameters_kp1[parameterName] = std::copysign(1.0, u_kRcpp[parameterName])*std::max(0.0, std::abs(u_kRcpp[parameterName]) - lambda_i/stepsize(0,0));
        }else{
          parameters_kp1[parameterName] = u_kRcpp[parameterName];
        }
      }

      try{
        // fit model with new parameters
        setParameterValues(parameters_kp1, parameterNames);

        // get initial fit and gradients
        computeAndFitKalman();
      }catch(...){
        // update step size
        stepsize = eta*stepsize;

        // update iteration counter
        k +=1;

        // skip rest
        continue;
      }

      if (indM2LL.has_nan()){
        // update step size
        stepsize = eta*stepsize;

        // update iteration counter
        k +=1;

        // skip rest
        continue;
      }

      m2LL_kp1 = m2LL;
      regM2LL_kp1 = m2LL_kp1 + computePenalty(parameters_kp1,
                                              regIndicators,
                                              lambda);

      if(!arma::is_finite(m2LL_kp1)){

        // update step size
        stepsize = eta*stepsize;

        // update iteration counter
        k +=1;

        // skip rest
        continue;
      }

      // break if line search condition is satisfied
      if(GISTLinesearchCriterion == "monotone"){
        arma::colvec parameterDiff = parameters_k - parameters_kp1;
        arma::mat sumSquared(1,1);
        sumSquared = sum(arma::pow(parameterDiff, 2));
        double comparisonValue = regM2LL_k - (sig/2) * stepsize(0,0) * sumSquared(0,0);
        breakCriterion = regM2LL_kp1 < comparisonValue;
      }else if(GISTLinesearchCriterion == "non-monotone"){
        Rcpp::stop("Non-monotone linesearch not implemented in C++");
        //int nBack = std::max(1,k_out-GISTNonMonotoneNBack);
        //arma::colvec parameterDiff = parameters_k - parameters_kp1;
        //breakCriterion <- regM2LL_kp1 <= max(regM2LL[nBack:k_out]) - (sig/2) * stepsize * sum((parameters_k - parameters_kp1)^2)
      }else{
        Rcpp::stop("Unknown GISTLinesearchCriterion. Possible are monotone and non-monotone.");
      }
      if(breakCriterion){
        break;
      }

      // update step size
      stepsize = eta*stepsize;

      // update iteration counter
      k +=1;

    } // end inner iteration

    // Compute gradients
    try{
      for(int i = 0; i < gradient_epsilons.length(); i++){
        gradient_epsilon = gradient_epsilons[i];
        gradients_kp1 = approxKalmanGradients(gradient_epsilon);
        tempComparison = Rcpp::any(Rcpp::is_na(gradients_kp1));
        if (!tempComparison){
          break;
        }
        if(i == gradient_epsilons.length()-1){
          Rcpp::stop("No gradients in outer iteration");
        }
      }
    }catch(...){
      Rcpp::stop("No gradients in outer iteration");
    }

    // update parameters for next iteration
    parameters_km1 = Rcpp::clone(parameters_k);
    parameters_k = Rcpp::clone(parameters_kp1);
    gradients_km1 = Rcpp::clone(gradients_k);
    gradients_k = Rcpp::clone(gradients_kp1);
    m2LL_k = m2LL_kp1;
    regM2LL_k = regM2LL_kp1;
    //regM2LLs.push_back(regM2LL_kp1);

    if(verbose == 1){
      Rcpp::Rcout << "Parameter in iteration " << k_out << ": " << parameters_k << std::endl;
    }

    //regM2LL[k_out] = regM2LL_k

    // break outer loop if stopping criterion is satisfied
    if(break_crit == "gradient"){

      // Gradient based break condition: Problem: gradients extremely dependent on the epsilon in the approximation

      subgradients = computeSubgradients(parameters_k, gradients_k, regIndicators, lambda);

      breakOuter = Rcpp::max(Rcpp::abs(subgradients)) < break_outer;

    }else if(break_crit == "parameterChange"){
      arma::colvec parameterDiff = parameters_k - parameters_km1;
      arma::colvec temppars = parameters_km1;
      arma::mat sumSquared1(1,1);
      arma::mat sumSquared2(1,1);
      sumSquared1 = sum(arma::pow(parameterDiff, 2));
      sumSquared2 = sum(arma::pow(temppars, 2));
      double comparisonValue = sqrt(sumSquared1(0,0))/sqrt(sumSquared2(0,0));
      breakOuter = comparisonValue < break_outer;

    }else{
      Rcpp::stop("Unknown breaking criterion. The value passed to break_crit has to be named (either 'parameterChange' or 'gradient') See ?controlGIST for details on break_outer");
    }
    if(breakOuter){
      //Rcpp::Rcout << "Break outer" << std::endl;
      convergence = true;
      break;
    }

  }// end outer loop

  if(k_out == maxIter_out){
    Rcpp::Rcout << "Warning: Maximal number if outer iterations used. Increase number of outer iterations or take smaller lambda-steps.";
  }

  subgradients = computeSubgradients(parameters_k, gradients_k, regIndicators, lambda);

  Rcpp::List retList = Rcpp::List::create(Named("par") = parameters_k , _["m2LL"] = m2LL, _["subgradients"] = subgradients, _["regM2LLs"] = regM2LLs, _["converged"] = convergence);

  return(retList);
}// end GIST


double cpptsemKalmanModel::computePenalty(Rcpp::NumericVector pars,
                                          Rcpp::StringVector regIndicators,
                                          Rcpp::NumericVector lambda){
  Rcpp::String currentParName;
  double penaltyValue = 0.0;

  for(int i = 0; i < regIndicators.length(); i++){
    currentParName = regIndicators[i];
    penaltyValue += lambda[currentParName]*std::abs(pars[currentParName]);
  }
  return(penaltyValue);
}

Rcpp::NumericVector cpptsemKalmanModel::computeSubgradients(Rcpp::NumericVector pars, Rcpp::NumericVector gradients, Rcpp::StringVector regIndicators, Rcpp::NumericVector lambdas){
  double absoluteValueOf, penaltyGradient, penaltyGradientLower, penaltyGradientUpper;
  bool setZero;
  // extract names
  Rcpp::String parameterNames = pars.names();
  Rcpp::String currentParLabel;
  // first part: derivative of Likelihood
  Rcpp::NumericVector subgradient = Rcpp::clone(gradients);

  // second part: derivative of penalty term
  for(int i = 0; i < regIndicators.length(); i++){
    currentParLabel = regIndicators(i);
    absoluteValueOf = pars[currentParLabel];

    if(!(absoluteValueOf == 0)){
      penaltyGradient = lambdas[currentParLabel]*std::copysign(1.0, absoluteValueOf);
      subgradient[currentParLabel] = subgradient[currentParLabel] + penaltyGradient;
    }else{
      penaltyGradientLower = -1.0*lambdas[currentParLabel];
      penaltyGradientUpper = lambdas[currentParLabel];
      // check if likelihood gradient is within interval:
      setZero = (subgradient[currentParLabel] > penaltyGradientLower) & (subgradient[currentParLabel] < penaltyGradientUpper);
      if(setZero){
        subgradient[currentParLabel] = 0;
      }else{
        subgradient[currentParLabel] = std::copysign(1.0, subgradient[currentParLabel])*(std::abs(subgradient[currentParLabel]) - lambdas[currentParLabel]);
      }
    } // end is zero
  } // end for
  return(subgradient);
}

RCPP_EXPOSED_CLASS(cpptsemKalmanModel)
  RCPP_MODULE(cpptsemKalmanModel_cpp){
    using namespace Rcpp;
    Rcpp::class_<cpptsemKalmanModel>( "cpptsemKalmanModel" )
      .constructor<std::string, Rcpp::List, Rcpp::DataFrame, bool, bool>("Creates a new model. Expects a model name, ctMatrixList, parameterTable")
      .field_readonly( "name", &cpptsemKalmanModel::name, "Name of the model")
      .field_readonly( "update", &cpptsemKalmanModel::update, "Are the latent scores updated in the Kalman filter?")
      .field_readonly( "ctMatrixList", &cpptsemKalmanModel::ctMatrixList, "List of ct matrices")
      .field_readonly( "parameterTable", &cpptsemKalmanModel::parameterTable, "Data frame of model parameters")
      .field_readonly( "discreteDRIFTUnique", &cpptsemKalmanModel::discreteDRIFTUnique, "discreteDRIFTUnique")
      .field_readonly( "discreteTRAITUnique", &cpptsemKalmanModel::discreteTRAITUnique, "discreteTRAITUnique")
      .field_readonly( "discreteCINTUnique", &cpptsemKalmanModel::discreteCINTUnique, "discreteCINTUnique")
      .field_readonly( "DRIFTHASHExponentialUnique", &cpptsemKalmanModel::DRIFTHASHExponentialUnique, "DRIFTHASHExponentialUnique")
      .field_readonly( "discreteDIFFUSIONUnique", &cpptsemKalmanModel::discreteDIFFUSIONUnique, "discreteDIFFUSIONUnique")
      .field_readonly( "DRIFTValues", &cpptsemKalmanModel::DRIFTValues, "DRIFTValues")
      .field_readonly( "DIFFUSIONValues", &cpptsemKalmanModel::DIFFUSIONValues, "DIFFUSIONValues")
      .field_readonly( "DIFFUSIONbaseValues", &cpptsemKalmanModel::DIFFUSIONbaseValues, "DIFFUSIONbaseValues")
      .field_readonly( "T0VARValues", &cpptsemKalmanModel::T0VARValues, "T0VARValues")
      .field_readonly( "T0MEANSValues", &cpptsemKalmanModel::T0MEANSValues, "T0MEANSValues")
      .field_readonly( "TRAITVARValues", &cpptsemKalmanModel::TRAITVARValues, "TRAITVARValues")
      .field_readonly( "asymptoticDIFFUSION", &cpptsemKalmanModel::asymptoticDIFFUSION, "asymptoticDIFFUSION")
      .field_readonly( "MANIFESTMEANSValues", &cpptsemKalmanModel::MANIFESTMEANSValues, "MANIFESTMEANSValues")
      .field_readonly( "MANIFESTVARValues", &cpptsemKalmanModel::MANIFESTVARValues, "MANIFESTVARValues")
      .field_readonly( "LAMBDAValues", &cpptsemKalmanModel::LAMBDAValues, "LAMBDAValues")
      .field_readonly( "m2LL", &cpptsemKalmanModel::m2LL, "-2 log likelihood")
      .field_readonly( "indM2LL", &cpptsemKalmanModel::indM2LL, "individual -2 log likelihood")
      .field_readonly( "kalmanData", &cpptsemKalmanModel::kalmanData, "Data for Kalman model")
      .field_readonly( "latentScores", &cpptsemKalmanModel::latentScores, "latentScores")
      .field_readonly( "predictedManifestValues", &cpptsemKalmanModel::predictedManifestValues, "predictedManifestValues")


    // methods
    .method( "setUpdate", &cpptsemKalmanModel::setUpdate, "Should the latent scores be updated in the Kalman Filter procedure? Expects a boolean.")
    .method( "setParameterValues", &cpptsemKalmanModel::setParameterValues, "Set the parameters. Expects a vector with parametervalues and a stringvector with labels")
    .method( "setDiscreteDRIFTUnique", &cpptsemKalmanModel::setDiscreteDRIFTUnique, "Set the number of discrete time drifts, corresponding dts and values")
    .method( "setDiscreteTRAITUnique", &cpptsemKalmanModel::setDiscreteTRAITUnique, "Set the number of discrete time traits, corresponding dts and values")
    .method( "setDRIFTHASHExponentialUnique", &cpptsemKalmanModel::setDRIFTHASHExponentialUnique, "Set the number of drift hash exponentials, corresponding dts and values")
    .method( "setDiscreteDIFFUSIONUnique", &cpptsemKalmanModel::setDiscreteDIFFUSIONUnique, "Set the number of discrete diffusions, corresponding dts and values")
    .method( "setDiscreteCINTUnique", &cpptsemKalmanModel::setDiscreteCINTUnique, "Set the number of discrete continuous time intercepts, corresponding dts and values")
    .method( "getParameterValues", &cpptsemKalmanModel::getParameterValues, "Get current parameter values")
    .method( "setKalmanMatrices", &cpptsemKalmanModel::setKalmanMatrices, "Set up the Kalman matrices. Requires numeric matrices for A, S, M, F and DataFrames with cpp compatible row and column indicators for A, S, M which specify where in the matrices the discrete time parameters go.")
    .method( "setKalmanData", &cpptsemKalmanModel::setKalmanData, "Set up the dataset for Kalman models. Expects a list with sampleSize, nObservedVariables, rawData, observedMeans, observedCov, expectedSelector, expectedMeans, expectedCovariance, m2LL for each missingness pattern")
    .method( "setKalmanGroupings", &cpptsemKalmanModel::setKalmanGroupings, "Set groups. Expects a vector with person-indices from 0 to nPersons and a vector with group indices")
    .method( "setDiscreteTimeParameterNames", &cpptsemKalmanModel::setDiscreteTimeParameterNames, "Set up the names of the discrete time parameters of the Kalman model")
    .method( "computeAndFitKalman", &cpptsemKalmanModel::computeAndFitKalman, "Computes the Kalman matrices")
    .method( "approxKalmanGradients", &cpptsemKalmanModel::approxKalmanGradients, "Returns a central approximation of the gradients.")
    .method( "approxKalmanGradient", &cpptsemKalmanModel::approxKalmanGradient, "Returns a central approximation of a single gradient.")
    .method( "GIST", &cpptsemKalmanModel::GIST, "Optimizes with GIST")
    .method( "computePenalty", &cpptsemKalmanModel::computePenalty, "Computes penalty value")
    .method( "computeSubgradients", &cpptsemKalmanModel::computeSubgradients, "Computes subgradients")
    ;
  }


