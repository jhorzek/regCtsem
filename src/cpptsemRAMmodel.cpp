#include "cpptsemRAMmodel.h"
#include "computeDiscreteParameters.h"
#include "getVarianceFromVarianceBase.h"
#include "fillRAMMatrices.h"
#include "computeRAMExpectations.h"
#include "computeRAMM2LL.h"
#include <RcppArmadillo.h>

// [[Rcpp :: depends ( RcppArmadillo )]]

using namespace Rcpp;

// Constructor:

cpptsemRAMmodel::cpptsemRAMmodel(std::string mName,
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

// general functions
void cpptsemRAMmodel::setDiscreteDRIFTUnique(Rcpp::List mDiscreteDRIFTUnique){
  // list with labels of the unique discrete drifts, corresponding time intervals and results
  discreteDRIFTUnique = mDiscreteDRIFTUnique;
  hasDiscreteDRIFTUnique = true;
}
void cpptsemRAMmodel::setDiscreteTRAITUnique(Rcpp::List mDiscreteTRAITUnique){
  // list with labels of the unique discrete traits, corresponding time intervals and results
  discreteTRAITUnique = mDiscreteTRAITUnique;
  hasDiscreteTRAITUnique = true;
}

void cpptsemRAMmodel::setDRIFTHASHExponentialUnique(Rcpp::List mDRIFTHASHExponentialUnique){
  // list with labels of the unique discrete traits, corresponding time intervals and results
  DRIFTHASHExponentialUnique = mDRIFTHASHExponentialUnique;
  hasDRIFTHASHExponentialUnique = true;
}

void cpptsemRAMmodel::setDiscreteDIFFUSIONUnique(Rcpp::List mDiscreteDIFFUSIONUnique){
  discreteDIFFUSIONUnique = mDiscreteDIFFUSIONUnique;
  hasDiscreteDIFFUSIONUnique = true;
}

void cpptsemRAMmodel::setDiscreteCINTUnique(Rcpp::List mDiscreteCINTUnique){
  discreteCINTUnique = mDiscreteCINTUnique;
  hasDiscreteCINTUnique = true;
}

void cpptsemRAMmodel::setParameterValues(Rcpp::NumericVector mParameters, Rcpp::StringVector parameterLabels){
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

// returns the parameter values of a cpptsemRAMmodel
Rcpp::NumericVector cpptsemRAMmodel::getParameterValues() {

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

// RAM specific functions

void cpptsemRAMmodel::setRAMMatrices(arma::mat mA,
                                     arma::mat mS,
                                     arma::mat mM,
                                     arma::mat mF,
                                     Rcpp::DataFrame mCppAParameterIndicators,
                                     Rcpp::DataFrame mCppSParameterIndicators,
                                     Rcpp::DataFrame mCppMParameterIndicators){
  A = mA;
  S = mS;
  M = mM;
  F = mF;
  cppAParameterIndicators = mCppAParameterIndicators;
  cppSParameterIndicators = mCppSParameterIndicators;
  cppMParameterIndicators = mCppMParameterIndicators;

}

void cpptsemRAMmodel::setRAMData(Rcpp::List mRAMdata){
  RAMdata = mRAMdata;
}

void cpptsemRAMmodel::computeRAM(){
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
    DIFFUSIONbaseValues = Rcpp::as<arma::mat>(DIFFUSIONbaseList["values"]);
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

  if( ctMatrixList.containsElementNamed("TRAITVARbase") ){
    // get TRAITVAR
    Rcpp::List TRAITVARbaseList = ctMatrixList["TRAITVARbase"];
    arma::mat TRAITVARbaseValues = TRAITVARbaseList["values"];
    TRAITVARValues = getVarianceFromVarianceBase(TRAITVARbaseValues);
  }

  if(stationaryT0VAR){
    // compute asymptotic diffusion
    arma::uword nrows = DIFFUSIONValues.n_rows;
    arma::uword ncols = DIFFUSIONValues.n_cols;
    arma::mat asymptoticDIFFUSION = -1*arma::inv(DRIFTHASH) * arma::vectorise(DIFFUSIONValues);
    T0VARValues = asymptoticDIFFUSION;
    T0VARValues.reshape(nrows, ncols);
  }

  if( ctMatrixList.containsElementNamed("MANIFESTVARbase") ){
    // MANIFESTVAR
    Rcpp::List MANIFESTVARbaseList = ctMatrixList["MANIFESTVARbase"];
    arma::mat MANIFESTVARbaseValues = MANIFESTVARbaseList["values"];
    MANIFESTVARValues = getVarianceFromVarianceBase(MANIFESTVARbaseValues);
  }

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

  // compute asymptotic diffusion
  //arma::mat asymptoticDIFFUSION = -1*arma::inv(DRIFTHASH) * arma::vectorise(DIFFUSIONValues);

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

  // Fill RAM Matrices

  // A
  Rcpp::List LAMBDAList = ctMatrixList["LAMBDA"];
  LAMBDAValues = Rcpp::as<arma::mat> (LAMBDAList["values"]);
  A = fillA(A,
            hasDiscreteDRIFTUnique, discreteDRIFTUnique,
            hasDiscreteTRAITUnique, discreteTRAITUnique,
            LAMBDAValues,
            cppAParameterIndicators);

  // S
  S = fillS(S, T0VARValues, MANIFESTVARValues,
            hasDiscreteTRAITUnique, TRAITVARValues,
            hasDiscreteDIFFUSIONUnique, discreteDIFFUSIONUnique,
            cppSParameterIndicators);

  // M
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

  M = fillM(M, T0MEANSValues, MANIFESTMEANSValues,
            hasDiscreteCINTUnique, discreteCINTUnique,
            cppMParameterIndicators);

  // Compute Expectations
  expectedCovariance = computeRAMExpectedCovariance(F, A, S);
  expectedMeans = computeRAMExpectedMeans(F, A, M);
}

// compute the -2 log likelihood
void cpptsemRAMmodel::fitRAM(){
  if(!computeWasCalled){
    Rcpp::Rcout << "Warning: fitRAM was called, but the RAM matrices have not been updated since the last change of paramters. Call computeRAM() first." << std::endl;
  }
  m2LL = computeRAMM2LL(RAMdata, expectedMeans, expectedCovariance);
}

// compute gradient vector. Returns the central gradient approximation
Rcpp::NumericVector cpptsemRAMmodel::approxRAMGradients(double epsilon){
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
    computeRAM();
    fitRAM();
    likelihoods(parameter,0) = m2LL;

    // reset parameter
    currentParameterValues(parameter) = currentParameterValues(parameter) - epsilon;
    setParameterValues(currentParameterValues, parameterNames);

    // backward step
    currentParameterValues(parameter) = currentParameterValues(parameter) - epsilon;
    setParameterValues(currentParameterValues, parameterNames);
    computeRAM();
    fitRAM();
    likelihoods(parameter,1) = m2LL;

    // reset parameter
    currentParameterValues(parameter) = currentParameterValues(parameter) + epsilon;
    setParameterValues(currentParameterValues, parameterNames);

    gradients(parameter) = (likelihoods(parameter,0) - likelihoods(parameter,1))/(2*epsilon);
  }

  return(Rcpp::clone(gradients));
}

// compute single gradient value. Returns the central gradient approximation
Rcpp::NumericVector cpptsemRAMmodel::approxRAMGradient(double epsilon, Rcpp::String parName){
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
  computeRAM();
  fitRAM();
  likelihoods(0,0) = m2LL;

  // reset parameter
  currentParameterValues(parName) = currentParameterValues(parName) - epsilon;
  setParameterValues(currentParameterValues, parameterNames);

  // backward step
  currentParameterValues(parName) = currentParameterValues(parName) - epsilon;
  setParameterValues(currentParameterValues, parameterNames);
  computeRAM();
  fitRAM();
  likelihoods(0,1) = m2LL;

  // reset parameter
  currentParameterValues(parName) = currentParameterValues(parName) + epsilon;
  setParameterValues(currentParameterValues, parameterNames);
  gradient(0) = (likelihoods(0,0) - likelihoods(0,1))/(2*epsilon);

  return(Rcpp::clone(gradient));
}


Rcpp::List cpptsemRAMmodel::GIST(
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
  computeRAM();
  fitRAM();

  if (!std::isfinite(m2LL)){
    Rcpp::stop("Infeasible initial values: Could not compute gradients.");
  }

  m2LL_k = m2LL;
  regM2LL_k = m2LL + computePenalty(parameters_k,
                                    regIndicators,
                                    lambda);
  try{
    for(int i = 0; i < gradient_epsilons.length(); i++){
      gradient_epsilon = gradient_epsilons[i];
      gradients_k = approxRAMGradients(gradient_epsilon);
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
        computeRAM();
        fitRAM();
      }catch(...){
        // update step size
        stepsize = eta*stepsize;

        // update iteration counter
        k +=1;

        // skip rest
        continue;
      }

      if (!std::isfinite(m2LL)){
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
        gradients_kp1 = approxRAMGradients(gradient_epsilon);
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


double cpptsemRAMmodel::computePenalty(Rcpp::NumericVector pars,
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

Rcpp::NumericVector cpptsemRAMmodel::computeSubgradients(Rcpp::NumericVector pars, Rcpp::NumericVector gradients, Rcpp::StringVector regIndicators, Rcpp::NumericVector lambdas){
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

RCPP_EXPOSED_CLASS(cpptsemRAMmodel)
RCPP_MODULE(cpptsemRAMmodel_cpp){
  using namespace Rcpp;
  Rcpp::class_<cpptsemRAMmodel>( "cpptsemRAMmodel" )
  .constructor<std::string, Rcpp::List, Rcpp::DataFrame, bool, bool>("Creates a new model. Expects a model name, ctMatrixList, parameterTable")
  .field_readonly( "cppAParameterIndicators", &cpptsemRAMmodel::cppAParameterIndicators, "cppAParameterIndicators")
  .field_readonly( "cppSParameterIndicators", &cpptsemRAMmodel::cppSParameterIndicators, "cppSParameterIndicators")
  .field_readonly( "cppMParameterIndicators", &cpptsemRAMmodel::cppMParameterIndicators, "cppMParameterIndicators")
  .field_readonly( "A", &cpptsemRAMmodel::A, "A matrix with directed effects in RAM notation")
  .field_readonly( "S", &cpptsemRAMmodel::S, "S matrix with covariances in RAM notation")
  .field_readonly( "M", &cpptsemRAMmodel::M, "M matrix with means in RAM notation")
  .field_readonly( "F", &cpptsemRAMmodel::F, "F matrix (filter matrix for separating manifest and latent variables) in RAM notation")
  .field_readonly( "expectedCovariance", &cpptsemRAMmodel::expectedCovariance, "Expected Covariance in RAM model")
  .field_readonly( "expectedMeans", &cpptsemRAMmodel::expectedMeans, "Expected Covariance in RAM model")
  .field_readonly( "RAMdata", &cpptsemRAMmodel::RAMdata, "Data for RAM model")
  .field_readonly( "name", &cpptsemRAMmodel::name, "Name of the model")
  .field_readonly( "ctMatrixList", &cpptsemRAMmodel::ctMatrixList, "List of ct matrices")
  .field_readonly( "parameterTable", &cpptsemRAMmodel::parameterTable, "Data frame of model parameters")
  .field_readonly( "discreteDRIFTUnique", &cpptsemRAMmodel::discreteDRIFTUnique, "discreteDRIFTUnique")
  .field_readonly( "discreteTRAITUnique", &cpptsemRAMmodel::discreteTRAITUnique, "discreteTRAITUnique")
  .field_readonly( "discreteCINTUnique", &cpptsemRAMmodel::discreteCINTUnique, "discreteCINTUnique")
  .field_readonly( "DRIFTHASHExponentialUnique", &cpptsemRAMmodel::DRIFTHASHExponentialUnique, "DRIFTHASHExponentialUnique")
  .field_readonly( "discreteDIFFUSIONUnique", &cpptsemRAMmodel::discreteDIFFUSIONUnique, "discreteDIFFUSIONUnique")
  .field_readonly( "DRIFTValues", &cpptsemRAMmodel::DRIFTValues, "DRIFTValues")
  .field_readonly( "DIFFUSIONValues", &cpptsemRAMmodel::DIFFUSIONValues, "DIFFUSIONValues")
  .field_readonly( "DIFFUSIONbaseValues", &cpptsemRAMmodel::DIFFUSIONbaseValues, "DIFFUSIONbaseValues")
  .field_readonly( "T0VARValues", &cpptsemRAMmodel::T0VARValues, "T0VARValues")
  .field_readonly( "T0MEANSValues", &cpptsemRAMmodel::T0MEANSValues, "T0MEANSValues")
  .field_readonly( "TRAITVARValues", &cpptsemRAMmodel::TRAITVARValues, "TRAITVARValues")
  .field_readonly( "MANIFESTMEANSValues", &cpptsemRAMmodel::MANIFESTMEANSValues, "MANIFESTMEANSValues")
  .field_readonly( "MANIFESTVARValues", &cpptsemRAMmodel::MANIFESTVARValues, "MANIFESTVARValues")
  .field_readonly( "LAMBDAValues", &cpptsemRAMmodel::LAMBDAValues, "LAMBDAValues")
  .field_readonly( "m2LL", &cpptsemRAMmodel::m2LL, "-2 log likelihood")

  // methods
  .method( "setParameterValues", &cpptsemRAMmodel::setParameterValues, "Set the parameters. Expects a vector with parametervalues and a stringvector with labels")
  .method( "setDiscreteDRIFTUnique", &cpptsemRAMmodel::setDiscreteDRIFTUnique, "Set the number of discrete time drifts, corresponding dts and values")
  .method( "setDiscreteTRAITUnique", &cpptsemRAMmodel::setDiscreteTRAITUnique, "Set the number of discrete time traits, corresponding dts and values")
  .method( "setDRIFTHASHExponentialUnique", &cpptsemRAMmodel::setDRIFTHASHExponentialUnique, "Set the number of drift hash exponentials, corresponding dts and values")
  .method( "setDiscreteDIFFUSIONUnique", &cpptsemRAMmodel::setDiscreteDIFFUSIONUnique, "Set the number of discrete diffusions, corresponding dts and values")
  .method( "setDiscreteCINTUnique", &cpptsemRAMmodel::setDiscreteCINTUnique, "Set the number of discrete continuous time intercepts, corresponding dts and values")
  .method( "getParameterValues", &cpptsemRAMmodel::getParameterValues, "Get current parameter values")
  .method( "setRAMMatrices", &cpptsemRAMmodel::setRAMMatrices, "Set up the RAM matrices. Requires numeric matrices for A, S, M, F and DataFrames with cpp compatible row and column indicators for A, S, M which specify where in the matrices the discrete time parameters go.")
  .method( "setRAMData", &cpptsemRAMmodel::setRAMData, "Set up the dataset for RAM models. Expects a list with sampleSize, nObservedVariables, rawData, observedMeans, observedCov, expectedSelector, expectedMeans, expectedCovariance, m2LL for each missingness pattern")
  .method( "computeRAM", &cpptsemRAMmodel::computeRAM, "Computes the RAM matrices")
  .method( "fitRAM", &cpptsemRAMmodel::fitRAM, "Fit the RAM model: outputs the -2log likelihood")
  .method( "approxRAMGradients", &cpptsemRAMmodel::approxRAMGradients, "Returns a central approximation of the gradients.")
  .method( "approxRAMGradient", &cpptsemRAMmodel::approxRAMGradient, "Returns a central approximation of a single gradient.")
  .method( "GIST", &cpptsemRAMmodel::GIST, "Optimizes with GIST")
  .method( "computePenalty", &cpptsemRAMmodel::computePenalty, "Computes penalty value")
  .method( "computeSubgradients", &cpptsemRAMmodel::computeSubgradients, "Computes subgradients")
  ;
}


