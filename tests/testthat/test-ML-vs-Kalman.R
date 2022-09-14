test_that("ML vs Kalman", {
  library(regCtsem)
  set.seed(17544)

  ## define the population model:

  # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
  ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)

  generatingModel<-ctsem::ctModel(Tpoints=15,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
                                  MANIFESTVAR=diag(0,2),
                                  LAMBDA=diag(1,2),
                                  DRIFT=ct_drift,
                                  DIFFUSION=matrix(c(.5,0,0,.5),2),
                                  CINT=matrix(c(0,0),nrow=2),
                                  T0MEANS=matrix(0,ncol=1,nrow=2),
                                  T0VAR=diag(1,2), type = "omx")

  # simulate a training data set
  traindata <- ctsem::ctGenerate(generatingModel,n.subjects = 30, wide = TRUE)

  myModel <- ctsem::ctModel(Tpoints=15,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
                            LAMBDA=diag(1,2),
                            MANIFESTVAR=diag(0,2),
                            CINT=matrix(c(0,0),nrow=2),
                            DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
                            T0MEANS=matrix(0,ncol=1,nrow=2),
                            T0VAR=diag(1,2),
                            type = "omx")

  # fit the model using ctsem:
  fit_myModel_ML <- ctsemOMX::ctFit(traindata, myModel, objective = "mxRAM")

  # select DRIFT values:
  regIndicators <- fit_myModel_ML$mxobj$DRIFT$labels[!diag(T,2)]

  regModel_ML <- regCtsem::regCtsem(ctsemObject = fit_myModel_ML,
                                        dataset = traindata,
                                        regIndicators = regIndicators)

  # and with Kalman filter
  fit_myModel_Kalman <- ctsemOMX::ctFit(traindata, myModel, objective = "Kalman")

  regModel_Kalman <- regCtsem::regCtsem(ctsemObject = fit_myModel_Kalman,
                                        dataset = traindata,
                                        regIndicators = regIndicators)

  testthat::expect_equal(all(round(regModel_ML$parameters - regModel_Kalman$parameters,3) == 0),
                         TRUE)


  # with standardization
  ## With asymptoticDiffusion
  regModel_ML <- regCtsem::regCtsem(ctsemObject = fit_myModel_ML,
                                    dataset = traindata,
                                    regIndicators = regIndicators,
                                    standardizeDrift = "asymptoticDiffusion")

  regModel_Kalman <- regCtsem::regCtsem(ctsemObject = fit_myModel_Kalman,
                                        dataset = traindata,
                                        regIndicators = regIndicators,
                                        standardizeDrift = "asymptoticDiffusion")

  testthat::expect_equal(all(round(regModel_ML$parameters - regModel_Kalman$parameters,3) == 0),
                         TRUE)

  ## With T0VAR
  regModel_ML <- regCtsem::regCtsem(ctsemObject = fit_myModel_ML,
                                    dataset = traindata,
                                    regIndicators = regIndicators,
                                    standardizeDrift = "T0VAR")

  regModel_Kalman <- regCtsem::regCtsem(ctsemObject = fit_myModel_Kalman,
                                        dataset = traindata,
                                        regIndicators = regIndicators,
                                        standardizeDrift = "T0VAR")

  testthat::expect_equal(all(round(regModel_ML$parameters - regModel_Kalman$parameters,3) == 0),
                         TRUE)

})
