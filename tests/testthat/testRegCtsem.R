### Check basic regCtsem features ###

library(regCtsem)

testthat::test_that(desc = "Testing basic features of regCtsem", code = {
  skip_on_cran()
  set.seed(17544)

  ## define the population model:

  # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
  ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)

  generatingModel<-ctsem::ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
                                  MANIFESTVAR=diag(0,2),
                                  LAMBDA=diag(1,2),
                                  DRIFT=ct_drift,
                                  DIFFUSION=matrix(c(.5,0,0,.5),2),
                                  CINT=matrix(c(0,0),nrow=2),
                                  T0MEANS=matrix(0,ncol=1,nrow=2),
                                  T0VAR=diag(1,2), type = "omx")

  # simulate a training data set
  traindata <- ctsem::ctGenerate(generatingModel,n.subjects = 100, wide = TRUE)
  testdata_1 <- ctsem::ctGenerate(generatingModel,n.subjects = 100, wide = TRUE)

  myModel <- ctsem::ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
                            LAMBDA=diag(1,2),
                            MANIFESTVAR=diag(0,2),
                            CINT=matrix(c(0,0),nrow=2),
                            DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
                            T0MEANS=matrix(0,ncol=1,nrow=2),
                            T0VAR=diag(1,2), type = "omx")

  # fit the model using ctsem:
  fit_myModel <- ctsemOMX::ctFit(traindata, myModel)
  fit_myModel_fortest <- ctsemOMX::ctFit(testdata_1, myModel, useOptimizer = FALSE)

  # select DRIFT values:
  regIndicators <- fit_myModel$mxobj$DRIFT$labels[!diag(T,2)]

  message("Testing GIST, GLMNET and approximate optimization with manual cross-validation.")
  regModel_BIC_GIST <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                              dataset = traindata,
                                              optimizer = "GIST",
                                              regIndicators = regIndicators,
                                              lambdasAutoLength = 10,
                                              standardizeDrift = "asymptoticDiffusion",
                                              cvSample = testdata_1))
  regModel_BIC_GLMNET <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                                dataset = traindata,
                                                optimizer = "GLMNET",
                                                regIndicators = regIndicators,
                                                lambdasAutoLength = 10,
                                                standardizeDrift = "asymptoticDiffusion",
                                                cvSample = testdata_1))
  regModel_BIC_Approx <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                                dataset = traindata,
                                                optimization = "approx",
                                                regIndicators = regIndicators,
                                                lambdasAutoLength = 10,
                                                standardizeDrift = "asymptoticDiffusion",
                                                cvSample = testdata_1))

  expect_equal(sum(round(regModel_BIC_GIST$fit - regModel_BIC_GLMNET$fit,1)),0)
  expect_equal(sum(round(regModel_BIC_GIST$parameterEstimatesRaw - regModel_BIC_GLMNET$parameterEstimatesRaw,1)),0)
  expect_equal(sum(round(regModel_BIC_GIST$parameterEstimatesRaw - regModel_BIC_Approx$parameterEstimatesRaw,1)),0)
  regModel_BIC_checkFI <- checkFI(mxObject = fit_myModel$mxobj, regCtsemObject = regModel_BIC_GIST, cvModel = fit_myModel_fortest$mxobj)
  expect_equal(regModel_BIC_checkFI, TRUE)
  regModel_BIC_checkFI <- checkFI(mxObject = fit_myModel$mxobj, regCtsemObject = regModel_BIC_Approx, cvModel = fit_myModel_fortest$mxobj)
  expect_equal(regModel_BIC_checkFI, TRUE)


  message("Testing GIST, GLMNET and approximate optimization with automatic cross-validation.")
  regModel_CV_GIST <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                             dataset = traindata,
                                             optimizer = "GIST",
                                             regIndicators = regIndicators,
                                             lambdasAutoLength = 10,
                                             autoCV = "kFold",
                                             k = 3,
                                             standardizeDrift = "asymptoticDiffusion"))
  regModel_autoCV_checkAutoCV <- checkAutoKFold(ctInit = fit_myModel$ctmodelobj, regCtsemObject = regModel_CV_GIST)
  expect_equal(regModel_autoCV_checkAutoCV, TRUE)

  regModel_CV_GLMNET <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                             dataset = traindata,
                                             optimizer = "GLMNET",
                                             regIndicators = regIndicators,
                                             lambdasAutoLength = 10,
                                             autoCV = "kFold",
                                             k = 3,
                                             standardizeDrift = "asymptoticDiffusion"))
  regModel_autoCV_checkAutoCV <- checkAutoKFold(ctInit = fit_myModel$ctmodelobj, regCtsemObject = regModel_CV_GLMNET)
  expect_equal(regModel_autoCV_checkAutoCV, TRUE)

  regModel_CV_Approx <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                               dataset = traindata,
                                               optimization = "approx",
                                               regIndicators = regIndicators,
                                               lambdasAutoLength = 10,
                                               autoCV = "kFold",
                                               k = 3,
                                               standardizeDrift = "asymptoticDiffusion"))
  regModel_autoCV_checkAutoCV <- checkAutoKFold(ctInit = fit_myModel$ctmodelobj, regCtsemObject = regModel_CV_Approx)
  expect_equal(regModel_autoCV_checkAutoCV, TRUE)

})
