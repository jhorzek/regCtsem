test_that(desc = "Testing Kalman filter", code = {
  library(regCtsem)
  set.seed(123)

  N <- 20
  ## define the population model:

  # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
  ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)

  generatingModel<-ctsem::ctModel(Tpoints=100,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
                                  MANIFESTVAR=diag(0,2),
                                  LAMBDA=diag(1,2),
                                  DRIFT=ct_drift,
                                  DIFFUSION=matrix(c(.5,0,0,.5),2),
                                  CINT=matrix(c(0,0),nrow=2),
                                  T0MEANS=matrix(0,ncol=1,nrow=2),
                                  T0VAR=diag(1,2), type = "omx")

  # simulate a training data set
  traindata <- ctsem::ctGenerate(generatingModel,n.subjects = N, wide = TRUE)
  testdata_1 <- ctsem::ctGenerate(generatingModel,n.subjects = N, wide = TRUE)

  myModel <- ctsem::ctModel(Tpoints=100,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
                            LAMBDA=diag(1,2),
                            MANIFESTVAR=diag(0,2),
                            CINT=matrix(c(0,0),nrow=2),
                            DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
                            T0MEANS=matrix(0,ncol=1,nrow=2),
                            T0VAR=diag(1,2), type = "omx")

  # fit the model using ctsem:
  fit_myModel <- suppressWarnings(ctsemOMX::ctFit(traindata, myModel, objective = "Kalman"))
  fit_myModel_fortest <- suppressWarnings(ctsemOMX::ctFit(testdata_1, myModel, useOptimizer = FALSE, objective = "Kalman"))

  # select DRIFT values:
  regIndicators <- fit_myModel$mxobj$DRIFT$labels[!diag(T,2)]

  # test optimizers without penalty:

  regModel <- try(regCtsem::regCtsem(ctsemObject = fit_myModel_fortest,
                                     dataset = traindata,
                                     optimizer = "GIST",
                                     regIndicators = regIndicators,
                                     lambdas = 0,
                                     standardizeDrift = "No",
                                     cvSample = testdata_1,
                                     control = controlGIST(numStart = 0,
                                                           forceCpptsem = TRUE,
                                                           approxFirst = FALSE)))
  testthat::expect_equal(round(regModel$fit["m2LL",] - fit_myModel$mxobj$fitfunction$result[[1]],4)[[1]], 0)
  testthat::expect_equal(round(regModel$fit["AIC",] - AIC(fit_myModel$mxobj),4)[[1]], 0)
  testthat::expect_equal(round(regModel$fit["BIC",] - (fit_myModel$mxobj$fitfunction$result[[1]] + log(N)*length(omxGetParameters(fit_myModel$mxobj))),4)[[1]], 0)
  testthat::expect_equal(
    all(
      abs(
        regModel$parameterEstimatesRaw[,1] - omxGetParameters(fit_myModel$mxobj)[rownames(regModel$parameterEstimatesRaw)]
      ) < .001), TRUE)


  regModel <- try(regCtsem::regCtsem(ctsemObject = fit_myModel_fortest,
                                     dataset = traindata,
                                     optimizer = "GLMNET",
                                     regIndicators = regIndicators,
                                     lambdas = 0,
                                     standardizeDrift = "No",
                                     cvSample = testdata_1,
                                     control = controlGLMNET(numStart = 0,
                                                             forceCpptsem = TRUE,
                                                             approxFirst = T)))
  testthat::expect_equal(round(regModel$fit["m2LL",] - fit_myModel$mxobj$fitfunction$result[[1]],4)[[1]], 0)
  testthat::expect_equal(round(regModel$fit["AIC",] - AIC(fit_myModel$mxobj),4)[[1]], 0)
  testthat::expect_equal(round(regModel$fit["BIC",] - (fit_myModel$mxobj$fitfunction$result[[1]] + log(N)*length(omxGetParameters(fit_myModel_fortest$mxobj))),4)[[1]], 0)

  testthat::expect_equal(
    all(
      abs(
        regModel$parameterEstimatesRaw[,1] - omxGetParameters(fit_myModel$mxobj)[rownames(regModel$parameterEstimatesRaw)]
      ) < .001), TRUE)

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

  testthat::expect_equal(sum(round(regModel_BIC_GIST$fit["regM2LL",] - regModel_BIC_GLMNET$fit["regM2LL",],3)),0)
  testthat::expect_equal(sum(round(regModel_BIC_GIST$parameterEstimatesRaw - regModel_BIC_GLMNET$parameterEstimatesRaw,3)),0)
  testthat::expect_equal(sum(round(regModel_BIC_GIST$parameterEstimatesRaw - regModel_BIC_Approx$parameterEstimatesRaw,2)),0)
  regModel_BIC_checkFI <- regCtsem:::checkFI(mxObject = fit_myModel$mxobj,
                                             regCtsemObject = regModel_BIC_GIST,
                                             cvModel = fit_myModel_fortest$mxobj,
                                             threshold = 1e-4,
                                             testIC = TRUE)
  testthat::expect_equal(regModel_BIC_checkFI, TRUE)
  regModel_BIC_checkFI <- regCtsem:::checkFI(mxObject = fit_myModel$mxobj,
                                             regCtsemObject = regModel_BIC_GLMNET,
                                             cvModel = fit_myModel_fortest$mxobj,
                                             threshold = 1e-4,
                                             testIC = TRUE)
  testthat::expect_equal(regModel_BIC_checkFI, TRUE)
  regModel_BIC_checkFI <- regCtsem:::checkFI(mxObject = fit_myModel$mxobj,
                                             regCtsemObject = regModel_BIC_Approx,
                                             cvModel = fit_myModel_fortest$mxobj,
                                             threshold = 1e-4,
                                             testIC = FALSE)
  testthat::expect_equal(regModel_BIC_checkFI, TRUE)


  message("Testing GIST, GLMNET and approximate optimization with automatic cross-validation.")
  regModel_CV_GIST <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                             dataset = traindata,
                                             optimizer = "GIST",
                                             regIndicators = regIndicators,
                                             lambdasAutoLength = 10,
                                             autoCV = "kFold",
                                             k = 3,
                                             standardizeDrift = "asymptoticDiffusion"))
  regModel_autoCV_checkAutoCV <- regCtsem:::checkAutoKFold(ctInit = fit_myModel$ctmodelobj,
                                                           regCtsemObject = regModel_CV_GIST,
                                                           threshold = 1e-4,
                                                           testIC = TRUE)
  testthat::expect_equal(regModel_autoCV_checkAutoCV, TRUE)

  regModel_CV_GLMNET <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                               dataset = traindata,
                                               optimizer = "GLMNET",
                                               regIndicators = regIndicators,
                                               lambdasAutoLength = 10,
                                               autoCV = "kFold",
                                               k = 3,
                                               standardizeDrift = "asymptoticDiffusion"))
  regModel_autoCV_checkAutoCV <- regCtsem:::checkAutoKFold(ctInit = fit_myModel$ctmodelobj,
                                                           regCtsemObject = regModel_CV_GLMNET,
                                                           threshold = 1e-4,
                                                           testIC = TRUE)
  testthat::expect_equal(regModel_autoCV_checkAutoCV, TRUE)

  regModel_CV_Approx <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                               dataset = traindata,
                                               optimization = "approx",
                                               regIndicators = regIndicators,
                                               lambdasAutoLength = 10,
                                               autoCV = "kFold",
                                               k = 3,
                                               standardizeDrift = "asymptoticDiffusion"))
  regModel_autoCV_checkAutoCV <- regCtsem:::checkAutoKFold(ctInit = fit_myModel$ctmodelobj,
                                                regCtsemObject = regModel_CV_Approx,
                                                threshold = 1e-4,
                                                testIC = FALSE)
  testthat::expect_equal(regModel_autoCV_checkAutoCV, TRUE)

  message("Testing GIST, GLMNET and approximate optimization with blocked cross-validation.")
  regModel_CV_GIST <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                             dataset = traindata,
                                             optimizer = "GIST",
                                             regIndicators = regIndicators,
                                             lambdasAutoLength = 10,
                                             autoCV = "Blocked",
                                             k = 3,
                                             standardizeDrift = "asymptoticDiffusion"))
  regModel_autoCV_checkAutoCV <- regCtsem:::checkAutoBlocked(ctInit = fit_myModel$ctmodelobj,
                                                             regCtsemObject = regModel_CV_GIST,
                                                             threshold = 1e-4,
                                                             testIC = TRUE)
  testthat::expect_equal(regModel_autoCV_checkAutoCV, TRUE)

  regModel_CV_GLMNET <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                               dataset = traindata,
                                               optimizer = "GLMNET",
                                               regIndicators = regIndicators,
                                               lambdasAutoLength = 10,
                                               autoCV = "kFold",
                                               k = 3,
                                               standardizeDrift = "asymptoticDiffusion"))
  regModel_autoCV_checkAutoCV <- regCtsem:::checkAutoKFold(ctInit = fit_myModel$ctmodelobj,
                                                           regCtsemObject = regModel_CV_GLMNET,
                                                           threshold = 1e-4,
                                                           testIC = TRUE)
  testthat::expect_equal(regModel_autoCV_checkAutoCV, TRUE)

  regModel_CV_Approx <- try(regCtsem::regCtsem(ctsemObject = fit_myModel,
                                               dataset = traindata,
                                               optimization = "approx",
                                               regIndicators = regIndicators,
                                               lambdasAutoLength = 10,
                                               autoCV = "kFold",
                                               k = 3,
                                               standardizeDrift = "asymptoticDiffusion"))
  regModel_autoCV_checkAutoCV <- regCtsem:::checkAutoKFold(ctInit = fit_myModel$ctmodelobj,
                                                regCtsemObject = regModel_CV_Approx,
                                                threshold = 1e-4,
                                                testIC = FALSE)
  testthat::expect_equal(regModel_autoCV_checkAutoCV, TRUE)

})
