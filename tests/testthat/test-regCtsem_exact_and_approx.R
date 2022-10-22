test_that("exact and approx", {
  library(regCtsem)
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

  myModel <- ctsem::ctModel(Tpoints=10,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
                            LAMBDA=diag(1,2),
                            MANIFESTVAR=diag(0,2),
                            CINT=matrix(c(0,0),nrow=2),
                            DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
                            T0MEANS=matrix(0,ncol=1,nrow=2),
                            T0VAR=diag(1,2), type = "omx")

  # fit the model using ctsem:
  fit_myModel <- suppressWarnings(ctsemOMX::ctFit(traindata, myModel))

  # select DRIFT values:
  regIndicators <- fit_myModel$mxobj$DRIFT$labels[!diag(T,2)]

  # test optimizers without penalty:

  regModel_GIST <- regCtsem::regCtsem(ctsemObject = fit_myModel,
                                     dataset = traindata,
                                     regIndicators = regIndicators,
                                     optimizer = "GIST"
                                     )
  testthat::expect_equal(all(regModel_GIST$parameterEstimatesRaw != 0),
                         FALSE)

  regModel_GLMNET <- regCtsem::regCtsem(ctsemObject = fit_myModel,
                                 dataset = traindata,
                                 regIndicators = regIndicators,
                                 optimizer = "GLMNET"
  )

  testthat::expect_equal(all(regModel_GLMNET$parameterEstimatesRaw != 0),
                         FALSE)

  testthat::expect_equal(all(round(regModel_GIST$parameters - regModel_GLMNET$parameters,3) == 0),
                         TRUE)

  regModel_Rsolnp <- regCtsem::regCtsem(ctsemObject = fit_myModel,
                                      dataset = traindata,
                                      regIndicators = regIndicators,
                                      optimization = "approx"
  )
  testthat::expect_equal(all(round(regModel_GIST$parameters - regModel_Rsolnp$parameters,2) == 0),
                         TRUE)

  testthat::expect_equal(all(regModel_Rsolnp$parameterEstimatesRaw != 0),
                         TRUE)
})
