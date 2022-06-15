### Check cpptsem Implementation ###

testthat::test_that(desc = "Testing standardization", code = {
  skip_on_cran()
  set.seed(17046)

  library(regCtsem)

  #### Example 1 ####
  ## Regularization with FIML objective function

  ## Population model:

  # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
  ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)

  generatingModel<-ctsem::ctModel(Tpoints=10,n.latent=2,n.TDpred=0,
                                  n.TIpred=0,n.manifest=2,
                                  MANIFESTVAR=diag(0,2),
                                  LAMBDA=diag(1,2),
                                  DRIFT=ct_drift,
                                  DIFFUSION=matrix(c(.5,0,0,.5),2),
                                  CINT=matrix(c(0,0),nrow=2),
                                  T0MEANS=matrix(0,ncol=1,nrow=2),
                                  T0VAR=diag(1,2), type = "omx")

  # simulate a training data set
  dat <- ctsem::ctGenerate(generatingModel, n.subjects = 100, wide = TRUE)

  ## Build the analysis model. Note that drift eta1_eta2 is freely estimated
  # although it is 0 in the population.
  myModel <- ctsem::ctModel(Tpoints=10,n.latent=2,n.TDpred=0,
                            n.TIpred=0,n.manifest=2,
                            LAMBDA=diag(1,2),
                            MANIFESTVAR=diag(0,2),
                            CINT=matrix(c(0,0),nrow=2),
                            DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
                            T0MEANS=matrix(0,ncol=1,nrow=2),
                            T0VAR="auto", type = "omx")

  # fit the model using ctsemOMX:
  fit_myModel <- ctsemOMX::ctFit(dat, myModel)

  # select DRIFT values for regularization:
  # Note: If you are unsure what the parameters are called in
  # your model, check: showParameters(fit_myModel)
  showParameters(fit_myModel)

  # regularize the cross-effects:
  regIndicators <- c("drift_eta2_eta1", "drift_eta1_eta2")

  # Optimize model using GIST with lasso penalty
  regModel <- regCtsem::regCtsem(ctsemObject = fit_myModel,
                                 dataset = dat,
                                 regIndicators = regIndicators,
                                 lambdas = "auto",
                                 lambdasAutoLength = 20)


  testthat::expect_equal(all(regModel$setup$adaptiveLassoWeights == 1) , TRUE)


  regModel <- regCtsem::regCtsem(ctsemObject = fit_myModel,
                                 dataset = dat,
                                 regIndicators = regIndicators,
                                 standardizeDrift = "asymptoticDiffusion",
                                 lambdas = "auto",
                                 lambdasAutoLength = 20)
  R <- diag(sqrt(diag(fit_myModel$mxobj$asymDIFFUSION$values)))
  DRIFT <- fit_myModel$mxobj$DRIFT$values
  DRIFTLabels <- fit_myModel$mxobj$DRIFT$labels

  weightedDrifts <- c()
  for(drl in DRIFTLabels){
    weightedDrifts <- c(weightedDrifts,
                        DRIFT[DRIFTLabels == drl] * regModel$setup$adaptiveLassoWeights[drl]
    )
  }
  testthat::expect_equal(all(abs(matrix(weightedDrifts,2,2) - solve(R)%*%DRIFT%*%R) < 1e-8), TRUE)

  regModel <- regCtsem::regCtsem(ctsemObject = fit_myModel,
                                 dataset = dat,
                                 regIndicators = regIndicators,
                                 standardizeDrift = "T0VAR",
                                 lambdas = "auto",
                                 lambdasAutoLength = 20)
  R <- diag(sqrt(diag(fit_myModel$mxobj$T0VAR$result)))
  DRIFT <- fit_myModel$mxobj$DRIFT$values
  DRIFTLabels <- fit_myModel$mxobj$DRIFT$labels

  weightedDrifts <- c()
  for(drl in DRIFTLabels){
    weightedDrifts <- c(weightedDrifts,
                        DRIFT[DRIFTLabels == drl] * regModel$setup$adaptiveLassoWeights[drl]
    )
  }
  testthat::expect_equal(all(abs(matrix(weightedDrifts,2,2) - solve(R)%*%DRIFT%*%R) < 1e-8), TRUE)

}
)



