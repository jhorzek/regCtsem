### Check cpptsem Implementation ###

library(regCtsem)

testthat::test_that(desc = "Testing implementation of cpptsem", code = {
  for(addCINT in c(TRUE,FALSE)){
    if(addCINT){
      CINT = matrix(c("cint1", "cint2"), nrow = 2, ncol = 1)
    }else{
      CINT = matrix(c(0,0), nrow = 2, ncol = 1)
    }
    stationary <- c('T0TRAITEFFECT','T0TIPREDEFFECT')

    ## ctsem model without trait
    AnomAuthmodel1 <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
                              Tpoints = 5, n.latent = 2, n.manifest = 2,
                              MANIFESTVAR=diag(0, 2),
                              TRAITVAR = NULL,
                              CINT = CINT)
    AnomAuthfit1 <- ctFit(AnomAuth, AnomAuthmodel1, useOptimizer = FALSE, stationary = stationary)
    gradientModel1 <- OpenMx::mxRun(OpenMx::mxModel(AnomAuthfit1$mxobj,
                                                    OpenMx::mxComputeSequence(steps=list(
                                                      OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
                                                                                    hessian = FALSE))
                                                    )))
    centralGrandients <- gradientModel1$compute$steps[[1]]$output[["gradient"]][,"central"]
    names(centralGrandients) <- rownames(gradientModel1$compute$steps[[1]]$output[["gradient"]])

    ## with cpptsem
    cpptsemmodel1 <- cpptsemFromCtsem(ctsemModel = AnomAuthfit1, wideData = AnomAuth)
    cpptsemmodel1$computeRAM()
    cpptsemmodel1$fitRAM()

    # check fit and gradients
    expect_equal(round(cpptsemmodel1$m2LL - AnomAuthfit1$mxobj$fitfunction$result[[1]],2) , 0)
    expect_equal(sum(round(cpptsemmodel1$approxRAMGradients((1.1 * 10^(-16))^(1/3))[names(centralGrandients)] - centralGrandients,2)), 0)

    # change parameter values
    AnomAuthfit1_1 <- ctFit(AnomAuth, AnomAuthmodel1, useOptimizer = TRUE, stationary = stationary)
    newParameters <- omxGetParameters(AnomAuthfit1_1$mxobj)
    cpptsemmodel1$setParameterValues(newParameters, names(newParameters))
    cpptsemmodel1$computeRAM()
    cpptsemmodel1$fitRAM()
    expect_equal(round(cpptsemmodel1$m2LL - AnomAuthfit1_1$mxobj$fitfunction$result[[1]],2), 0)

    ## ctsem model with trait
    AnomAuthmodel2 <- ctModel(LAMBDA = matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2),
                              Tpoints = 5, n.latent = 2, n.manifest = 2,
                              MANIFESTVAR=diag(0, 2),
                              TRAITVAR = "auto",
                              CINT = CINT)
    AnomAuthfit2 <- ctFit(AnomAuth, AnomAuthmodel2, useOptimizer = FALSE, stationary = stationary)
    AnomAuthfit2$mxobj$fitfunction$result[[1]]
    gradientModel2 <- OpenMx::mxRun(OpenMx::mxModel(AnomAuthfit2$mxobj,
                                                    OpenMx::mxComputeSequence(steps=list(
                                                      OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
                                                                                    hessian = FALSE))
                                                    )))
    centralGrandients <- gradientModel2$compute$steps[[1]]$output[["gradient"]][,"central"]
    names(centralGrandients) <- rownames(gradientModel2$compute$steps[[1]]$output[["gradient"]])

    ## with cpptsem
    cpptsemmodel2 <- cpptsemFromCtsem(AnomAuthfit2, wideData = AnomAuth)
    cpptsemmodel2$computeRAM()
    cpptsemmodel2$fitRAM()

    # check fit and gradients
    expect_equal(round(cpptsemmodel2$m2LL - AnomAuthfit2$mxobj$fitfunction$result[[1]],2) , 0)
    expect_equal(sum(round(cpptsemmodel2$approxRAMGradients((1.1 * 10^(-16))^(1/3))[names(centralGrandients)] - centralGrandients,2)), 0)

    ## Example 3: Kalman Filter
    set.seed(175446)
    ## define the population model:

    # set the drift matrix. Note that drift eta_1_eta2 is set to equal 0 in the population.
    ct_drift <- matrix(c(-.3,.2,0,-.5), ncol = 2)

    generatingModel<-ctsem::ctModel(Tpoints=20,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
                                    MANIFESTVAR=diag(0,2),
                                    LAMBDA=diag(1,2),
                                    DRIFT=ct_drift,
                                    DIFFUSION=matrix(c(.5,0,0,.5),2),
                                    CINT=matrix(0,nrow = 2, ncol = 1),
                                    T0MEANS=matrix(0,ncol=1,nrow=2),
                                    T0VAR=diag(1,2), type = "omx")

    # simulate a training data and testing data set
    traindata <- ctsem::ctGenerate(generatingModel,n.subjects = 20, wide = TRUE)

    ## Build the analysis model. Note that drift eta1_eta2 is freely estimated
    # although it is 0 in the population.
    myModel <- ctsem::ctModel(Tpoints=20,n.latent=2,n.TDpred=0,n.TIpred=0,n.manifest=2,
                              LAMBDA=diag(1,2),
                              MANIFESTVAR=diag(0,2),
                              CINT=CINT,
                              DIFFUSION=matrix(c('eta1_eta1',0,0,'eta2_eta2'),2),
                              T0MEANS=matrix(0,ncol=1,nrow=2),
                              T0VAR="auto", type = "omx")
    myModel <- ctFit(myModel, dat = traindata, objective = "Kalman", useOptimizer = F, stationary = stationary)

    gradientModel3 <- OpenMx::mxRun(OpenMx::mxModel(myModel$mxobj,
                                                    OpenMx::mxComputeSequence(steps=list(
                                                      OpenMx::mxComputeNumericDeriv(checkGradient = FALSE,
                                                                                    hessian = FALSE))
                                                    )))
    centralGrandients <- gradientModel3$compute$steps[[1]]$output[["gradient"]][,"central"]
    names(centralGrandients) <- rownames(gradientModel3$compute$steps[[1]]$output[["gradient"]])

    KalmanScores <- mxKalmanScores(myModel$mxobj)


    ## with cpptsem
    cpptsemmodel3 <- cpptsemFromCtsem(ctsemModel = myModel,wideData = traindata)
    cpptsemmodel3$computeAndFitKalman()


    # check fit and gradients
    expect_equal(round(cpptsemmodel3$m2LL - myModel$mxobj$fitfunction$result[[1]],2) , 0)
    expect_equal(sum(round(cpptsemmodel3$approxKalmanGradients((1.1 * 10^(-16))^(1/3))[names(centralGrandients)] - centralGrandients,2)), 0)
    expect_equal(sum(round(KalmanScores$xUpdated[2:21,] - matrix(cpptsemmodel3$latentScores[1,], ncol = 2, byrow = TRUE),2)), 0)

  }
})



