#' createKalmanMultiSubjectModel
#'
#' Creates an mxModel Object for N >= 1 individuals using the Kalman filter
#' @param ctsemObject object of type ctsemInit from ctModel
#' @param dataset data set in wide format compatible to ctsem
#' @param useOptimizer Boolean: should the model be optimized
#' @export
createKalmanMultiSubjectModel <- function(ctsemObject, dataset, useOptimizer, silent = FALSE, KalmanStartValues = NULL){


  suppressMessages(invisible(capture.output(fit_kalmanModels <- ctsemOMX::ctFit(ctmodelobj = ctsemObject,
                                                                             dat = dataset,
                                                                             objective = "Kalman",
                                                                             fit = useOptimizer))))
  mxObject <- fit_kalmanModels$mxobj
  if(!is.null(KalmanStartValues)){
    parameterLabels <- names(OpenMx::omxGetParameters(mxObject))
    if(!all(names(KalmanStartValues)%in%parameterLabels)){
      stop("KalmanStartValues must have the same labels as the parameters in the model.")
    }
    mxObject <- OpenMx::omxSetParameters(mxObject, labels = parameterLabels, values = KalmanStartValues[parameterLabels])
  }
  return(mxObject)

  ######## Not used #######
  # create individual models
  # Note that we assume all persons to have the same parameter values
  individualModels <- vector("list", length = nrow(dataset))
  individualModelNames <- paste0("person", 1:nrow(dataset))

  for(person in 1:nrow(dataset)){
    if(!silent){
    cat('\r',paste0("Setting up the Kalman model: ", person, " of ", nrow(dataset)))
    flush.console()}
    suppressMessages(invisible(capture.output(individualModels[[person]] <- OpenMx::mxModel(name = individualModelNames[person],
                                                  ctsemOMX::ctFit(dat = t(as.matrix(dataset[person,])),
                                                               ctmodelobj = ctsemObject,
                                                               objective = 'Kalman',
                                                               useOptimizer = FALSE)$mxobj))))
  }

  pointersToMatricesAndAlgebras <- vector("list", length = (length(individualModels[[1]]$algebras) + length(individualModels[[1]]$matrices)))
  namesOfMatricesAndAlgebras <- c(names(individualModels[[1]]$algebras), names(individualModels[[1]]$matrices))
  names(pointersToMatricesAndAlgebras) <- namesOfMatricesAndAlgebras
  for(i in 1:length(namesOfMatricesAndAlgebras)){
    pointersToMatricesAndAlgebras[[i]] <- mxAlgebraFromString(name = namesOfMatricesAndAlgebras[i], algString = paste0("person1.", namesOfMatricesAndAlgebras[i]))
  }

  mxIndividualModels <- OpenMx::mxModel(name = "MultiModel",
                                        submodels  = individualModels,
                                        # OpenMx::mxData(dataset, type = "raw"),
                                        OpenMx::mxFitFunctionMultigroup(individualModelNames),
                                        # for easy access to the matrices and algebras:
                                        pointersToMatricesAndAlgebras
  )

  mxIndividualModels <- OpenMx::omxAssignFirstParameters(mxIndividualModels)

  fit.mxIndividualModels <- mxRun(mxIndividualModels, silent = TRUE, useOptimizer = useOptimizer)
  cat("\n")
  return(fit.mxIndividualModels)
}




