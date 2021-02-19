prepareKalmanMatrices <- function(nlatent, nmanifest, Tpoints, sampleSize){
  latentScores <- matrix(-9999999, nrow = sampleSize, ncol = nlatent * Tpoints)
  predictedManifestValues <- matrix(-9999999, nrow = sampleSize, ncol = nmanifest * Tpoints)

  return(list(
    "latentScores" = latentScores,
    "predictedManifestValues" = predictedManifestValues
  ))
}

