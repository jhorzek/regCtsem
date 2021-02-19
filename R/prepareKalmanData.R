prepareKalmanData <- function(dataset, nlatent, nmanifest, dtIndicators, Tpoints){

  sampleSize <- nrow(dataset)

  dataList <- vector("list", Tpoints)

  return(list("dataset" = dataset,
              "sampleSize" = sampleSize,
              "nlatent" = nlatent,
              "nmanifest" = nmanifest,
              "dtIndicators" = dtIndicators,
              "Tpoints" = Tpoints))
}

