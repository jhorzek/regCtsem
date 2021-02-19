prepareFMatrix <- function(nlatent, nmanifest, Tpoints){

  Fmatrix <- matrix(0,
                    nrow = nmanifest*Tpoints,
                    ncol = Tpoints*nlatent + nlatent + Tpoints*nmanifest)

  diag(Fmatrix[,(Tpoints*nlatent + nlatent + 1):(Tpoints*nlatent + nlatent + Tpoints*nmanifest)]) <- 1

  return(Fmatrix)
}

