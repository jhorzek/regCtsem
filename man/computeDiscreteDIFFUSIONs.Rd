% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{computeDiscreteDIFFUSIONs}
\alias{computeDiscreteDIFFUSIONs}
\title{computeDiscreteDIFFUSIONs}
\usage{
computeDiscreteDIFFUSIONs(
  DRIFTHASHInverse,
  DIFFUSION,
  DRIFTHASHExponentialUnique,
  discreteDIFFUSIONUnique
)
}
\arguments{
\item{DRIFTHASHInverse}{inverse of drift matrix}

\item{DIFFUSION}{matrix with diffusion values}

\item{DRIFTHASHExponentialUnique}{list with expm(drifthash)-names, dTs and results}

\item{discreteDIFFUSIONUnique}{list with diffusion-names, dTs and results}
}
\description{
Computes the discrete time diffusion matrices given the DRIFTHASHInverse, DIFFUSION,
DRIFTHASHExponentialUnique (list with expm(drifthash)-names, dTs and results), and discreteDIFFUSIONUnique (list with diffusion-names, dTs and results)
The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
}
\keyword{internal}
