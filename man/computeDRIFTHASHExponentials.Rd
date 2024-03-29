% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{computeDRIFTHASHExponentials}
\alias{computeDRIFTHASHExponentials}
\title{computeDRIFTHASHExponentials}
\usage{
computeDRIFTHASHExponentials(DRIFTHASH, DRIFTHASHExponentialUnique)
}
\arguments{
\item{DRIFTHASH}{matrix with drift hash}

\item{DRIFTHASHExponentialUnique}{list with discrete trait names, dTs and results}
}
\value{
list with drift hash exponentials
}
\description{
Computes the DRIFTHASHExponential for the discrete time diffusion given the drifthash and DRIFTHASHExponential (list with discrete trait names, dTs and results)
The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
}
\keyword{internal}
