% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exact_optimization.R
\name{exact_getFitIndices}
\alias{exact_getFitIndices}
\title{exact_getFitIndices}
\usage{
exact_getFitIndices(parameterLabels, fitAndParameters, lambdas, sampleSize)
}
\arguments{
\item{parameterLabels}{labels of optimized parameters}

\item{fitAndParameters}{table with fit and parameter values}

\item{lambdas}{vector of penalty values (tuning parameter). E.g., seq(0,1,.01)}

\item{sampleSize}{sample size}
}
\description{
computes fit indices for optimization = "exact"
}
\details{
NOTE: Function located in file exact_optimization.R
}
\author{
Jannik Orzek
}
\keyword{internal}
