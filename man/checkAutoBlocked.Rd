% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testthat_helper.R
\name{checkAutoBlocked}
\alias{checkAutoBlocked}
\title{checkAutoBlocked}
\usage{
checkAutoBlocked(ctInit, regCtsemObject, threshold, testIC)
}
\arguments{
\item{ctInit}{init object from ctsem}

\item{regCtsemObject}{object from regCtsem}

\item{threshold}{how close to zero should differences be to be treated as zero?}

\item{testIC}{should information criteria be tested?}
}
\description{
used in testthat to check the automatic k-fold blocked cross-validation feature of regCtsem
}
\keyword{internal}
