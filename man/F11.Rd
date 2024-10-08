% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{F11}
\alias{F11}
\title{Auxiliar function for hyper Poisson}
\usage{
F11(z, c, maxiter_series = 10000, tol = 1e-10)
}
\arguments{
\item{z, c}{values for F11.}

\item{maxiter_series}{maximum value to obtain F11.}

\item{tol}{this is the tolerance of the infinite sum.}
}
\value{
returns the value for the F11 function.
}
\description{
This function is used inside density function of Hyper Poisson.
}
\keyword{internal}
