% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{logLik_GGEO}
\alias{logLik_GGEO}
\title{logLik function for GGEO}
\usage{
logLik_GGEO(param = c(0, 0), x)
}
\arguments{
\item{param}{vector with parameters in log and logit scale.}

\item{x}{vector with the response variable.}
}
\value{
returns the loglikelihood given the parameters and random sample.
}
\description{
Calculates logLik for GGEO distribution.
}
\keyword{internal}
