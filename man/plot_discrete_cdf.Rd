% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_discrete_cdf.R
\name{plot_discrete_cdf}
\alias{plot_discrete_cdf}
\title{Draw the CDF for a discrete random variable}
\usage{
plot_discrete_cdf(x, fx, col = "blue", lwd = 3, ...)
}
\arguments{
\item{x}{vector with the values of the random variable \eqn{X}.}

\item{fx}{vector with the probabilities of \eqn{X}.}

\item{col}{color for the line.}

\item{lwd}{line width.}

\item{...}{further arguments and graphical parameters.}
}
\value{
A plot with the cumulative distribution function.
}
\description{
Draw the CDF for a discrete random variable
}
\examples{
# Example 1
# for a particular distribution

x <- 1:6
fx <- c(0.19, 0.21, 0.4, 0.12, 0.05, 0.03)
plot_discrete_cdf(x, fx, las=1, main="")

# Example 2
# for a Poisson distribution
x <- 0:10
fx <- dpois(x, lambda=3)
plot_discrete_cdf(x, fx, las=1,
                  main="CDF for Poisson")
}
\author{
Freddy Hernandez, \email{fhernanb@unal.edu.co}
}
