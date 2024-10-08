% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DMOLBE.R
\name{DMOLBE}
\alias{DMOLBE}
\title{The DMOLBE family}
\usage{
DMOLBE(mu.link = "log", sigma.link = "log")
}
\arguments{
\item{mu.link}{defines the mu.link, with "log" link as the default for the mu parameter.}

\item{sigma.link}{defines the sigma.link, with "log" link as the default for the sigma.}
}
\value{
Returns a \code{gamlss.family} object which can be used
to fit a DMOLBE distribution
in the \code{gamlss()} function.
}
\description{
The function \code{DMOLBE()} defines the Discrete Marshall-Olkin Length Biased
Exponential distribution, a two parameter
distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
using the function \code{gamlss()}.
}
\details{
The DMOLBE distribution with parameters \eqn{\mu} and \eqn{\sigma}
has a support 0, 1, 2, ... and mass function given by

\eqn{f(x | \mu, \sigma) = \frac{\sigma ((1+x/\mu)\exp(-x/\mu)-(1+(x+1)/\mu)\exp(-(x+1)/\mu))}{(1-(1-\sigma)(1+x/\mu)\exp(-x/\mu)) ((1-(1-\sigma)(1+(x+1)/\mu)\exp(-(x+1)/\mu))}}

with \eqn{\mu > 0} and \eqn{\sigma > 0}
}
\examples{
# Example 1
# Generating some random values with
# known mu and sigma
set.seed(1234)
y <- rDMOLBE(n=100, mu=10, sigma=7)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=DMOLBE,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what='mu'))
exp(coef(mod1, what='sigma'))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ DMOLBE
gendat <- function(n) {
  x1 <- runif(n, min=0.4, max=0.6)
  x2 <- runif(n, min=0.4, max=0.6)
  mu    <- exp(1.21 - 3 * x1) # 0.75 approximately
  sigma <- exp(1.26 - 2 * x2) # 1.30 approximately
  y <- rDMOLBE(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1,x2=x2)
}

set.seed(123)
dat <- gendat(n=350)

# Fitting the model
mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=DMOLBE, data=dat,
                 control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod2)
}
\references{
\insertRef{Aljohani2023}{DiscreteDists}
}
\seealso{
\link{dDMOLBE}.
}
\author{
Olga Usuga, \email{olga.usuga@udea.edu.co}
}
