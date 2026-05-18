#' New Poisson-generalised Lindley distribution
#'
#' @author Tomas Mesa, \email{tomas.mesaz@udea.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Poisson-generalised Lindley (NPGL)
#' distribution
#' with parameters \eqn{\mu} and \eqn{\sigma}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param sigma vector of the sigma parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @references
#' Altun, E.
#' A new two-parameter discrete poisson-generalized Lindley distribution with
#' properties and applications to healthcare data sets. Comput Stat 36,
#' 2841–2861 (2021). https://doi.org/10.1007/s00180-021-01097-0
#'
#' @seealso \link{NPGL}.
#'
#' @details
#' The Poisson-generalised Lindley distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has support \eqn{x = 0, 1, 2, \ldots} and probability mass function given by
#'
#' \eqn{f(x \mid \mu, \sigma)=\frac{\mu^2+\frac{\mu^{\sigma}(\mu+1)^{1-\sigma}\Gamma(x+\sigma)}{\Gamma(\sigma)\Gamma(x+1)}}{(\mu+1)^{x+2}}}
#'
#' with \eqn{\mu > 0} and \eqn{\sigma > 0}.
#'
#' This distribution is useful for modeling over-dispersed count data.
#'
#' Note: in this implementation we changed the original parameters \eqn{\theta} and \eqn{\alpha}
#' for \eqn{\mu} and \eqn{\sigma} respectively, we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dNPGL} gives the density, \code{pNPGL} gives the distribution
#' function, \code{qNPGL} gives the quantile function, \code{rNPGL}
#' generates random deviates.
#'
#' @example examples/examples_dNPGL.R
#'
#' @export
#' @useDynLib DiscreteDists
#' @importFrom Rcpp sourceCpp
dNPGL <- function(x, mu=0.1, sigma=2, log=FALSE){
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")

  # Ensure same length vector
  ly <- max(length(x), length(mu), length(sigma))
  xx <- rep(x, length=ly)
  mu <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)

  # Temporal change for invalid x's
  xx[x < 0] <- 0
  xx[is.infinite(x)] <- 1
  xx[!is.na(x) & abs(x - round(x)) > .Machine$double.eps^0.5] <- 2 # No integers

  # pdf in log-scale
  part1 <- - 2*log(mu+1) - xx*log(mu+1)
  part2a <- log(((mu^2)*gamma(sigma)*gamma(xx+1))+(mu^sigma)*((mu+1)^(1-sigma))*gamma(xx+sigma))
  part2b <- - log(gamma(sigma)) - log(gamma(xx+1))
  part2 <- part2a + part2b
  p <- part1 + part2

  # Assign values for invalid x's
  p[x < 0] <- -Inf
  p[is.infinite(x)] <- -Inf
  p[!is.na(x) & abs(x - round(x)) > .Machine$double.eps^0.5] <- -Inf

  if (log == TRUE)
    return(p)
  else
    return(exp(p))
}
#' @export
#' @rdname dNPGL
pNPGL <- function(q, mu=0.1, sigma=2, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")

  # Ensure same length vector
  ly <- max(length(q), length(mu), length(sigma))
  qq <- rep(q, length = ly)
  mu <- rep(mu, length = ly)
  sigma <- rep(sigma, length = ly)

  # Temporal change for invalid q's
  qq[q < 0] <- 0
  qq[q == Inf] <- 0

  # For non-integer q's the cumulative is the same as the lower integer
  qq <- as.integer(qq)

  # Auxiliary function
  fn <- function(q, mu, sigma) sum(dNPGL(x=0:q, mu=mu, sigma=sigma))
  Vec_fn <- Vectorize(fn)

  # The cumulative
  cdf <- Vec_fn(q=qq, mu=mu, sigma=sigma)

  # Assign values for invalid q's
  cdf[q < 0] <- 0
  cdf[q == Inf] <- 1

  if (lower.tail == FALSE)
    cdf <- 1 - cdf
  if (log.p == TRUE)
    cdf <- log(cdf)

  return(cdf)
}
#' @importFrom stats runif
#' @export
#' @rdname dNPGL
rNPGL <- function(n, mu=0.1, sigma=2){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))

  n <- ceiling(n)
  u <- runif(n=n)

  p_exp <- mu / (mu + 1)

  lambda <- ifelse(u < p_exp,
                   rexp(n, rate=mu),
                   rgamma(n, shape=sigma, rate=mu))

  x <- rpois(n=n, lambda=lambda)

  return(x)
}
#' @export
#' @rdname dNPGL
qNPGL <- function(p, mu = 0.1, sigma = 2, lower.tail = TRUE, log.p = FALSE) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")

  if(log.p == TRUE)
    p <- exp(p)
  if(lower.tail == FALSE)
    p <- 1 - p

  # Ensure same length vector
  ly <- max(length(p), length(mu), length(sigma))
  pp <- rep(p, length = ly)
  mu <- rep(mu, length = ly)
  sigma <- rep(sigma, length = ly)

  # Temporal change for invalid p's
  pp[p < 0] <- 0.1
  pp[p > 1] <- 0.1
  pp[p == 1] <- 0.1
  pp[p == 0] <- 0.1

  # Begin auxiliary function
  one_quantile_npgl <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      i <- 0
      prob <- dNPGL(x=i, mu=mu, sigma=sigma, log=FALSE)
      F <- prob
      while (p >= F) {
        i <- i + 1
        prob <- dNPGL(x=i, mu=mu, sigma, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_npgl <- Vectorize(one_quantile_npgl)
  # End auxiliary function

  # The quantile
  q <- one_quantile_npgl(p=pp, mu=mu, sigma=sigma)

  # To deal with invalid p's
  q[p < 0] <- NaN
  q[p > 1] <- NaN
  q[p == 1] <- Inf
  q[p == 0] <- 0

  return(q)
}
