#' The Poisson–transmuted record type exponential distribution
#'
#' @author Rebeca Isabel Rodriguez Gonzalez, \email{rebeca.rodriguez@udea.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the
#' Poisson–transmuted record type exponential (PTRTE) distribution
#' with parameters \eqn{\mu} and \eqn{\sigma}.
#'
#' This distribution was proposed by Erbayram and Akdogan (2025)
#' as a new discrete distribution obtained from a mixed Poisson model,
#' where the Poisson parameter follows a transmuted record type exponential
#' distribution.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param sigma vector of the sigma parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @references
#' Erbayram, T., & Akdogan, Y. (2025).
#' A new discrete model generated from mixed Poisson transmuted
#' record type exponential distribution.
#' Ricerca di Matematica, 74, 1225–1247.
#'
#' @seealso \link{PTRTE}.
#'
#' @details
#' The Poisson–transmuted record type exponential distribution with
#' parameters \eqn{\mu} and \eqn{\sigma} has support
#' \eqn{x = 0,1,2,\dots} and probability mass function given by
#'
#' \deqn{f(x | \mu, \sigma) = \frac{\mu(\sigma x\mu + 1 + \mu - \sigma)}{(1+\mu)^{x+2}}}
#'
#' with \eqn{\mu > 0} and \eqn{0 < \sigma < 1}.
#'
#' @return
#' \code{dPTRTE} gives the density,
#' \code{pPTRTE} gives the distribution function,
#' \code{qPTRTE} gives the quantile function,
#' \code{rPTRTE} generates random deviates.
#'
#' @example examples/examples_dPTRTE.R
#'
#' @export
#'
dPTRTE <- function(x, mu, sigma, log=FALSE) {
  if (any(mu <= 0)) stop("parameter mu has to be positive!")
  if (any(sigma <= 0 | sigma >= 1)) stop("parameter sigma must be in (0,1)!")

  # Asegurar vectores del mismo largo
  ly    <- max(length(x), length(mu), length(sigma))
  xx    <- rep(x, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)

  # Validaciones temporales
  xx[x < 0] <- 0
  xx[is.infinite(x)] <- 1
  xx[!is.na(x) & abs(x - round(x)) > .Machine$double.eps^0.5] <- 2 # No enteros

  # pdf en escala log
  p <- log(mu) + log(sigma * xx * mu + 1 + mu - sigma) - (xx + 2) * log(1 + mu)

  # Asignar -Inf a valores inválidos
  p[x < 0] <- -Inf
  p[is.infinite(x)] <- -Inf
  p[!is.na(x) & abs(x - round(x)) > .Machine$double.eps^0.5] <- -Inf

  if (log == FALSE)
    p <- exp(p)
  return(p)
}
#' @export
#' @rdname dPTRTE
pPTRTE <- function(q, mu, sigma, lower.tail=TRUE, log.p=FALSE) {
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(sigma <= 0 | sigma >= 1)) stop("parameter sigma must be in (0,1)!")

  # Ensure same length vector
  ly    <- max(length(q), length(mu), length(sigma))
  qq    <- rep(q, length=ly)
  mu    <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)

  # Temporal change for invalid x's
  qq[q < 0] <- 0
  qq[q == Inf] <- 0

  # For non-integer x's, the cumulative is the same as the lower integer
  qq <- as.integer(qq)

  # The cumulative
  cdf <- 1 - (1 + mu + sigma * mu * (qq + 1)) / ( (1 + mu)^(qq + 2) )


  # Assign values for invalid x's
  cdf[q < 0] <- 0
  cdf[q == Inf] <- 1

  if (lower.tail == FALSE)
    cdf <- 1 - cdf
  if (log.p == TRUE)
    cdf <- log(cdf)
  return(cdf)
}
#' @export
#' @rdname dPTRTE
qPTRTE <- function(p, mu, sigma, lower.tail=TRUE, log.p=FALSE) {
  if (any(mu <= 0)) stop("parameter mu has to be positive!")
  if (any(sigma <= 0 | sigma >= 1)) stop("parameter sigma must be in (0,1)!")

  if (log.p == TRUE)
    p <- exp(p)
  if (lower.tail == FALSE)
    p <- 1 - p

  # Ensure same length vector
  ly <- max(length(p), length(mu), length(sigma))
  pp <- rep(p, length=ly)
  mu <- rep(mu, length=ly)
  sigma <- rep(sigma, length=ly)

  # Temporal change for invalid p's
  pp[p < 0]  <-  0.1
  pp[p > 1]  <-  0.1
  pp[p == 1] <-  0.1
  pp[p == 0] <-  0.1

  # Begin auxiliar function
  one_quantile_dPTRTE <- function(p, mu,sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      i <- 0
      prob <- dPTRTE(x=i,mu = mu, sigma = sigma, log=FALSE)
      F <- prob
      while (p >= F) {
        i <- i + 1
        prob <- dPTRTE(x=i,mu = mu, sigma = sigma, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  # Vectorizar la función auxiliar
  one_quantile_dPTRTE <- Vectorize(one_quantile_dPTRTE)
  # End auxiliary function

  # The quantile
  q <- one_quantile_dPTRTE(p=pp, mu=mu, sigma=sigma)

  # To deal with invalid p's
  q[p <  0] <- NaN
  q[p >  1] <- NaN
  q[p == 1] <- Inf
  q[p == 0] <- 0

  return(q)
}
#' @importFrom stats runif
#' @export
#' @rdname dPTRTE
rPTRTE <- function(n, mu, sigma) {
  # Validaciones
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(sigma <= 0) | any(sigma >= 1))  stop("parameter sigma must be in (0, 1)")
  if (any(n <= 0))    stop(paste("n must be a positive integer", "\n", ""))

  # Generar números uniformes
  n <- ceiling(n)
  u <- runif(n=n)
  x <- qPTRTE(p=u, mu=mu, sigma=sigma)
  return(x)
}
