#' The Discrete Perks distribution
#'
#' @author Veronica Seguro Varela, \email{vseguro@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Discrete Perks, DPERKS(),
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
#' Tyagi, A., Choudhary, N., & Singh, B. (2020). A new discrete
#' distribution: Theory and applications to discrete failure
#' lifetime and count data. J. Appl. Probab. Statist, 15, 117-143.
#'
#' @seealso \link{DPERKS}.
#'
#' @details
#' The discrete Perks distribution with parameters \eqn{\mu > 0} and \eqn{\sigma > 0}
#' has a support 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\mu(1+\mu)(e^\sigma-1)e^{\sigma x}}{(1+\mu e^{\sigma x})(1+\mu e^{\sigma(x+1)})} }
#'
#' Note: in this implementation we changed the original parameters
#' \eqn{\lambda} for \eqn{\mu} and \eqn{\beta} for \eqn{\sigma},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dDPERKS} gives the density, \code{pDPERKS} gives the distribution
#' function, \code{qDPERKS} gives the quantile function, \code{rDPERKS}
#' generates random deviates.
#'
#' @example examples/examples_dDPERKS.R
#'
#' @export
#'
dDPERKS <- function(x, mu=0.5, sigma=0.5, log=FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(x < 0))       stop(paste("x must be >=0", "\n", ""))
  p <- log(mu) + log(1+mu) + log(exp(sigma) - 1) + sigma*x - log(1+mu*exp(sigma*x)) - log(1 + mu*exp(sigma*(x+1)))
  if(log){
    return(p)}
  else{
    return(exp(p))}
}
#' @export
#' @rdname dDPERKS
pDPERKS <- function(q, mu=0.5, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(q < 0))       stop(paste("x must be >=0", "\n", ""))
  cdf <- mu*(exp(sigma*(q+1))- 1)/(1+mu*exp(sigma*(q+1)))
  if (lower.tail == TRUE)
    cdf <- cdf
  else cdf = 1 - cdf
  if (log.p == FALSE)
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @importFrom stats runif
#' @export
#' @rdname dDPERKS
rDPERKS <- function(n, mu=0.5, sigma=0.5) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))
  ## Begin auxiliar function
  one_perks <- function(u, mu, sigma){
    if (any(mu <= 0)) stop("parameter mu has to be positive!")
    if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
    p = mu * (exp(sigma) - 1)/(1 + mu * exp(sigma)) #cumulative distribution function
    acum = p
    i = 0
    while(u >= acum){
      #i = i + 1
      p = exp(sigma) * (1 + mu*exp(sigma*i))/(1+ mu*exp(sigma*(i+2))) * p #recurrence relation for probabilities
      acum = acum + p
      i = i + 1
     }
     return(i)
    }
  one_perks <- Vectorize(one_perks)
  ## End auxiliar function
  one_perks(runif(n), mu, sigma)
}
#' @export
#' @rdname dDPERKS
qDPERKS <- function(p, mu=0.5, sigma=0.5, lower.tail=TRUE, log.p=FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1.0001)) stop(paste("p must be between 0 and 1", "\n", ""))
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  # Begin aux fun
  one_quantile <- function(p, mu, sigma){
    quan <- ceiling(((1/sigma)*log((p + mu)/(mu*(1-p)))) - 1)
    return(quan)
  }
  one_quantile <- Vectorize(one_quantile)
  # End aux fun
  one_quantile(p, mu, sigma)
}

