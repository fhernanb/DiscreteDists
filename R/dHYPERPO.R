#' The hyper-Poisson distribution
#'
#' @description
#' Those functions define the density, distribution function, quantile
#' function and random generation for the Poisson, PO(), distribution
#' with parameters \code{mu} and \code{sigma}.
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param mu vector of positive values of this parameter.
#' @param sigma vector of positive values of this parameter.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#'
#' @return
#' \code{dHYPERPO} gives the density, \code{pHYPERPO} gives the distribution
#' function, \code{qHYPERPO} gives the quantile function, \code{rHYPERPO}
#' generates random deviates.
#'
#' @example  examples/examples_dHYPERPO.R
#'
#' @export
#'
dHYPERPO <- function(x, mu=1, sigma=1, log=FALSE){
  if (any(sigma <= 0))  stop("parameter gamma has to be positive!")
  if (any(mu <= 0)) stop("parameter lambda has to be positive!")
  if (any(x < 0)) stop(paste("x must be >=0", "\n", ""))
  # Begin auxiliar function
  F11 <- function(x, a, c, z) {
    p1 <- gamma(a+x) / gamma(a)
    p2 <- gamma(c+x) / gamma(c)
    p1 * z^x / (p2 * factorial(x))
  }
  # End auxiliar function
  p1 <- lgamma(sigma) + x * log(mu) - lgamma(sigma + x)
  f11 <- add(f=F11, lower=0, upper=99, a=1, c=sigma, z=mu)$value
  p2 <- log(f11)
  res <- p1 - p2
  if(log)
    return(res)
  else
    return(exp(res))
}
