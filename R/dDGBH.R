#' The Discrete Generalized Burr Hatke distribution
#'
#' @author Sofia Cuartas García, \email{scuartasg@unal.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Discrete Generalized Burr Hatke distribution
#' with parameter \eqn{\mu}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param sigma vector of the sigma parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \code{P[X > x]}.
#'
#' @references
#' \insertRef{yousof2021}{DiscreteDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @seealso \link{DGBH}.
#'
#' @details
#' The Discrete Burr-Hatke distribution with parameters \eqn{\mu} has a support
#' 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu) = (\frac{1}{x+1}-\frac{\mu}{x+2})\mu^{x}}
#'
#'
#' The pmf is log-convex for all values of \eqn{0 < \mu < 1}, where \eqn{\frac{f(x+1;\mu)}{f(x;\mu)}}
#' is an increasing function in \eqn{x} for all values of the parameter \eqn{\mu}.
#'
#' Note: in this implementation we changed the original parameters \eqn{\lambda} for \eqn{\mu},
#' we did it to implement this distribution within gamlss framework.
#'
#' @return
#' \code{dDGBH} gives the density, \code{pDGBH} gives the distribution
#' function, \code{qDGBH} gives the quantile function, \code{rDGBH}
#' generates random deviates.
#'
#' @example  examples/examples_dDGBH.R
#'
#' @importFrom dplyr between
#' @export
#'
dDGBH <- function(x, mu = 0.999, sigma = 100, log = FALSE){
  if (any(x %% 1 != 0))       stop(paste("x must be integer", "\n", ""))
  if (any(!between(mu,0,1) )) stop(paste("mu must  be between 0 and 1", "\n", ""))
  if (any(sigma <= 0))        stop(paste("sigma must be positive", "\n", ""))
  if(x==0) p <- 1 - mu/2
  else p <- mu^(x^sigma)/(x+1)-mu^((x+1)^sigma)/(x+2)
  if(log == FALSE)
    return(p)
  else
    return(exp(p))
}
dDGBH <- Vectorize(dDGBH)
#' @importFrom dplyr between
#' @export
#' @rdname dDGBH
pDGBH <- function(q, mu= 0.999, sigma= 100, lower.tail = TRUE, log.p = FALSE){
  if (any(q %% 1 != 0))         stop(paste("q must be integer", "\n", ""))
  if (any(!between(mu,0,1) ))   stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0))          stop(paste("sigma must be positive", "\n", ""))

  cdf <- 1 -  (mu^(q+1)^sigma)/(q+2) #FDA
  if (lower.tail == TRUE) {
    cdf <- cdf
  }
  else {
    cdf <- 1 - cdf
  }
  if (log.p == FALSE){
    cdf <- cdf}
  else {cdf <- log(cdf)}
  return(cdf)
}
#' @importFrom dplyr between
#' @export
#' @rdname dDGBH
qDGBH <- function(p, mu = mu, sigma = sigma, lower.tail = TRUE) {
  if (any(!between(mu,0,1) ))       stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0))              stop(paste("sigma must be positive", "\n", ""))
  if (any(p < 0) | any(p > 1.0001)) stop(paste("p must be between 0 and 1", "\n", ""))

  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p

  one_quantile_dgbh <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {

      prob <- dDGBH(0, mu=mu, sigma=sigma)
      F <- prob
      i <- 0
      while (p >= F) {
        i <- i + 1
        prob <- dDGBH(x=i, mu=mu, sigma)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_dgbh <- Vectorize(one_quantile_dgbh)

  one_quantile_dgbh(p=p, mu=mu, sigma=sigma)
}
#' @importFrom dplyr between
#' @export
#' @rdname dDGBH
rDGBH <- function(n,  mu = 0.999 , sigma = 100){
  runif(n=n)
}
rDGBH <- Vectorize(rDGBH)
