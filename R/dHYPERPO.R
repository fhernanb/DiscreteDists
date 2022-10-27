#' The hyper-Poisson distribution
#'
#' @description
#' Those functions define the density, distribution function, quantile
#' function and random generation for the hyper-Poisson, HYPERPO(), distribution
#' with parameters \code{mu} and \code{sigma}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of positive values of this parameter.
#' @param sigma vector of positive values of this parameter.
#' @param n number of random values to return
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
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
#' @export
#' @rdname dHYPERPO
pHYPERPO <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(sigma <= 0))  stop("parameter gamma has to be positive!")
  if (any(mu <= 0)) stop("parameter lambda has to be positive!")
  if (any(q < 0)) stop(paste("x must be >=0", "\n", ""))
  ly <- max(length(q), length(mu), length(sigma))
  q <- rep(q, length = ly)
  mu <- rep(mu, length = ly)
  sigma <- rep(sigma, length = ly)
  # Begin auxiliar function
  aux_func <- function(q, mu, sigma) {
    cdf <- numeric(length(q))
    for (i in 1:length(q)) {
      res <- dHYPERPO(x=0:q[i], mu=mu[i], sigma=sigma[i], log=FALSE)
      cdf[i] <- sum(res)
    }
    cdf
  }
  # End auxiliar function
  cdf <- aux_func(q=q, mu=mu, sigma=sigma)
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
#' @rdname dHYPERPO
rHYPERPO <- function(n, mu=1, sigma=1) {
  if (any(sigma <= 0))  stop("parameter gamma has to be positive!")
  if (any(mu <= 0)) stop("parameter lambda has to be positive!")
  if (any(n <= 0)) stop(paste("n must be a positive integer", "\n", ""))
  # Begin auxiliar function
  one_random_hyperpo <- function(u, mu, sigma) {
    p <- dHYPERPO(x=0, mu=mu, sigma=sigma, log=FALSE)
    F <- p
    i <- 0
    while (u >= F) {
      i <- i + 1
      p <- dHYPERPO(x=i, mu=mu, sigma=sigma, log=FALSE)
      F <- F + p
    }
    return(i)
  }
  one_random_hyperpo <- Vectorize(one_random_hyperpo)
  # End auxiliar function
  one_random_hyperpo(u=runif(n), mu, sigma)
}
#' @export
#' @rdname dHYPERPO
qHYPERPO <- function(p, mu = 1, sigma = 1, lower.tail = TRUE,
                     log.p = FALSE) {
  if (any(mu <= 0))
    stop(paste("mu must be greater than 0 ", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be greater than 0 ", "\n", ""))
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))
  if (log.p == TRUE)
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE)
    p <- p
  else p <- 1 - p
  # Begin auxiliar function
  one_quantile_hyperpo <- function(p, mu, sigma) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      prob <- dHYPERPO(x=0, mu=mu, sigma=sigma, log=FALSE)
      F <- prob
      i <- 0
      while (p >= F) {
        i <- i + 1
        prob <- dHYPERPO(x=i, mu=mu, sigma, log=FALSE)
        F <- F + prob
      }
    }
    return(i)
  }
  one_quantile_hyperpo <- Vectorize(one_quantile_hyperpo)
  # End auxiliar function
  one_quantile_hyperpo(p=p, mu=mu, sigma=sigma)
}
