#' The hyper-Poisson distribution (with mu as mean)
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the hyper-Poisson in
#' the second parameterization with parameters \eqn{\mu} (as mean) and
#' \eqn{\sigma} as the dispersion parameter.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of positive values of this parameter.
#' @param sigma vector of positive values of this parameter.
#' @param n number of random values to return
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \code{P[X > x]}.
#'
#' @references
#' \insertRef{saez2013hyperpo}{DiscreteDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @seealso \link{HYPERPO2}, \link{HYPERPO}.
#'
#' @details
#' The hyper-Poisson distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ...
#'
#' Note: in this implementation the parameter \eqn{\mu} is the mean
#' of the distribution and \eqn{\sigma} corresponds to
#' the dispersion parameter. If you fit a model with this parameterization,
#' the time will increase because an internal procedure to convert \eqn{\mu}
#' to \eqn{\lambda} parameter.
#'
#' @return
#' \code{dHYPERPO2} gives the density, \code{pHYPERPO2} gives the distribution
#' function, \code{qHYPERPO2} gives the quantile function, \code{rHYPERPO2}
#' generates random deviates.
#'
#' @example  examples/examples_dHYPERPO2.R
#'
#' @export
#'
dHYPERPO2 <- function(x, mu=1, sigma=1, log=FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(x < 0))       stop(paste("x must be >=0", "\n", ""))
  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda(media=mu, gamma=sigma)
  p1 <- x * log(mu) - log(AR(a=sigma, r=x)) # AR is an util function
  f11 <- F11(c=sigma, z=mu) # F11 is an util function
  p2 <- log(f11)
  res <- p1 - p2
  if(log)
    return(res)
  else
    return(exp(res))
}
#' @export
#' @rdname dHYPERPO2
pHYPERPO2 <- function(q, mu=1, sigma=1, lower.tail = TRUE, log.p = FALSE){
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(q < 0))       stop(paste("x must be >=0", "\n", ""))
  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda(media=mu, gamma=sigma)
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
#' @rdname dHYPERPO2
rHYPERPO2 <- function(n, mu=1, sigma=1) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(n <= 0))      stop(paste("n must be a positive integer", "\n", ""))
  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda(media=mu, gamma=sigma)
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
#' @rdname dHYPERPO2
qHYPERPO2 <- function(p, mu = 1, sigma = 1, lower.tail = TRUE,
                     log.p = FALSE) {
  if (any(sigma <= 0))  stop("parameter sigma has to be positive!")
  if (any(mu <= 0))     stop("parameter mu has to be positive!")
  if (any(p < 0) | any(p > 1.0001))
    stop(paste("p must be between 0 and 1", "\n", ""))
  # To obtain the mu in the older parameterization
  mu <- obtaining_lambda(media=mu, gamma=sigma)
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
