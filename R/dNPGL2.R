#' New Poisson-generalised Lindley distribution (with mu as mean)
#'
#' @author Tomás Mesa, \email{tomas.mesaz@udea.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Poisson-generalised Lindley in second
#' parametrization with parameters \eqn{\mu} (as mean) and \eqn{\sigma}.
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
#' properties and applications to healthcare data sets. Comput Stat 36, 2841–2861 (2021). https://doi.org/10.1007/s00180-021-01097-0
#'
#' @seealso \link{NPGL}, \link{NPGL2}.
#'
#' @details
#' The Poisson-generalised Lindley distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has support \eqn{x = 0, 1, 2, \ldots} and probability mass function given by
#'
#' Note: in this implementation the parameter \eqn{\mu} is the mean
#' of the distribution and \eqn{\sigma} (equivalent to \eqn{\alpha}
#' in the original NPGL parametrization).
#'
#' @return
#' \code{dNPGL2} gives the density, \code{pNPGL2} gives the distribution
#' function, \code{qNPGL2} gives the quantile function, \code{rNPGL2}
#' generates random deviates.
#'
#' @example examples/examples_dNPGL2.R
#'
#' @export
#' @useDynLib DiscreteDists
#' @importFrom Rcpp sourceCpp
dNPGL2 <- function(x, mu = 2, sigma = 2, log = FALSE) {
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")

  # Ensure same length vector
  ly  <- max(length(x), length(mu), length(sigma))
  xx  <- rep(x,     length = ly)
  mu  <- rep(mu,    length = ly)
  sig <- rep(sigma, length = ly)

  discriminante <- 4 * sig * mu + mu^2 - 2 * mu + 1
  theta <- (sqrt(discriminante) - mu + 1) / (2 * mu)

  p <- dNPGL(x = xx, mu = theta, sigma = sig, log = TRUE)

  if (log == TRUE)
    return(p)
  else
    return(exp(p))
}
#' @export
#' @rdname dNPGL2
pNPGL2 <- function(q, mu = 2, sigma = 2, lower.tail = TRUE, log.p = FALSE) {
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")

  # Ensure same length vector
  ly <- max(length(q), length(mu), length(sigma))
  qq  <- rep(q,     length = ly)
  mu  <- rep(mu,    length = ly)
  sig <- rep(sigma, length = ly)

  # Temporal change for invalid q's
  qq[q < 0]   <- 0
  qq[q == Inf] <- 0
  qq <- as.integer(qq)

  # Calcular theta desde mu (media) y sigma (alfa)
  discriminante <- 4 * sig * mu + mu^2 - 2 * mu + 1
  theta <- (sqrt(discriminante) - mu + 1) / (2 * mu)

  # Auxiliary function: delega en pNPGL original
  fn <- function(q, theta, sig) sum(dNPGL(x = 0:q, mu = theta, sigma = sig))
  Vec_fn <- Vectorize(fn)

  cdf <- Vec_fn(q = qq, theta = theta, sig = sig)

  cdf[q < 0]   <- 0
  cdf[q == Inf] <- 1

  if (lower.tail == FALSE) cdf <- 1 - cdf
  if (log.p == TRUE)       cdf <- log(cdf)

  return(cdf)
}
#' @importFrom stats runif
#' @export
#' @rdname dNPGL2
rNPGL2 <- function(n, mu = 2, sigma = 2) {
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  if (any(mu <= 0))    stop("parameter mu has to be positive!")
  if (any(n <= 0))     stop(paste("n must be a positive integer", "\n", ""))

  n <- ceiling(n)

  # Calcular theta desde mu (media) y sigma (alfa)
  discriminante <- 4 * sigma * mu + mu^2 - 2 * mu + 1
  theta <- (sqrt(discriminante) - mu + 1) / (2 * mu)

  # Delegar en rNPGL original con la parametrización (theta, alfa)
  x <- rNPGL(n = n, mu = theta, sigma = sigma)

  return(x)
}
#' @export
#' @rdname dNPGL2
qNPGL2 <- function(p, mu = 2, sigma = 2, lower.tail = TRUE, log.p = FALSE) {
  if (any(sigma <= 0)) stop("parameter sigma has to be positive!")
  if (any(mu <= 0))    stop("parameter mu has to be positive!")

  if (log.p == TRUE)       p <- exp(p)
  if (lower.tail == FALSE) p <- 1 - p

  # Ensure same length vector
  ly <- max(length(p), length(mu), length(sigma))
  pp  <- rep(p,     length = ly)
  mu  <- rep(mu,    length = ly)
  sig <- rep(sigma, length = ly)

  # Temporal change for invalid p's
  pp[p < 0]  <- 0.1
  pp[p > 1]  <- 0.1
  pp[p == 1] <- 0.1
  pp[p == 0] <- 0.1

  # Calcular theta desde mu (media) y sigma (alfa)
  discriminante <- 4 * sig * mu + mu^2 - 2 * mu + 1
  theta <- (sqrt(discriminante) - mu + 1) / (2 * mu)

  # Begin auxiliary function
  one_quantile_npgl2 <- function(p, theta, sig) {
    if (p + 1e-09 >= 1)
      i <- Inf
    else {
      i <- 0
      prob <- dNPGL(x = i, mu = theta, sigma = sig, log = FALSE)
      F    <- prob
      while (p >= F) {
        i    <- i + 1
        prob <- dNPGL(x = i, mu = theta, sigma = sig, log = FALSE)
        F    <- F + prob
      }
    }
    return(i)
  }
  one_quantile_npgl2 <- Vectorize(one_quantile_npgl2)

  q <- one_quantile_npgl2(p = pp, theta = theta, sig = sig)

  q[p < 0]  <- NaN
  q[p > 1]  <- NaN
  q[p == 1] <- Inf
  q[p == 0] <- 0

  return(q)
}
