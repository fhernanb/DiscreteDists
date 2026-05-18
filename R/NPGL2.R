#' Mean-parametrized New Poisson-generalised Lindley distribution (NPGL2)
#'
#' @author Tomas Mesa, \email{tomas.mesaz@udea.edu.co}
#'
#' @description
#' The function \code{NPGL2()} defines the mean-parametrized version of the
#' Poisson-generalised Lindley distribution, a two-parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS fitting using
#' the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for
#'   the mu parameter (the mean of the distribution).
#' @param sigma.link defines the sigma.link, with "log" link as the default
#'   for the sigma parameter (parameter alpha of the original NPGL).
#'
#' @references
#' Altun, E.
#' A new two-parameter discrete poisson-generalized Lindley distribution with
#' properties and applications to healthcare data sets. Comput Stat 36,
#' 2841-2861 (2021). https://doi.org/10.1007/s00180-021-01097-0
#'
#' @seealso \link{dNPGL2}, \link{NPGL}.
#'
#' @details
#' This family uses the mean-parametrized version of the NPGL distribution
#' proposed in Section 6 of Altun (2021). The new parameters are:
#' \itemize{
#'   \item \eqn{\mu > 0}: the mean of the distribution, i.e. \eqn{E[X] = \mu}.
#'   \item \eqn{\sigma > 0}: the shape parameter (equivalent to \eqn{\alpha}
#'         in the original NPGL parametrization).
#' }
#'
#' The reparametrization links \eqn{\mu} (mean) and \eqn{\sigma} (\eqn{\alpha})
#' to the original parameter \eqn{\theta} via the positive root of the mean
#' equation \eqn{E[X] = (\alpha + \theta) / (\theta (1 + \theta)) = \mu}:
#'
#' \deqn{\theta = \frac{\sqrt{4\sigma\mu + \mu^2 - 2\mu + 1} - \mu + 1}{2\mu}}
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used to fit a
#' mean-parametrized NPGL distribution in the \code{gamlss()} function.
#'
#' @example examples/examples_NPGL2.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
NPGL2 <- function(mu.link = "log", sigma.link = "log") {

  mstats <- checklink("mu.link", "Mean-parametrized Poisson-generalised Lindley",
                      substitute(mu.link), c("log", "inverse", "own"))
  dstats <- checklink("sigma.link", "Mean-parametrized Poisson-generalised Lindley",
                      substitute(sigma.link), c("log", "inverse", "own"))

  structure(
    list(
      family     = c("NPGL2", "Mean-parametrized Poisson-generalised Lindley"),
      parameters = list(mu = TRUE, sigma = TRUE),
      nopar      = 2,
      type       = "Discrete",

      mu.link    = as.character(substitute(mu.link)),
      sigma.link = as.character(substitute(sigma.link)),

      mu.linkfun    = mstats$linkfun,
      sigma.linkfun = dstats$linkfun,

      mu.linkinv    = mstats$linkinv,
      sigma.linkinv = dstats$linkinv,

      mu.dr    = mstats$mu.eta,
      sigma.dr = dstats$mu.eta,

      # First derivatives

      dldm = function(y, mu, sigma) {
        dm <- gamlss::numeric.deriv(
          expr  = dNPGL2(x = y, mu = mu, sigma = sigma, log = TRUE),
          theta = "mu",
          delta = 1e-04
        )
        as.vector(attr(dm, "gradient"))
      },

      dldd = function(y, mu, sigma) {
        dd <- gamlss::numeric.deriv(
          expr  = dNPGL2(x = y, mu = mu, sigma = sigma, log = TRUE),
          theta = "sigma",
          delta = 1e-04
        )
        as.vector(attr(dd, "gradient"))
      },

      # Second derivatives (outer-product approximation, standard in gamlss)

      d2ldm2 = function(y, mu, sigma) {
        dm <- gamlss::numeric.deriv(
          expr  = dNPGL2(x = y, mu = mu, sigma = sigma, log = TRUE),
          theta = "mu",
          delta = 1e-04
        )
        dldm <- as.vector(attr(dm, "gradient"))
        -dldm * dldm
      },

      d2ldmdd = function(y, mu, sigma) {
        dm <- gamlss::numeric.deriv(
          expr  = dNPGL2(x = y, mu = mu, sigma = sigma, log = TRUE),
          theta = "mu",
          delta = 1e-04
        )
        dldm <- as.vector(attr(dm, "gradient"))

        dd <- gamlss::numeric.deriv(
          expr  = dNPGL2(x = y, mu = mu, sigma = sigma, log = TRUE),
          theta = "sigma",
          delta = 1e-04
        )
        dldd <- as.vector(attr(dd, "gradient"))

        -dldm * dldd
      },

      d2ldd2 = function(y, mu, sigma) {
        dd <- gamlss::numeric.deriv(
          expr  = dNPGL2(x = y, mu = mu, sigma = sigma, log = TRUE),
          theta = "sigma",
          delta = 1e-04
        )
        dldd <- as.vector(attr(dd, "gradient"))
        -dldd * dldd
      },

      G.dev.incr = function(y, mu, sigma, pw = 1, ...) {
        -2 * dNPGL2(x = y, mu = mu, sigma = sigma, log = TRUE)
      },

      rqres = expression(
        rqres(pfun = "pNPGL2", type = "Discrete",
              ymin = 0, y = y, mu = mu, sigma = sigma)
      ),

      mu.initial = expression(
        mu <- rep(mean(y), length(y))
      ),

      sigma.initial = expression(
        sigma <- rep(estim_mu_sigma_NPGL2(y)[2], length(y))
      ),

      mu.valid    = function(mu)    all(mu > 0),
      sigma.valid = function(sigma) all(sigma > 0),

      y.valid = function(y) all(y >= 0),

      mean = function(mu, sigma) {
        mu
      },

      variance = function(mu, sigma) {
        disc  <- 4 * sigma * mu + mu^2 - 2 * mu + 1
        theta <- (sqrt(disc) - mu + 1) / (2 * mu)
        (theta * (sigma^2 + theta * (sigma + theta + 2) + 2) + sigma) /
          (theta^2 * (1 + theta)^2)
      }
    ),
    class = c("gamlss.family", "family")
  )
}
#'
#' estim_mu_sigma_NPGL2
#'
#' @param y vector with the random sample
#' @examples
#' y <- rNPGL2(n = 100, mu = 3, sigma = 2)
#' estim_mu_sigma_NPGL2(y = y)
#' @importFrom stats optim
#' @export
estim_mu_sigma_NPGL2 <- function(y) {
  mod <- optim(
    par     = c(log(mean(y)), 0),
    fn      = logLik_NPGL2,
    method  = "Nelder-Mead",
    control = list(fnscale = -1, maxit = 100000),
    x       = y
  )

  theta_hat <- exp(mod$par[1])
  alpha_hat <- exp(mod$par[2])

  mu_hat <- (alpha_hat + theta_hat) / (theta_hat * (1 + theta_hat))

  res <- c(mu_hat = mu_hat, sigma_hat = alpha_hat)
  return(res)
}
#'
#' logLik_NPGL2
#'
#' Auxiliary function to evaluate the log-likelihood of the mean-parametrized
#' NPGL2 distribution.
#'
#' @param param vector with the values for mu and sigma
#' @param x vector with the data
#' @examples
#' y <- rNPGL2(n = 100, mu = 3, sigma = 2)
#' logLik_NPGL2(param = c(log(3), log(2)), x = y)
#' @export
logLik_NPGL2 <- function(param = c(log(3), log(2)), x) {
  return(sum(
    dNPGL2(x     = x,
           mu    = exp(param[1]),
           sigma = exp(param[2]),
           log   = TRUE)
  ))
}

################################################################################
#'
