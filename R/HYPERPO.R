#' The hyper Poisson family
#'
#' @description
#' The function \code{HYPERPO()} defines the hyper Poisson distribution, a two parameter
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#'
#' @details
#' The hyper Poisson with parameters \code{mu} and \code{sigma}
#' has mass function given by
#'
#' \eqn{f(x) = bla bla bla}
#'
#' for x = 0, 1, 2, ....
#'
#' @example examples/examples_HYPERPO.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
HYPERPO <- function (mu.link="log", sigma.link="log") {

  mstats <- checklink("mu.link", "Generalized exponential-gaussian",
                      substitute(mu.link), c("identity", "own"))
  dstats <- checklink("sigma.link", "Generalized exponential-gaussian",
                      substitute(sigma.link), c("log", "own"))

  structure(list(family=c("HYPERPO", "hyper Poisson"),
                 parameters=list(mu=TRUE, sigma=TRUE),
                 nopar=2,
                 type="Discrete",

                 mu.link    = as.character(substitute(mu.link)),
                 sigma.link = as.character(substitute(sigma.link)),

                 mu.linkfun    = mstats$linkfun,
                 sigma.linkfun = dstats$linkfun,

                 mu.linkinv    = mstats$linkinv,
                 sigma.linkinv = dstats$linkinv,

                 mu.dr    = mstats$mu.eta,
                 sigma.dr = dstats$mu.eta,

                 # Primeras derivadas, por ahora son computacionales

                 dldm = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },

                 # Segundas derivadas, por ahora son computacionales

                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE), "mu", delta=1e-04)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dHYPERPO(y, mu, sigma, log=TRUE), "sigma", delta=1e-04)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, ...) -2*dHYPERPO(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pHYPERPO", type="Discrete", ymin = 0, y = y, mu = mu, sigma = sigma)),

                 mu.initial    = expression({ mu <- (y + mean(y))/2 }),
                 sigma.initial = expression(sigma <- rep(1, length(y))),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),

                 y.valid = function(y) TRUE
  ),
  class=c("gamlss.family", "family"))
}

