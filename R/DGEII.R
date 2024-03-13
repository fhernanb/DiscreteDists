#' The DGEII family
#'
#' @author Valentina Hurtado Sepúlveda, \email{vhurtados@unal.edu.co}
#'
#' @description
#' The function \code{DGEII()} defines the Discrete generalized exponential distribution,
#' Second type, a two parameter
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "logit" link as the default for the sigma. Other links are "probit" and "cloglog"'(complementary log-log)
#'
#' @references
#' \insertRef{nekoukhou2013}{DiscreteDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @seealso \link{dDGEII}.
#'
#' @details
#' The DGEII distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = (1-\sigma^{x+1})^{\mu}-(1-\sigma^x)^{\mu}}
#'
#' with \eqn{\mu > 0} and \eqn{0 < \sigma < 1}
#'
#' @example examples/examples_DGEII.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
DGEII <- function (mu.link="log", sigma.link="logit") {

  mstats <- checklink("mu.link", "DGEII",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "DGEII",
                      substitute(sigma.link), c("logit", "probit", "cloglog", "cauchit"))

  structure(list(family=c("DGEII", "Discrete generalized exponential distribution of a second type II"),
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

                 # Primeras derivadas por ahora computacionales

                 dldm = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },


                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },

                 # Segundas derivadas por ahora computacionales

                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dDGEII(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dDGEII(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pDGEII", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),

                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_DGEII(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_DGEII(y)[2], length(y)) ),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(0 < sigma & sigma < 1),

                 y.valid = function(y) all(y >= 0)

  ),
  class=c("gamlss.family", "family"))
}

