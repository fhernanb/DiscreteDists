#' The GGEO family
#'
#' @author Valentina Hurtado Sepúlveda, \email{vhurtados@unal.edu.co}
#'
#' @description
#' The function \code{GGEO()} defines the Generalized Geometric distribution,
#' a two parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "logit" link as the default for the sigma. Other links are "probit" and "cloglog"'(complementary log-log)
#'
#' @references
#' \insertRef{gomez2010}{DiscreteDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @seealso \link{dGGEO}.
#'
#' @details
#' The GGEO distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\mu \sigma^x (1-sigma)}{(1-(1-mu) \sigma^{x+1})(1-(1-mu) \sigma^{x})}}
#'
#' with \eqn{\mu > 0} and \eqn{0 < \sigma < 1}
#'
#' @example examples/examples_GGEO.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
GGEO <- function (mu.link="log", sigma.link="logit") {

  mstats <- checklink("mu.link", "GGEO",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "GGEO",
                      substitute(sigma.link), c("logit", "probit", "cloglog", "cauchit"))

  structure(list(family=c("GGEO", "Generalized geometric"),
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

                 # Primeras derivadas manuales

                 dldm = function(y, mu, sigma) {
                   dldm<- (-sigma^(y+1)/(1-(1-mu)*sigma^(y+1)))-(sigma^y/(1-(1-mu)*sigma^y))+1/mu
                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   dldd <- (1-mu)*y*sigma^(y-1)/(1-(1-mu)*sigma^y)+(1-mu)*(y+1)*sigma^y/((1-(1-mu)*sigma^(y+1)))+y/sigma-1/(1-sigma)
                   dldd
                 },

                 # Segundas derivadas manuales

                 d2ldm2 = function(y, mu, sigma) {
                   dldm <- (-sigma^(y+1)/(1-(1-mu)*sigma^(y+1)))-(sigma^y/(1-(1-mu)*sigma^y))+1/mu
                   d2ldm2 <- - dldm * dldm
                   d2ldm2
                 },
                 d2ldmdd = function(y, mu, sigma) {
                   dldm <- (-sigma^(y+1)/(1-(1-mu)*sigma^(y+1)))-(sigma^y/(1-(1-mu)*sigma^y))+1/mu

                   dldd <- (1-mu)*y*sigma^(y-1)/(1-(1-mu)*sigma^y)+(1-mu)*(y+1)*sigma^y/((1-(1-mu)*sigma^(y+1)))+y/sigma-1/(1-sigma)

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd
                 },
                 d2ldd2  = function(y, mu, sigma) {
                   dldd <- (1-mu)*y*sigma^(y-1)/(1-(1-mu)*sigma^y)+(1-mu)*(y+1)*sigma^y/((1-(1-mu)*sigma^(y+1)))+y/sigma-1/(1-sigma)
                   d2ldd2 <- - dldd * dldd
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dGGEO(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pGGEO", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),

                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_GGEO(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_GGEO(y)[2], length(y)) ),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(0 < sigma & sigma < 1),

                 y.valid = function(y) all(y >= 0)

  ),
  class=c("gamlss.family", "family"))
}

