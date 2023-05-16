#' The Discrete Generalized Burr Hatke family
#'
#' @author Sofia Cuartas García, \email{scuartasg@unal.edu.co}
#'
#' @description
#' The function \code{DBH()} defines the Discrete Burr Hatke distribution, one-parameter
#' discrete distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "logit" link as the default for the mu parameter. Other links are "probit" and "cloglog"'(complementary log-log)
#'
#' @references
#' \insertRef{yousof2021}{DiscreteDists}
#'
#' @importFrom Rdpack reprompt
#'
#' @seealso \link{dDGBH}.
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
#'
#' Note: in this implementation we changed the original parameters \eqn{\lambda} for \eqn{\mu},
#' we did it to implement this distribution within gamlss framework.
#'
#' @example  examples/examples_DGBH.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
DGBH <- function (mu.link="logit") {
  mstats <- checklink("mu.link", "DBH",
                      substitute(mu.link), c("logit", "probit", "cloglog", "cauchit"))

  structure(list(family=c("DBH", "Burr Hatke"),
                 parameters=list(mu=TRUE),
                 nopar=1,
                 type="Discrete",
                 mu.link    = as.character(substitute(mu.link)),
                 mu.linkfun = mstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 mu.dr      = mstats$mu.eta,

                 # Primeras derivadas, por ahora son computacionales

                 dldm = function(y, mu) {
                   dm   <- gamlss::numeric.deriv(dDBH(y, mu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },

                 # Segundas derivadas, por ahora son computacionales

                 d2ldm2 = function(y, mu) {
                   dm   <- gamlss::numeric.deriv(dDBH(y, mu, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 G.dev.incr = function(y, mu, pw = 1, ...) -2*dDBH(y, mu, log=TRUE),
                 rqres      = expression(rqres(pfun="pDBH", type="Discrete",
                                               ymin = 0, y = y, mu = mu)),


                 mu.initial = expression(mu <- rep(estim_mu_DBH(y), length(y)) ),

                 mu.valid = function(mu) all(0 < mu & mu < 1),

                 y.valid = function(y) all(y >= 0),

                 mean = function(mu) {
                   -log(1-mu)/mu - 1
                 },

                 variance = function(mu) {
                   media <- -log(1-mu)/mu - 1
                   varianza <- mu/2*((2*(mu-3))/(mu*(mu-1))+(6*log(1-mu))/mu^2) - media^2
                   return(varianza)
                 }

  ),
  class=c("gamlss.family", "family"))
}
