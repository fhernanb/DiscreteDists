#' Discrete power-Ailamujia distribution
#'
#' @author Maria Camila Mena Romana, \email{mamenar@unal.edu.co}
#'
#' @description
#' The function \code{DsPA()} defines the discrete power-Ailamujia distribution, a two parameter
#' distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "logit" link as the default for the sigma parameter. Other links are "probit" and "cloglog"'(complementary log-log).
#'
#' @references
#' Alghamdi, A. S., Ahsan-ul-Haq, M., Babar, A., Aljohani, H. M., Afify, A. Z., & Cell, Q. E. (2022). The discrete power-Ailamujia distribution: properties, inference, and applications. AIMS Math, 7(5), 8344-8360.
#'
#' @seealso \link{dDsPA}.
#'
#' @details
#' The discrete power-Ailamujia distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and density given by
#'
#' \eqn{f(x | \mu, \sigma) = (\sigma^x)^\mu (1-x^\mu \log(\lambda)) -  (\sigma^{(x+1)})^\mu (1-(x+1)^\mu \log(\lambda))}
#'
#' Note: in this implementation we changed the original parameters \eqn{\beta} and \eqn{\lambda}
#' for \eqn{\mu} and \eqn{\sigma} respectively, we did it to implement this distribution within gamlss framework.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a discrete power-Ailamujia distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_DsPA.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
DsPA <- function (mu.link="log", sigma.link="logit") {

  mstats <- checklink("mu.link", "DsPA",
                      substitute(mu.link), c("log", "inverse", "own"))
  dstats <- checklink("sigma.link", "DsPA",
                      substitute(sigma.link), c("logit", "probit", "cloglog", "own"))

  structure(list(family=c("DsPA", "discrete power-Ailamujia"),
                 parameters=list(mu=TRUE, sigma=TRUE),
                 nopar=2,
                 type="Discrete",

                 mu.link   =as.character(substitute(mu.link)),
                 sigma.link=as.character(substitute(sigma.link)),

                 mu.linkfun   =mstats$linkfun,
                 sigma.linkfun=dstats$linkfun,

                 mu.linkinv   =mstats$linkinv,
                 sigma.linkinv=dstats$linkinv,

                 mu.dr   =mstats$mu.eta,
                 sigma.dr=dstats$mu.eta,

                 # First derivates

                 dldm = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDsPA(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dDsPA(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },

                 # Second derivates

                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDsPA(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dDsPA(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.00001)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dDsPA(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dDsPA(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.00001)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },

                 G.dev.incr=function(y, mu, sigma, pw=1, ...) -2*dDsPA(y, mu, sigma, log=TRUE),
                 rqres     =expression(rqres(pfun="pDsPA", type="Discrete",
                                             ymin=0, y=y, mu=mu, sigma=sigma)),

                 mu.initial   =expression(mu    <- rep(estim_mu_sigma_DsPA(y)[1], length(y)) ),
                 sigma.initial=expression(sigma <- rep(estim_mu_sigma_DsPA(y)[2], length(y)) ),

                 mu.valid   =function(mu)       all(mu > 0),
                 sigma.valid=function(sigma)    all(sigma > 0) && all(sigma < 1),

                 y.valid=function(y)  TRUE

  ),
  class=c("gamlss.family", "family"))
}
#'
#' estim_mu_sigma_DsPA
#'
#' This function generates initial values for DsPA distribution
#'
#' @param y vector with the random sample
#' @examples
#' y <- rDsPA(n=100, mu=0.5 , sigma=0.5)
#' estim_mu_sigma_DsPA(y=y)
#' @importFrom stats optim
#' @export
estim_mu_sigma_DsPA <- function(y){
  mod <- optim(par=c (0, 0),
               fn=logLik_DsPA,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)

  logit_inv <- function(x) exp(x) / (exp(x) + 1)

  res <- c(mu_hat       = exp(mod$par[1]),
           sigma_hat    = logit_inv(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#'
#' logLik_DsPA
#' This is an auxiliar function to obtain the logLik for DsPA
#'
#' @param param  vector with the values for mu and sigma
#' @param  x vector with the data
#' @examples
#' y <- rDsPA(n=100, mu=0.5 , sigma=0.5)
#' logLik_DsPA(param=c(0, 0), x=y)
#' @importFrom stats optim
#' @export
logLik_DsPA <- function(param=c (0, 0), x){
  logit_inv <- function(x) exp(x) / (exp(x) + 1)
  return(sum(dDsPA(x,
                   mu    = exp(param[1]),
                   sigma = logit_inv(param[2]),
                   log=TRUE)))
}
