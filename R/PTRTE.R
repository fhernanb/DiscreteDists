#' Poisson-transmuted record type exponential distribution
#'
#' @author Rebeca Isabel Rodriguez Gonzalez, \email{rebeca.rodriguez@udea.edu.co}
#'
#' @description
#' The function \code{PTRTE()} defines the Poisson-transmuted record type exponential distribution,
#' a two-parameter discrete distribution, for a \code{gamlss.family} object to be used in GAMLSS fitting
#' using the function \code{gamlss()}.
#'
#' @param mu.link defines the link function for the mu parameter, with "log" as the default.
#' @param sigma.link defines the link function for the sigma parameter, with "logit" as the default.
#'
#' @references
#' Erbayram, T., & Akdogan, Y. (2025).
#' A new discrete model generated from mixed Poisson transmuted
#' record type exponential distribution.
#' Ricerca di Matematica, 74, 1225–1247.
#'
#' @seealso \link{dPTRTE}.
#'
#' @details
#' The Poisson-transmuted record type exponential distribution with parameters \eqn{\mu}
#' and \eqn{\sigma} has support \eqn{x = 0,1,2,\dots} and probability mass function given by
#'
#' \deqn{f(x | \mu, \sigma) = \frac{\mu(\sigma x\mu + 1 + \mu - \sigma)}{(1+\mu)^{x+2}}}
#'
#' Parameter restrictions:
#' \eqn{\mu > 0}  and \eqn{0 < \sigma < 1}.
#'
#' Note: we renamed the original parameters \eqn{\theta} and \eqn{p} to \eqn{\mu} and \eqn{\sigma} respectively
#' to implement this distribution within the \code{gamlss} framework.
#'
#' @return
#' Returns a \code{gamlss.family} object that can be used to fit the Poisson-transmuted record type exponential
#' distribution using the \code{gamlss()} function.
#'
#' @example examples/examples_PTRTE.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @importFrom stats dnorm pnorm
#' @export
PTRTE <- function (mu.link="log", sigma.link="logit") {

  mstats <- checklink("mu.link", "PTRTE",
                      substitute(mu.link), c("log","inverse", "own"))
  dstats <- checklink("sigma.link", "PTRTE",
                      substitute(sigma.link), c("logit", "probit", "cloglog", "own"))

  structure(list(family=c("PTRTE", "Poisson-Transmuted Record Type Exponential"),
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


                 # Primeras derivadas

                 dldm = function(y, mu, sigma) {
                   term1 <- 1/mu
                   term2 <- (sigma * y + 1 ) / (sigma * mu * y + 1 + mu - sigma)
                   term3 <- (y + 2) / (1 + mu)

                   dldm <- term1 + term2 - term3
                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   numerator <- mu * y - 1
                   denominator <- sigma * mu * y + 1 + mu - sigma

                   dldd <- numerator / denominator
                   dldd
                 },

                 # Segundas derivadas

                 d2ldm2 = function(y, mu, sigma) {
                   term1 <- 1/mu
                   term2 <- (sigma * y + 1) / (sigma * mu * y + 1 + mu - sigma)
                   term3 <- (y + 2) / (1 + mu)
                   dldm <- term1 + term2 - term3
                   d2ldm2 <- -dldm * dldm
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {

                   #Derivada de mu
                   term1 <- 1/mu
                   term2 <- (sigma * y + 1) / (sigma * mu * y + 1 + mu - sigma)
                   term3 <- (y + 2) / (1 + mu)
                   dldm <- term1 + term2 - term3

                   #Derivada de sigma
                   numerator <- mu * y - 1
                   denominator <- sigma * mu * y + 1 + mu - sigma
                   dldd <- numerator / denominator

                   d2ldmdd <- -dldm*dldd

                   d2ldmdd

                 },

                 d2ldd2  = function(y, mu, sigma) {
                   numerator <- mu * y - 1
                   denominator <- sigma * mu * y + 1 + mu - sigma
                   dldd <- numerator / denominator

                   d2ldd2 <- - dldd * dldd
                   d2ldd2
                 },


                 G.dev.incr = function(y, mu, sigma, pw=1, ...) -2*dPTRTE(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pPTRTE", type="Discrete",ymin=0, y=y, mu=mu, sigma=sigma)),

                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_PTRTE(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_PTRTE(y)[2], length(y)) ),


                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0) && all(sigma < 1) ,

                 y.valid = function(y) all( y>=0)
  ),
  class=c("gamlss.family", "family"))
}
#' Initial values for PTRTE
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the MLE estimations.
#' @keywords internal
#' @keywords noRd
#' @export
#' @importFrom stats optim
estim_mu_sigma_PTRTE <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_PTRTE,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  logit_inv <- function(x) exp(x) / (exp(x)+1)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = logit_inv(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
logLik_PTRTE <- function(param=c(0, 0), x) {
  logit_inv <- function(x) exp(x) / (exp(x)+1)
  return(sum(dPTRTE(x,
                    mu    = exp(param[1]),
                    sigma = logit_inv(param[2]),
                    log=TRUE)))
}
