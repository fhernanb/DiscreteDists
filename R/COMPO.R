#' The COMPO family
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @description
#' The function \code{COMPO()} defines the
#' Conway-Maxwell-Poisson distribution,
#' a two parameter distribution, for a \code{gamlss.family}
#' object to be used in GAMLSS fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#'
#' @references
#' Shmueli, G., Minka, T. P., Kadane, J. B., Borle, S., & Boatwright, P. (2005).
#' A useful distribution for fitting discrete data: revival of the
#' Conway–Maxwell–Poisson distribution. Journal of the Royal Statistical
#' Society Series C: Applied Statistics, 54(1), 127-142.
#'
#' @seealso \link{dCOMPO}.
#'
#' @details
#' The COMPO distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has a support 0, 1, 2, ... and mass function given by
#'
#' \eqn{f(x | \mu, \sigma) = \frac{\mu^x}{(x!)^{\sigma} Z(\mu, \sigma)} }
#'
#' with \eqn{\mu > 0}, \eqn{\sigma \geq 0} and
#'
#' \eqn{Z(\mu, \sigma)=\sum_{j=0}^{\infty} \frac{\mu^j}{(j!)^\sigma}}.
#'
#' The proposed functions here are based on the functions from
#' the COMPoissonReg package.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a COMPO distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_COMPO.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
COMPO <- function (mu.link="log", sigma.link="log") {

  mstats <- checklink("mu.link", "COMPO",
                      substitute(mu.link), c("log"))
  dstats <- checklink("sigma.link", "COMPO",
                      substitute(sigma.link), c("log"))

  structure(list(family=c("COMPO", "Comway-Maxwell-Poisson"),
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
                   dm   <- gamlss::numeric.deriv(dCOMPO(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dldm
                 },

                 dldd = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dCOMPO(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   dldd
                 },

                 # Segundas derivadas

                 d2ldm2 = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dCOMPO(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   d2ldm2 <- - dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2, -1e-15)
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {
                   dm   <- gamlss::numeric.deriv(dCOMPO(y, mu, sigma, log=TRUE),
                                                 theta="mu",
                                                 delta=0.01)
                   dldm <- as.vector(attr(dm, "gradient"))
                   dd   <- gamlss::numeric.deriv(dCOMPO(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldmdd <- - dldm * dldd
                   d2ldmdd <- ifelse(d2ldmdd < -1e-15, d2ldmdd, -1e-15)
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {
                   dd   <- gamlss::numeric.deriv(dCOMPO(y, mu, sigma, log=TRUE),
                                                 theta="sigma",
                                                 delta=0.01)
                   dldd <- as.vector(attr(dd, "gradient"))
                   d2ldd2 <- - dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2, -1e-15)
                   d2ldd2
                 },

                 G.dev.incr = function(y, mu, sigma, pw = 1, ...) -2*dCOMPO(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pCOMPO", type="Discrete",
                                               ymin = 0, y = y, mu = mu, sigma = sigma)),

                 #mu.initial    = expression(mu    <- rep(1, length(y)) ),
                 #sigma.initial = expression(sigma <- rep(1, length(y)) ),

                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_COMPO(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_COMPO(y)[2], length(y)) ),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma >= 0),

                 y.valid = function(y) all(y >= 0),

                 mean = function(mu, sigma) {mu^(1/sigma) - (sigma-1)/(2*sigma)},
                 variance = function(mu, sigma) {mu^(1/sigma) / sigma}

  ),
  class=c("gamlss.family", "family"))
}
#' Initial values for COMPO
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @return returns a vector with the moments estimations.
#' @keywords internal
#' @export
#' @importFrom nleqslv nleqslv
#' @importFrom stats var
estim_mu_sigma_COMPO <- function(y) {
  # Star aux fun
  fun <- function(theta, y) {
    mu    <- theta[1]
    sigma <- theta[2]
    approx_mean <- mu^(1/sigma) - (sigma-1)/(2*sigma)
    approx_var  <- mu^(1/sigma) / sigma
    z <- numeric(2) # contiene las ecuaciones
    z[1] <- approx_mean - mean(y)
    z[2] <- approx_var  - var(y)
    z
  }
  # End aux fun

  #library(nleqslv)
  res <- nleqslv(x=c(1, 1), fn=fun, method="Newton",
                 control=list(btol=0.01), y=y)
  res$x
}

# estim_mu_sigma_COMPO <- function(y) {
#   fit <- glm.cmp(y ~ 1)
#   res <- exp(fit$opt.res$par)
#   return(res)
# }


