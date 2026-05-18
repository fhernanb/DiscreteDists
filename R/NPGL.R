#' Poisson-generalised Lindley distribution
#'
#' @author Tomás Mesa, \email{tomas.mesaz@udea.edu.co}
#'
#' @description
#' The function \code{NPGL()} defines the
#' Poisson-generalised Lindley distribution,
#' a two parameter distribution,
#' for a \code{gamlss.family} object to be used in GAMLSS
#' fitting using the function \code{gamlss()}.
#'
#' @param mu.link defines the mu.link, with "log" link as the default for the mu parameter.
#' @param sigma.link defines the sigma.link, with "log" link as the default for the sigma.
#'
#' @references
#' Altun, E.
#' A new two-parameter discrete poisson-generalized Lindley distribution with
#' properties and applications to healthcare data sets. Comput Stat 36, 2841–2861 (2021). https://doi.org/10.1007/s00180-021-01097-0
#'
#' @seealso \link{dNPGL}.
#'
#' @details
#' The Poisson-generalised Lindley distribution with parameters \eqn{\mu} and \eqn{\sigma}
#' has support \eqn{x = 0, 1, 2, \ldots} and probability mass function given by
#'
#' \eqn{f(x \mid \mu, \sigma)=\frac{\mu^2+\frac{\mu^{\sigma}(\mu+1)^{1-\sigma}\Gamma(x+\sigma)}{\Gamma(\sigma)\Gamma(x+1)}}{(\mu+1)^{x+2}}}
#'
#' with \eqn{\mu > 0} and \eqn{\sigma > 0}.
#'
#' This distribution is useful for modeling over-dispersed count data.
#'
#' Note: in this implementation we changed the original parameters \eqn{\theta} and \eqn{\alpha}
#' for \eqn{\mu} and \eqn{\sigma} respectively, we did it to implement this distribution within gamlss framework.
#'
#' @return
#' Returns a \code{gamlss.family} object which can be used
#' to fit a NPGL distribution
#' in the \code{gamlss()} function.
#'
#' @example examples/examples_NPGL.R
#'
#' @importFrom gamlss.dist checklink
#' @importFrom gamlss rqres.plot
#' @export
NPGL <- function (mu.link="log", sigma.link="log") {

  mstats <- checklink("mu.link", "Poisson-generalised Lindley",
                      substitute(mu.link), c("log", "inverse", "own"))
  dstats <- checklink("sigma.link", "Poisson-generalised Lindley",
                      substitute(sigma.link), c("log", "inverse", "own"))

  structure(list(family=c("NPGL", "Poisson-generalised Lindley"),
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

                 # First derivatives

                 dldm = function(y, mu, sigma) {

                   A <- (mu+1)*(mu^sigma)*gamma(y+sigma)
                   B <- (mu^2)*((mu+1)^sigma)*(gamma(sigma))*gamma(y+1)

                   part1 <- A*(mu-sigma+mu*y)+B*(mu*y-2)
                   part2 <- (mu*(mu+1))*(A+B)

                   dldm <- -(part1/part2)
                   dldm
                 },

                 dldd = function(y, mu, sigma) {

                   A <- (mu+1)*(mu^sigma)*gamma(y+sigma)
                   B <- (mu^2)*((mu+1)^sigma)*(gamma(sigma))*gamma(y+1)

                   part1 <- A*(log(mu)-log(mu+1)-digamma(sigma)+digamma(y+sigma))
                   part2 <- A+B

                   dldd <- part1/part2
                   dldd
                 },

                 # Second derivatives

                 d2ldm2 = function(y, mu, sigma) {

                   A <- (mu+1)*(mu^sigma)*gamma(y+sigma)
                   B <- (mu^2)*((mu+1)^sigma)*(gamma(sigma))*gamma(y+1)

                   part1 <- A*(mu-sigma+mu*y)+B*(mu*y-2)
                   part2 <- (mu*(mu+1))*(A+B)

                   dldm <- -(part1/part2)

                   d2ldm2 <- - dldm * dldm
                   d2ldm2
                 },

                 d2ldmdd = function(y, mu, sigma) {

                   A <- (mu+1)*(mu^sigma)*gamma(y+sigma)
                   B <- (mu^2)*((mu+1)^sigma)*(gamma(sigma))*gamma(y+1)

                   part1 <- A*(mu-sigma+mu*y)+B*(mu*y-2)
                   part2 <- (mu*(mu+1))*(A+B)

                   dldm <- -(part1/part2)

                   A <- (mu+1)*(mu^sigma)*gamma(y+sigma)
                   B <- (mu^2)*((mu+1)^sigma)*(gamma(sigma))*gamma(y+1)

                   part1 <- A*(log(mu)-log(mu+1)-digamma(sigma)+digamma(y+sigma))
                   part2 <- A+B

                   dldd <- part1/part2

                   d2ldmdd <- - dldm * dldd
                   d2ldmdd
                 },

                 d2ldd2  = function(y, mu, sigma) {

                   A <- (mu+1)*(mu^sigma)*gamma(y+sigma)
                   B <- (mu^2)*((mu+1)^sigma)*(gamma(sigma))*gamma(y+1)

                   part1 <- A*(log(mu)-log(mu+1)-digamma(sigma)+digamma(y+sigma))
                   part2 <- A+B

                   dldd <- part1/part2

                   d2ldd2 <- - dldd * dldd
                   d2ldd2
                 },


                 G.dev.incr = function(y, mu, sigma, pw=1, ...) -2*dNPGL(y, mu, sigma, log=TRUE),
                 rqres      = expression(rqres(pfun="pNPGL", type="Discrete",
                                               ymin=0, y=y, mu=mu, sigma=sigma)),

                 mu.initial    = expression(mu    <- rep(estim_mu_sigma_NPGL(y)[1], length(y)) ),
                 sigma.initial = expression(sigma <- rep(estim_mu_sigma_NPGL(y)[2], length(y)) ),

                 mu.valid    = function(mu)    all(mu > 0),
                 sigma.valid = function(sigma) all(sigma > 0),

                 y.valid = function(y) all(y >= 0),

                 mean = function(mu, sigma) {
                   (sigma + mu) / (mu * (1 + mu))
                 },
                 variance = function(mu, sigma) {
                   (mu * (sigma^2 + mu*(sigma + mu + 2) + 2) + sigma) / (mu^2 * (1 + mu)^2)
                 }
  ),
  class=c("gamlss.family", "family"))
}
#'
#' estim_mu_sigma_NPGL
#'
#' This function generates initial values for NPGL distribution.
#'
#' @param y vector with the random sample
#' @examples
#' y <- rNPGL(n=100, mu=0.1, sigma=2)
#' estim_mu_sigma_NPGL(y=y)
#' @importFrom stats optim
#' @export
estim_mu_sigma_NPGL <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_NPGL,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  #res <- c(0, 1, 1, 1) # esto se lo puse para que no tenga en cuenta optim
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#'
#' logLik_NPGL
#'
#' This is an auxiliar function to obtain the logLik for NPGL.
#'
#' @param param vector with the values for mu and sigma
#' @param x vector with the data
#' @examples
#' y <- rNPGL(n=100, mu=0.1, sigma=2)
#' logLik_NPGL(param=c(0, 0), x=y)
#' @importFrom stats optim
#' @export
logLik_NPGL <- function(param=c(0, 0), x){
  return(sum(dNPGL(x,
                  mu    = exp(param[1]),
                  sigma = exp(param[2]),
                  log=TRUE)))
}

################################################################################
#'
#' estim_mu_sigma_NPGL_MM
#'
#' This function generates initial values for NPGL distribution
#' using the method of moments estimators given in equations (19) and (20)
#' of Altun (2021).
#'
#' @param y vector with the random sample
#' @examples
#' y <- rNPGL(n=100, mu=0.1, sigma=2)
#' estim_mu_sigma_NPGL_MM(y=y)
#' @export
estim_mu_sigma_NPGL_MM <- function(y) {
  m1 <- mean(y)
  m2 <- mean(y^2)

  m2 <- var(y) + m1^2

  discriminant <- m1^4 - 6*m1^3 - m1^2*(2*m2 + 3) + 2*m1*m2 + m2^2 + m1 + m2

  mu_hat <- (sqrt(discriminant) - m1^2) / (2 * m1^2)
  sigma_hat <- m1 * mu_hat * (1 + mu_hat) - mu_hat

  res <- c(mu_hat    = mu_hat,
           sigma_hat = sigma_hat)
  return(res)
}
