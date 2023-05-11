#' logLik function for hyper Poisson
#' @description Calculates logLik for hyper Poisson distribution.
#' @param logparam vector with parameters in log scale.
#' @param x vector with the response variable.
#' @keywords internal
#' @export
logLik_HYPERPO <- function(logparam=c(0, 0), x){
  return(sum(dHYPERPO(x     = x,
                      mu    = exp(logparam[1]),
                      sigma = exp(logparam[2]),
                      log=TRUE)))
}
#' Initial values for hyper Poisson
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @keywords internal
#' @export
#' @importFrom stats optim
estim_mu_sigma_HYPERPO <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_HYPERPO,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#' Auxiliar function for hyper Poisson
#' @description This function is used inside density function of Hyper Poisson.
#' @param c,z values for F11.
#' @keywords internal
#' @export
F11 <- function(c, z) {
  k <- 0:999
  res <- lgamma(c) + k*log(z) - lgamma(c+k)
  res <- exp(res)
  sum(res)
}
F11 <- Vectorize(F11)
#' Auxiliar function for hyper Poisson
#' @description This function is used to calculate (a)r.
#' @param a first value.
#' @param r second value.
#' @keywords internal
#' @export
AR <- function(a, r) {
  res <- gamma(a+r) / gamma(a)
  res
}
#' logLik function for hyper Poisson in second parameterization
#' @description Calculates logLik for hyper Poisson distribution.
#' @param logparam vector with parameters in log scale.
#' @param x vector with the response variable.
#' @keywords internal
#' @export
logLik_HYPERPO2 <- function(logparam=c(0, 0), x){
  return(sum(dHYPERPO2(x     = x,
                       mu    = exp(logparam[1]),
                       sigma = exp(logparam[2]),
                       log=TRUE)))
}
#' Initial values for hyper Poisson in second parameterization
#' @description This function generates initial values for the parameters.
#' @param y vector with the response variable.
#' @keywords internal
#' @export
#' @importFrom stats optim
estim_mu_sigma_HYPERPO2 <- function(y) {
  mod <- optim(par=c(0, 0),
               fn=logLik_HYPERPO2,
               method="Nelder-Mead",
               control=list(fnscale=-1, maxit=100000),
               x=y)
  res <- c(mu_hat    = exp(mod$par[1]),
           sigma_hat = exp(mod$par[2]))
  names(res) <- c("mu_hat", "sigma_hat")
  return(res)
}
#' Auxiliar function to obtain lambda from E(X)
#' @description This function implements the procedure given in page 150.
#' @param media the value for the mean or E(X).
#' @param gamma the value for the gamma parameter.
#' @keywords internal
#' @export
obtaining_lambda <- function(media, gamma) {
  # Begin aux function
  fun <- function(x) x-(gamma-1)*(F11(gamma,x)-1)/F11(gamma,x)-media
  fun <- Vectorize(fun)
  # End aux function
  if (gamma == 1)
    result <- media
  else {
    res <- uniroot(f=fun,
                   lower=min(media, max(media+gamma-1, gamma*media)),
                   upper=max(media, min(media+gamma-1, gamma*media)))
    result <- res$root
  }
  result
}
obtaining_lambda <- Vectorize(obtaining_lambda)

#' logLik function for Discrete Burr Hatke
#' @description Calculates logLik for Discrete Burr Hatke  distribution.
#' @param param value for mu.
#' @param x vector with the response variable.
#' @keywords internal
#' @export
logLik_DBH <- function(param=0.5, x){
  mu <- param
  minus_ll <- -sum(dDBH(x=x, mu=mu, log=TRUE))
  return(minus_ll)
}
#' Initial values for Discrete Burr Hatke
#' @description This function generates initial values for the parameter mu.
#' @param y vector with the response variable.
#' @keywords internal
#' @export
#' @importFrom stats optimize
estim_mu_DBH <- function(y){
  mod <- optimize(f=logLik_DBH, interval=c(0, 1), x=y)
  res <- mod$minimum
  return(res)
}
