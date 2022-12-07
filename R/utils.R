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
#' @param a,c,z values for F11.
#' @keywords internal
#' @export
F11 <- function(a, c, z) {
  r <- 0:99
  num <- AR(a=a, r=r) * z^r
  den <- AR(a=c, r=r) * factorial(r)
  sum(num / den)
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
  # Aux function
  fun <- function(x) x-(gamma-1)*(F11(1,gamma,x)-1)/F11(1,gamma,x)-media
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
