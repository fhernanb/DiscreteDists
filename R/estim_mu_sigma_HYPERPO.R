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
