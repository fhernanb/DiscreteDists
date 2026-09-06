#' Create a Link for GAMLSS families
#'
#' The function \code{make.link.gamlss()} is used with \code{gamlss.family}
#' distributions in package \pkg{gamlss()}.Given a link, it returns a link
#' function, an inverse link  function, the derivative dpar/deta where 'par'
#' is the appropriate distribution parameter and a function for checking
#' the domain. It differs from the usual \code{make.link} of \code{glm()} by
#' having extra links as the \code{logshifto1}, and the \code{own}.
#' For the use of the \code{own} link see the example bellow.
#' \code{show.link} provides a way in which the user can identify the link
#' functions available for each gamlss distribution. If your required link
#' function is not available for any of the gamlss distributions you can add it in.
#'
#' @param link character or numeric; one of \code{"logit"}, \code{"probit"},
#' \code{"cloglog"}, \code{"identity"}, \code{"log"},  \code{"sqrt"},
#' \code{"1/mu^2"}, \code{"inverse"}, \code{"logshifted"},
#' \code{"logitshifted"}, or number, say lambda resulting in power link
#' \eqn{\mu^\lambda}{mu^lambda}.
#'
#' @details The \code{own} link function is added to allow the user greater
#' flexibility. In order to used the own link function for any of the
#' parameters of the distribution the \code{own} link should  appear in the
#' available links for this parameter. You can check this using the function
#' \code{show.link}. If the \code{own} do not appear in the list you can
#' create a new function for the distribution in which \code{own} is added
#' in the list. For example the first line of the code  of the binomial
#' distribution, \code{BI}, has change from
#'
#' "mstats <- checklink("mu.link", "Binomial", substitute(mu.link),
#' c("logit", "probit", "cloglog", "log")), in version 1.0-0 of gamlss, to
#'
#' "mstats <- checklink("mu.link", "Binomial", substitute(mu.link),
#' c("logit", "probit", "cloglog", "log", "own"))
#'
#' in version 1.0-1. Given that the parameter has \code{own} as an option the user needs also
#' to define the following four new functions in order to used an
#' \code{own} link.
#'
#' i) own.linkfun
#'
#' ii) own.linkinv
#'
#' iii) own.mu.eta and
#'
#' iv) own.valideta.
#'
#' An example is given below.
#'
#' Only one parameter of the distribution at a time is allowed to have
#' its \code{own} link, (unless the same four \code{own} functions above
#' are suitable for more that one parameter of the distribution).
#'
#' Note that from \pkg{gamlss} version 1.9-0 the user can introduce its own
#' link function by define an appropriate function, (see the example below).
#'
#' @return
#' For the \code{make.link.gamlss} a list with components
#'
#' linkfun: Link function \code{function(parameter)}
#'
#' linkinv: Inverse link function \code{function(eta)}
#'
#' mu.eta: Derivative \code{function(eta)} dparameter/deta
#'
#' valideta: \code{function(eta)} TRUE if all of eta is in the domain of \code{linkinv}.
#'
#' For the \code{show.link} a list with components the available links for the distribution parameters
#'
#' @note
#' For the links  involving parameters as in \code{logshifted} and
#' \code{logitshifted} the parameters can be passed in the definition of the
#' distribution by calling the  \code{ checklink} function, for example in the
#' definition of the \code{tau} parameter in BCPE distribution the following
#' call is made: \code{tstats <- checklink("tau.link",
#' "Box Cox Power Exponential", substitute(tau.link),
#' c("logshifted", "log", "identity"), par.link = c(1))  }
#'
#' @note This function is the same \code{make.link.gamlss} fuction from gamlss.dist package.
#'
#' @importFrom stats power qlogis plogis dlogis qnorm qcauchy pcauchy dcauchy
#' @export
make.link.gamlss <- function (link) {
  if (is.character(link) && length(grep("^power", link) > 0)) {
    warning("calling make.link(\"power(z)\") is deprecated",
            domain = NA)
    return(eval(parse(text = link)))
  }
  else if (!is.character(link) && !is.na(lambda <- as.numeric(link))) {
    warning("calling make.link(number) is deprecated", domain = NA)
    return(power(lambda))
  }
  else switch(link, logit = {
    linkfun <- function(mu) qlogis(mu)
    linkinv <- function(eta) {
      thresh <- -qlogis(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      plogis(eta)
    }
    mu.eta <- function(eta) pmax(dlogis(eta), .Machine$double.eps)
    valideta <- function(eta) TRUE
  }, probit = {
    linkfun <- function(mu) qnorm(mu)
    linkinv <- function(eta) {
      thresh <- -qnorm(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      pnorm(eta)
    }
    mu.eta <- function(eta) pmax(dnorm(eta), .Machine$double.eps)
    valideta <- function(eta) TRUE
  }, cauchit = {
    linkfun <- function(mu) qcauchy(mu)
    linkinv <- function(eta) {
      thresh <- -qcauchy(.Machine$double.eps)
      eta <- pmin(pmax(eta, -thresh), thresh)
      pcauchy(eta)
    }
    mu.eta <- function(eta) pmax(dcauchy(eta), .Machine$double.eps)
    valideta <- function(eta) TRUE
  }, cloglog = {
    linkfun <- function(mu) log(-log(1 - mu))
    linkinv <- function(eta) pmax(pmin(-expm1(-exp(eta)),
                                       1 - .Machine$double.eps), .Machine$double.eps)
    mu.eta <- function(eta) {
      eta <- pmin(eta, 700)
      pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
    }
    valideta <- function(eta) TRUE
  }, identity = {
    linkfun <- function(mu) mu
    linkinv <- function(eta) eta
    mu.eta <- function(eta) rep(1, length(eta))
    valideta <- function(eta) TRUE
  }, log = {
    linkfun <- function(mu) log(mu)
    linkinv <- function(eta) pmax(exp(eta), .Machine$double.eps)
    mu.eta <- function(eta) pmax(exp(eta), .Machine$double.eps)
    valideta <- function(eta) TRUE
  }, sqrt = {
    linkfun <- function(mu) mu^0.5
    linkinv <- function(eta) eta^2
    mu.eta <- function(eta) 2 * eta
    valideta <- function(eta) all(eta > 0)
  }, `1/mu^2` = {
    linkfun <- function(mu) 1/mu^2
    linkinv <- function(eta) 1/eta^0.5
    mu.eta <- function(eta) -1/(2 * eta^1.5)
    valideta <- function(eta) all(eta > 0)
  }, `mu^2` = {
    linkfun <- function(mu) mu^2
    linkinv <- function(eta) eta^0.5
    mu.eta <- function(eta) 0.5 * (eta^-0.5)
    valideta <- function(eta) all(eta > 0)
  }, logshiftto1 = {
    linkfun <- function(mu) log(mu - 1 + 1e-05)
    linkinv <- function(eta) 1 + pmax(.Machine$double.eps,
                                      exp(eta))
    mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
    valideta <- function(eta) TRUE
  }, logshiftto2 = {
    linkfun <- function(mu) log(mu - 2 + 1e-05)
    linkinv <- function(eta) 2 + pmax(.Machine$double.eps,
                                      exp(eta))
    mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
    valideta <- function(eta) TRUE
  }, logshiftto0 = {
    linkfun <- function(mu) {
      log(mu - 1e-05)
    }
    linkinv <- function(eta) {
      1e-05 + pmax(.Machine$double.eps, exp(eta))
    }
    mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
    valideta <- function(eta) TRUE
  }, Slog = {
    linkfun <- function(mu) {
      log(mu - 1e-05)
    }
    linkinv <- function(eta) {
      1e-05 + pmax(.Machine$double.eps, exp(eta))
    }
    mu.eta <- function(eta) pmax(.Machine$double.eps, exp(eta))
    valideta <- function(eta) TRUE
  }, `[-1,1]` = {
    linkfun <- function(mu) {
      delta <- 1e-10
      shift <- c(-1 - delta, 1 + delta)
      log((mu - shift[1])/(shift[2] - mu))
    }
    linkinv <- function(eta) {
      delta <- 1e-10
      shift <- c(-1 - delta, 1 + delta)
      thresh <- -log(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      (shift[2] * exp(eta) + shift[1])/(1 + exp(eta))
    }
    mu.eta <- function(eta) {
      delta <- 1e-10
      shift <- c(-1 - delta, 1 + delta)
      thresh <- -log(.Machine$double.eps)
      res <- rep(.Machine$souble.eps, length(eta))
      res[abs(eta) < thresh] <- (shift[2] * exp(eta))/(1 +
                                                         exp(eta))[abs(eta) < thresh] - (exp(eta) * (shift[2] *
                                                                                                       exp(eta) + shift[1]))/((1 + exp(eta))^2)[abs(eta) <
                                                                                                                                                  thresh]
      res
    }
    valideta <- function(eta) TRUE
  }, `(0,2]` = {
    linkfun <- function(mu) {
      delta <- 1e-10
      shift <- c(0, 2 + delta)
      log((mu - shift[1])/(shift[2] - mu))
    }
    linkinv <- function(eta) {
      delta <- 1e-10
      shift <- c(0, 2 + delta)
      thresh <- -log(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      (shift[2] * exp(eta) + shift[1])/(1 + exp(eta))
    }
    mu.eta <- function(eta) {
      delta <- 1e-10
      shift <- c(0, 2 + delta)
      thresh <- -log(.Machine$double.eps)
      res <- rep(.Machine$souble.eps, length(eta))
      res[abs(eta) < thresh] <- (shift[2] * exp(eta))/(1 +
                                                         exp(eta))[abs(eta) < thresh] - (exp(eta) * (shift[2] *
                                                                                                       exp(eta) + shift[1]))/((1 + exp(eta))^2)[abs(eta) <
                                                                                                                                                  thresh]
      res
    }
    valideta <- function(eta) TRUE
  }, `(0,5]` = {
    linkfun <- function(mu) {
      delta <- 1e-10
      shift <- c(0, 5 + delta)
      log((mu - shift[1])/(shift[2] - mu))
    }
    linkinv <- function(eta) {
      delta <- 1e-10
      shift <- c(0, 5 + delta)
      thresh <- -log(.Machine$double.eps)
      eta <- pmin(thresh, pmax(eta, -thresh))
      (shift[2] * exp(eta) + shift[1])/(1 + exp(eta))
    }
    mu.eta <- function(eta) {
      delta <- 1e-10
      shift <- c(0, 5 + delta)
      thresh <- -log(.Machine$double.eps)
      res <- rep(.Machine$souble.eps, length(eta))
      res[abs(eta) < thresh] <- (shift[2] * exp(eta))/(1 +
                                                         exp(eta))[abs(eta) < thresh] - (exp(eta) * (shift[2] *
                                                                                                       exp(eta) + shift[1]))/((1 + exp(eta))^2)[abs(eta) <
                                                                                                                                                  thresh]
      res
    }
    valideta <- function(eta) TRUE
  }, own = {
    linkfun <- function(mu) eval(body(get("own.linkfun",
                                          envir = globalenv())))
    linkinv <- function(eta) eval(body(get("own.linkinv",
                                           envir = globalenv())))
    mu.eta <- function(eta) eval(body(get("own.mu.eta",
                                          envir = globalenv())))
    valideta <- function(eta) eval(body(get("own.valideta",
                                            envir = globalenv())))
  }, inverse = {
    linkfun <- function(mu) 1/mu
    linkinv <- function(eta) 1/eta
    mu.eta <- function(eta) -1/(eta^2)
    valideta <- function(eta) all(eta != 0)
  }, stop(sQuote(link), " link not recognised"))
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
                 valideta = valideta, name = link), class = "link-gamlss")
}
