#' Set the Right Link Function for Specified Parameter and Distribution
#'
#' This function is used within the distribution family specification of a
#' GAMLSS model to define the right link for each of the parameters of the
#' distribution. This function should not be called by the user unless
#' he/she specify a new distribution family or wishes to change existing
#' link functions in the parameters.
#'
#' @param which.link which parameter link e.g. \code{which.link="mu.link"}
#' @param which.dist which distribution family e.g. \code{which.dist="Cole.Green"}
#' @param link a repetition of \code{which.link} e.g. \code{link=substitute(mu.link)}
#' @param link.List what link function are required e.g. \code{link.List=c("inverse", "log", "identity")}
#'
#' @details This function is the same \code{checklink} fuction from gamlss.dist package.
#'
#' @return Defines the right link for each parameter.
#'
#' @export
checklink <- function (which.link = NULL, which.dist = NULL,
                       link = NULL, link.List = NULL) {
  if (is.null(which.link))
    stop(paste("The parameter link name has not been defined."))
  if (is.null(which.dist))
    stop(paste("The distribution has not been defined."))
  if (is.null(link))
    stop(paste("The link has not been defined."))
  if (is.null(link.List))
    stop(paste("The list of links has not been defined."))
  linktemp <- link
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  }
  if (linktemp %in% link.List)
    stats <- make.link.gamlss(linktemp)
  else if (is.character(link)) {
    stats <- make.link.gamlss(link)
    linktemp <- link
  }
  else {
    if (inherits(eval(link), "link-gamlss")) {
      stats <- eval(link)
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else {
      stop(gettextf("\"%s link \"%s\" not available, for \"%s family; available links are %s",
                    which.link, linktemp, which.dist, paste(sQuote(link.List),
                                                            collapse = ", ")), domain = NA)
    }
  }
  link.result <- c(linktemp, stats)
  link.result
}
