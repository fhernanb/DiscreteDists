#' Drae the CDF for a discrete random variables
#'
#' @author Freddy Hernandez, \email{fhernanb@unal.edu.co}
#'
#' @param x vector with the values of the random variable.
#' @param fx vector with the probabilities of x.
#' @param ... further arguments and graphical parameters.
#'
#' @return
#' A plot with the cumulative density function.
#'
#' @example  examples/examples_plot_cdf.R
#' @importFrom graphics abline points
#' @export
#'
plot_discrete_cdf <- function(x, fx, ...) {
  Fx <- cumsum(fx)
  n <- length(x)

  plot(x = NA, y = NA, pch = NA,
       xlim = c(0, max(x)),
       ylim = c(0, 1), ...)

  abline(h=c(0, 1), col="grey", lty=2) # horizontal lines

  points(x=x[-n], y=Fx[-1], pch=19)
  points(x=x[-1], y=Fx[-1], pch=1)
  for(i in 1:(n-1))
    points(x=x[i+0:1], y=Fx[c(i,i)+1], type="l")
}
