# Example 1
# Generating some random values with
# known mu
y <- rDBH(n=1000, mu=0.74)

library(gamlss)
mod1 <- gamlss(y~1, family=DBH,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse logit function
inv_logit <- function(x) exp(x) / (1+exp(x))
inv_logit(coef(mod1, parameter="mu"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ DBH
gendat <- function(n) {
  x1 <- runif(n)
  mu    <- inv_logit(-3 + 5 * x1)
  y <- rDBH(n=n, mu=mu)
  data.frame(y=y, x1=x1)
}

datos <- gendat(n=150)

mod2 <- NULL
mod2 <- gamlss(y~x1, family=DBH, data=datos,
               control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod2)

# Example 3
# Number of carious teeth among the four deciduous molars.
# Taken from EL-MORSHEDY (2020) page 74364.

y <- rep(0:4, times=c(64, 17, 10, 6, 3))

mod3 <- gamlss(y~1, family=DBH,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
inv_logit <- function(x) 1/(1 + exp(-x))
inv_logit(coef(mod3, what="mu"))

# Example 4
# Counts of cysts of kidneys using steroids.
# Taken from EL-MORSHEDY (2020) page 74365.

y <- rep(0:11, times=c(65, 14, 10, 6, 4, 2, 2, 2, 1, 1, 1, 2))

mod4 <- gamlss(y~1, family=DBH,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
inv_logit <- function(x) 1/(1 + exp(-x))
inv_logit(coef(mod4, what="mu"))

