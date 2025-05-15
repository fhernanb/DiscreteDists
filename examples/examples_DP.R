# Example 1
# Generating some random values with
# known mu and sigma
set.seed(123)
y <- rDP(n=1000, mu=2.5, sigma=0.4)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=DP,
               control=gamlss.control(n.cyc=500, trace=TRUE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ DP
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(-1.6 + 5 * x1)
  sigma <- exp(1.7 - 5 * x2)
  y <- rDP(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(123)
datos <- gendat(n=1000)

mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=DP, data=datos,
                 control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod2)

# Example 3
y <- rep(0:4, times=c(46, 76, 24, 9, 1))

# Fitting the model
library(gamlss)
mod3 <- gamlss(y~1, family=DP)

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod3, what="mu"))
exp(coef(mod3, what="sigma"))

# Example 4
y <- rep(0:9, times=c(1, 2, 11, 20, 29, 23, 10, 3, 1, 0))

# Fitting the model
library(gamlss)
mod4 <- gamlss(y~1, family=DP)

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod4, what="mu"))
exp(coef(mod4, what="sigma"))

# Example 5
y <- c(1, 3, 3, 6, 7, 7, 10, 12, 14, 15, 18,
       19, 22, 26, 28, 29, 34, 40, 48, 49)

# Fitting the model
library(gamlss)
mod4 <- gamlss(y~1, family=DP)

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod4, what="mu"))
exp(coef(mod4, what="sigma"))


