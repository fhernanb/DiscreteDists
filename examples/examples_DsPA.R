# Example 1
# Generating some random values with
# known mu and sigma
y <- rDsPA(n=100, mu=1.2, sigma=0.5)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=DsPA,
               control=gamlss.control(n.cyc=500, trace=TRUE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
inv_logit <- function(x) 1/(1 + exp(-x))

exp(coef(mod1, what="mu"))
inv_logit(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ GGEO
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  x3 <- runif(n)
  x4 <- runif(n)
  mu    <- exp(1 + 1.2*x1 + 0.2*x2)
  sigma <- inv_logit(2 + 1.5*x3 + 1.5*x4)
  y <- rDsPA(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2, x3=x3, x4=x4)
}

set.seed(123)
datos <- gendat(n=100)

mod2 <- gamlss(y~x1+x2, sigma.fo=~x3+x4, family=DsPA, data=datos,
               control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod2)

# Example 3
# failure times for a sample of 15 electronic components in an acceleration life test
# Taken from
# Alghamdi et. al (202) page 8354.

y <- c(1.0, 5.0, 6.0, 11.0, 12.0, 19.0, 20.0, 22.0,
       23.0, 31.0, 37.0, 46.0, 54.0, 60.0, 66.0)

mod3 <- gamlss(y~1, family=DsPA,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
inv_logit <- function(x) 1/(1 + exp(-x))

exp(coef(mod3, what="mu"))
inv_logit(coef(mod3, what="sigma"))

# Example 4
# number of fires in Greece from July 1, 1998 to August 31, 1998.
# Taken from
# Alghamdi et. al (202) page 8354.

y <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
       1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
       2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
       3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4,
       4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
       5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6,
       6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7,
       8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9,
       9, 9, 9, 9, 10, 10, 10, 11, 11,
       11, 11, 12, 12, 12, 12, 12, 12,
       15, 15, 15, 15, 16, 20, 43)

mod4 <- gamlss(y~1, family=DsPA,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
inv_logit <- function(x) 1/(1 + exp(-x))

exp(coef(mod4, what="mu"))
inv_logit(coef(mod4, what="sigma"))

