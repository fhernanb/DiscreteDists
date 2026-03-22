# Example 1
# Generating some random values with known mu and sigma
# logit_inv function
logit_inv <- function(x) exp(x) / (exp(x)+1)
y <- rPTRTE(n=100, mu=0.2, sigma=0.5)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=PTRTE)

# Extracting the fitted values for mu and sigma
exp(coef(mod1, what="mu"))
logit_inv(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model
# A function to simulate a data set with Y ~ PTRTE

gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(2 + 1 * x1) # 12 approximately
  sigma <- logit_inv(2 - 2 * x2) # 0.73 approximately
  y <- rPTRTE(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(1234)
dat <- gendat(n=1000)

# Fitting the model
mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=PTRTE, data=dat,
               control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod2)

# Example 3 (Second data set of the article)
# European corn-borer count data reported by McGuire et al. (1957).
# The observed and fitted frequencies are given in Table 11 of
# Erbayram and Akdogan (2025), where the P-TRTE distribution is
# illustrated using this data set.

values <- 0:5
freq <- c(188, 83, 36, 14, 2, 1)

y <- rep(x=values, times=freq)

mod3 <- gamlss(y~1, sigma.fo=~1, family=PTRTE(),
               control=gamlss.control(n.cyc=500, trace=TRUE))

exp(coef(mod3, what="mu"))
logit_inv(coef(mod3, what="sigma"))
