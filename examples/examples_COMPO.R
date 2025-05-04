# Example 1
# Generating some random values with
# known mu and sigma
set.seed(1234)
y <- rCOMPO(n=100, mu=5, sigma=1.5)

estim_mu_sigma_COMPO(y)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=COMPO,
               control=gamlss.control(n.cyc=500, trace=TRUE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model

\donttest{
# A function to simulate a data set with Y ~ COMPO
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(2 + 1 * x1) # 12 approximately
  sigma <- exp(2 - 2 * x2) # 2.71 approximately
  y <- rCOMPO(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(123)
dat <- gendat(n=100)

# Fitting the model
mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=COMPO, data=dat,
               control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod2)
}

