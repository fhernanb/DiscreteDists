# Example 1
# Generating some random values with
# known mu and sigma
y <- rCOMPO2(n=500, mu=10, sigma=-1)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=COMPO2,
               control=gamlss.control(n.cyc=500, trace=TRUE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what="mu"))
coef(mod1, what="sigma")

# Example 2
# Generating random values under some model

\dontrun{
# A function to simulate a data set with Y ~ COMPO2
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(2 + 1 * x1) # 12 approximately
  sigma <- 1 - 2 * x2      # 0 approximately
  y <- rCOMPO2(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(123)
dat <- gendat(n=200)

# Fitting the model
mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=COMPO2, data=dat)

summary(mod2)
}
