# Example 1
# Generating some random values with
# known mu and sigma
set.seed(150)
y <- rDIKUM(1000, mu=1, sigma=5)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y ~ 1, sigma.fo = ~1, family=DIKUM,
               control = gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model

library(gamlss)

# A function to simulate a data set with Y ~ DIKUM
gendat <- function(n) {
  x1 <- runif(n, min=0.4, max=0.6)
  x2 <- runif(n, min=0.4, max=0.6)
  mu    <- exp(1.21 - 3 * x1) # 0.75 approximately
  sigma <- exp(1.26 - 2 * x2) # 1.30 approximately
  y <- rDIKUM(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

dat <- gendat(n=1000)

# Fitting the model
mod2 <- gamlss(y ~ x1, sigma.fo = ~x2, family = "DIKUM", data=dat,
               control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod2)
