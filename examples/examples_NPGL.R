# Example 1
# Generating some random values with
# known mu and sigma

set.seed(123)
y <- rNPGL(n=100, mu=20, sigma=2)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=NPGL,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu and sigma
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ NPGL
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(1.7 - 2.8 * x1) # Approx 1.35
  sigma <- exp(0.73 + 1 * x2)  # Approx 3.42
  y <- rNPGL(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(1234)
datos <- gendat(n=200)

mod2 <- gamlss(y~x1, sigma.fo=~x2, family=NPGL, data=datos,
               control=gamlss.control(n.cyc=800, trace=FALSE))

summary(mod2)
