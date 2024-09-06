# Example 1
# Generating some random values with
# known mu
y <- rDLD(n=100, mu=0.3)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=DLD,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod1, what='mu'))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ DLD
gendat <- function(n) {
  x1 <- runif(n)
  mu    <- exp(2 - 4 * x1)
  y <- rDLD(n=n, mu=mu)
  data.frame(y=y, x1=x1)
}

set.seed(1235)
datos <- gendat(n=150)

mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=DLD, data=datos,
                 control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod2)
