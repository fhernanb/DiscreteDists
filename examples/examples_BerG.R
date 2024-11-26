# Example 1
# Generating some random values with
# known mu and sigma
y <- rBerG(n=500, mu=0.75, sigma=0.5)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=BerG,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu and sigma
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ GGEO
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  x3 <- runif(n)
  x4 <- runif(n)
  mu    <- exp(1 + 1.2*x1 + 0.2*x2)
  sigma <- exp(2 + 1.5*x3 + 1.5*x4)
  y <- rBerG(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2, x3=x3, x4=x4)
}

set.seed(16494786)
datos <- gendat(n=500)

mod2 <- gamlss(y~x1+x2, sigma.fo=~x3+x4, family=BerG, data=datos,
               control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod2)

# Example using the dataset grazing from the bergreg package
# https://github.com/rdmatheus/bergreg

# This example corresponds to example 5.1
# presented by Bourguignon & Medeiros (2022)
# A simple and useful regression model for fitting count data

data("grazing")
hist(grazing$birds)

mod3 <- gamlss(birds ~ when + grazed,
               sigma.fo=~1,
               family=BerG, data=grazing,
               control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod3)

