# Example 1
# Generating some random values with
# known mu
y <- rDBH(n=1000, mu=0.74)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=DBH,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse logit function
inv_logit <- function(x) exp(x) / (1+exp(x))
inv_logit(coef(mod1, parameter = 'mu'))

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
               control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod2)
