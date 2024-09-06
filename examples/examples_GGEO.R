# Example 1
# Generating some random values with
# known mu and sigma
set.seed(123)
y <- rGGEO(n=200, mu=0.95, sigma=1.5)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=GGEO,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
inv_logit <- function(x) 1/(1 + exp(-x))

inv_logit(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ GGEO
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- inv_logit(1.7 - 2.8*x1)
  sigma <- exp(0.73 + 1*x2)
  y <- rGGEO(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(78353)
datos <- gendat(n=100)

mod2 <- gamlss(y~x1, sigma.fo=~x2, family=GGEO, data=datos,
               control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod2)

# Example 3
# Number of accidents to 647 women working on H. E. Shells
# for 5 weeks. Taken from Gomez-Deniz (2010) page 411.

y <- rep(x=0:5, times=c(447, 132, 42, 21, 3, 2))

mod3 <- gamlss(y~1, family=GGEO,
               control=gamlss.control(n.cyc=500, trace=TRUE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
inv_logit <- function(x) 1/(1 + exp(-x))

inv_logit(coef(mod3, what="mu"))
exp(coef(mod3, what="sigma"))
