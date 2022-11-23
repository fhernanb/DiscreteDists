# Example 1
# Generating some random values with
# known mu and sigma
y <- rHYPERPO(n=300, mu=5, sigma=1.5)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=HYPERPO,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what='mu'))
exp(coef(mod1, what='sigma'))

# Example 2
# Generating random values under some model

#\dontrun{
n <- 200
x1 <- runif(n)
x2 <- runif(n)
mu <- exp(1.21 - 3 * x1)
sigma <- exp(1.26 - 2 * x2)
x <- rHYPERPO(n=n, mu, sigma)

mod2 <- gamlss(x~x1, sigma.fo=~x2, family=HYPERPO,
              control=gamlss.control(n.cyc=5000, trace=FALSE))

coef(mod2, what="mu")
coef(mod2, what="sigma")
#}
