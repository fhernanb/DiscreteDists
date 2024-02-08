# Example 1
# Generating some random values with
# known mu and sigma
set.seed(123)
y <- rGGEO(n=5000, mu=1.7, sigma=0.3)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=GGEO, 
               control=gamlss.control(n.cyc=500, trace=TRUE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
inv_logit <- function(x) 1/(1 + exp(-x))
exp(coef(mod1, what='mu'))
inv_logit(coef(mod1, what='sigma'))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ GGEO
gendat <- function(n) { #----------------------------
  x1 <- runif(n, min=0.40, max=0.6)
  x2 <- runif(n,min=0.40, max=0.6)
  mu    <- exp(1.7 - 2.8*x1) 
  sigma <- inv_logit(0.73 - 3*x2) 
  y <- rGGEO(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1,x2=x2)
}

set.seed(7835)

datos <- gendat(n=2000)

mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=GGEO, data=datos,
               control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod2)
