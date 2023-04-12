# Example 1
# Generating some random values with
# known mu and sigma
y <- rHYPERPO2(n=100, mu=5, sigma=1.5)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=HYPERPO2,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what='mu'))
exp(coef(mod1, what='sigma'))

# Example 2
# Generating random values under some model

\dontrun{
  # A function to simulate a data set with Y ~ HYPERPO
  gendat <- function(n) {
    x1 <- runif(n)
    x2 <- runif(n)
    mu    <- exp(1.21 - 3 * x1) # 0.75 en promedio
    sigma <- exp(1.26 - 2 * x2) # 1.30 en promedio
    y <- rHYPERPO2(n=n, mu=mu, sigma=sigma)
    data.frame(y=y, x1=x1,x2=x2)
  }

  set.seed(1235)
  datos <- gendat(n=150)

  mod2 <- NULL
  mod2 <- gamlss(y~x1, sigma.fo=~x2, family=HYPERPO2, data=datos,
                 control=gamlss.control(n.cyc=500, trace=TRUE))

  summary(mod2)
}
