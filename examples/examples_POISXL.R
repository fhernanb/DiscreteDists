# Example 1
# Generating some random values with
# known mu
y <- rPOISXL(n=1000, mu=1)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=POISXL,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod1, what='mu'))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ POISXL
gendat <- function(n) {
  x1 <- runif(n, min=0.4, max=0.6)
  mu    <- exp(1.21 - 3 * x1) # 0.75 en promedio
  y <- rPOISXL(n=n, mu=mu)
  data.frame(y=y, x1=x1)
}

dat <- gendat(n=1500)

# Fitting the model
mod2 <- NULL
mod2 <- gamlss(y~x1, family=POISXL, data=dat,
               control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod2)
