# Example 1
# Generating some random values with
# known mu
y <- rPOISXL(n=1000, mu=1)

estim_mu_POISXL(y)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=POISXL,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod1, what='mu'))

