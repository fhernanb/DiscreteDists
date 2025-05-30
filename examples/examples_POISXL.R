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
exp(coef(mod1, what="mu"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ POISXL
gendat <- function(n) {
  x1 <- runif(n, min=0.4, max=0.6)
  mu <- exp(1.21 - 3 * x1) # 0.75 approximately
  y <- rPOISXL(n=n, mu=mu)
  data.frame(y=y, x1=x1)
}

dat <- gendat(n=1500)

# Fitting the model
mod2 <- NULL
mod2 <- gamlss(y~x1, family=POISXL, data=dat,
               control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod2)

# Example 3
# The counts the number of borers per hill of corn in an
# experiment conducted randomly on 8 hills in 15 replications.
# Taken from Ahsan-ul-Haq et al (2022) page 10.

y <- rep(x=0:8, times=c(43, 35, 17, 11, 5, 4, 1, 2, 2))

mod3 <- gamlss(y~1, family=POISXL,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
exp(coef(mod3, what="mu"))

# Example 4
# The number of mammalian cytogenetic dosimetry lesions produced
# by streptogramin exposure in rabbit.
# Taken from Ahsan-ul-Haq et al (2022) page 10.

y <- rep(x=0:4, times=c(200, 57, 30, 7, 6))

mod4 <- gamlss(y~1, family=POISXL,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
exp(coef(mod4, what="mu"))


