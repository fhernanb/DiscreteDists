# Example 1
# Generating some random values with
# known mu and sigma
set.seed(1234)
y <- rDMOLBE(n=100, mu=10, sigma=7)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=DMOLBE,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ DMOLBE
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(1.21 - 3 * x1) # 0.75 approximately
  sigma <- exp(1.26 - 2 * x2) # 1.30 approximately
  y <- rDMOLBE(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(123)
dat <- gendat(n=200)

# Fitting the model
mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=DMOLBE, data=dat,
                 control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod2)


# Example 3
# Data Set I (death due to coronavirus in China). The first data set is the number
# of deaths due to coronavirus in China from 23 January to 28 March.
# The data sets used in the paper was collected from 2020 year. The data set
# is reported in https://www.worldometers.info/coronavirus/country/china/.
# The data are:

y <- c(8, 16, 15, 24, 26, 26, 38, 43, 46, 45, 57, 64, 65, 73, 73, 86, 89, 97,
       108, 97, 146, 121, 143, 142, 105, 98, 136, 114, 118, 109, 97, 150, 71,
       52, 29, 44, 47, 35, 42, 31, 38, 31, 30, 28, 27, 22, 17, 22, 11, 7,
       13, 10, 14, 13, 11, 8, 3, 7, 6, 9, 7, 4, 6, 5, 3, 5)

# Fitting the model
mod3 <- gamlss(y~1, sigma.fo=~1, family=DMOLBE,
               control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod3)

# Extracting the fitted values for mu and sigma
# using the inverse link function
mu_hat <- exp(coef(mod3, what="mu"))
mu_hat
sigma_hat <- exp(coef(mod3, what="sigma"))
sigma_hat

# Example 4
# Data Set II (daily death due to coronavirus in Pakistan). The second data
# set is the daily deaths due to coronavirus in Pakistan from 18 March
# to 30 June. The data sets used in the paper was collected from 2020 year.
# The data is reported in
# https://www.worldometers.info/coronavirus/country/Pakistan.
# The data are:

y <- c(1, 6, 6, 4, 4, 4, 1, 20, 5, 2, 3, 15, 17, 7, 8, 25, 8, 25, 11,
       25, 16, 16, 12, 11, 20, 31, 42, 32, 23, 17, 19, 38, 50, 21, 14,
       37, 23, 47, 31, 24, 9, 64, 39, 30, 36, 46, 32, 50, 34, 32, 34,
       30, 28, 35, 57, 78, 88, 60, 78, 67, 82, 68, 97, 67, 65, 105,
       83, 101, 107, 88, 178, 110, 136, 118, 136, 153, 119, 89, 105,
       60, 148, 59, 73, 83, 49, 137, 91)

# Fitting the model
mod4 <- gamlss(y~1, sigma.fo=~1, family=DMOLBE,
               control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod4)

# Extracting the fitted values for mu and sigma
# using the inverse link function
mu_hat <- exp(coef(mod4, what="mu"))
mu_hat
sigma_hat <- exp(coef(mod4, what="sigma"))
sigma_hat
