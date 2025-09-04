# Example 1
# Generating some random values with
# known mu and sigma
\dontrun{
set.seed(12)
y <- rCOMPO(n=100, mu=10, sigma=3)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, sigma.fo=~1, family=COMPO,
               control=gamlss.control(n.cyc=500, trace=TRUE))

# Extracting the fitted values for mu and sigma
# using the inverse link function
exp(coef(mod1, what="mu"))
exp(coef(mod1, what="sigma"))
}

# Example 2
# Generating random values under some model

\dontrun{
# A function to simulate a data set with Y ~ COMPO
gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(2 + 1 * x1) # 12 approximately
  sigma <- exp(2 - 2 * x2) # 2.71 approximately
  y <- rCOMPO(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

set.seed(123)
dat <- gendat(n=100)

# Fitting the model
mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=COMPO, data=dat,
               control=gamlss.control(n.cyc=500, trace=TRUE))

summary(mod2)
}

# Example 3
# Using the data from Shmueli et al. (2005) page 134
# The dataset consists of quarterly sales of a well-known brand of a
# particular article of clothing at stores of a large national retailer.
\dontrun{
values <- 0:30
freq <- c(514, 503, 457, 423, 326, 233, 195, 139, 101, 77, 56, 40,
          37, 22, 9, 7, 10, 9, 3, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1)

y <- rep(x=values, times=freq)

mod3 <- gamlss(y~1, sigma.fo=~1, family=COMPO,
               control=gamlss.control(n.cyc=500, trace=TRUE))

exp(coef(mod3, what="mu"))
exp(coef(mod3, what="sigma"))

estim_mu_sigma_COMPO(y)

library(COMPoissonReg)
fit <- glm.cmp(y ~ 1)
res <- exp(fit$opt.res$par)
res

}

# Example 4
# Using the data from Shmueli et al. (2005) page 134
# The dataset contains lengths of words (numbers of syllables)
# in a Hungarian dictionary

\dontrun{
# Slovak dictionary
y <- rep(x=1:5, times=c(7, 33, 49, 22, 6))

# Hungarian dictionary
y <- rep(x=1:9, times=c(1421, 12333, 20711, 15590, 5544, 1510, 289, 60, 1))

mod4 <- gamlss(y~1, sigma.fo=~1, family=COMPO,
               control=gamlss.control(n.cyc=500, trace=TRUE))

exp(coef(mod4, what="mu"))
exp(coef(mod4, what="sigma"))

estim_mu_sigma_COMPO(y)

library(COMPoissonReg)
fit <- glm.cmp(y ~ 1)
res <- exp(fit$opt.res$par)
res

}
