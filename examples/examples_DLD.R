# Example 1
# Generating some random values with
# known mu
y <- rDLD(n=100, mu=0.3)

# Fitting the model
library(gamlss)
mod1 <- gamlss(y~1, family=DLD,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod1, what="mu"))

# Example 2
# Generating random values under some model

# A function to simulate a data set with Y ~ DLD
gendat <- function(n) {
  x1 <- runif(n)
  mu    <- exp(2 - 4 * x1)
  y <- rDLD(n=n, mu=mu)
  data.frame(y=y, x1=x1)
}

set.seed(1235)
datos <- gendat(n=150)

mod2 <- NULL
mod2 <- gamlss(y~x1, sigma.fo=~x2, family=DLD, data=datos,
                 control=gamlss.control(n.cyc=500, trace=FALSE))

summary(mod2)

# Example 3
# Survival times in days of 72 guinea pigs.
# Taken from Bakouch et al (2014) page 26.

y <- c(12, 15, 22, 24, 24, 32, 32, 33, 34, 38, 38, 43, 44, 48,
       52, 53, 54, 54, 55, 56, 57, 58, 58, 59, 60, 60, 60, 60,
       61, 62, 63, 65, 65, 67, 68, 70, 70, 72, 73, 75, 76, 76,
       81, 83, 84, 85, 87, 91, 95, 96, 98, 99, 109, 110, 121,
       127, 129, 131, 143, 146, 146, 175, 175, 211, 233, 258,
       258, 263, 297, 341, 341, 376)

mod3 <- gamlss(y~1, family=DLD,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod3, what="mu"))


# Example 4
# Remission times in weeks for 20 leukaemia patients.
# Taken from Bakouch et al (2014) page 26.

y <- c(1, 3, 3, 6, 7, 7, 10, 12, 14, 15, 18, 19,
       22, 26, 28, 29, 34, 40, 48, 49)

mod4 <- gamlss(y~1, family=DLD,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod4, what="mu"))

# Example 5
# Numbers of fires in Greece for the period from 1
# July 1998 to 31 August of the same year .
# Taken from Bakouch et al (2014) page 26.

y <- c(2, 4, 4, 3, 3, 1, 2, 4, 3, 1, 1, 0, 5, 5, 0, 3, 1,
       1, 0, 1, 0, 2, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,
       1, 4, 2, 2, 1, 2, 1, 2, 0, 2, 2, 1, 0, 3, 2, 1, 2,
       2, 7, 3, 5, 2, 5, 4, 5, 6, 5, 4, 3, 8, 43, 8, 4, 4,
       3, 10, 5, 4, 5, 12, 3, 8, 12, 10, 11, 6, 1, 8, 9,
       12, 9, 4, 8, 12, 11, 8, 6, 4, 7, 9, 15, 12, 15, 15,
       12, 9, 16, 7, 11, 9, 11, 6, 5, 20, 9, 8, 8, 5, 7, 10,
       6, 6, 5, 5, 15, 6, 8, 5, 6)

mod5 <- gamlss(y~1, family=DLD,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod5, what="mu"))

# Example 6
# Final examination marks of space students in
# mathematics in the Indian Institute of Technology at Kanpur.
# Taken from Bakouch et al (2014) page 26.

y <- c(29, 25, 50, 15, 13, 27, 15, 18, 7, 7, 8, 19, 12,
       18, 5, 21, 15, 86, 21, 15, 14, 39, 15, 14, 70,
       44, 6, 23, 58, 19, 50, 23, 11, 6, 34, 18, 28, 34,
       12, 37, 4, 60, 20, 23, 40, 65, 19, 31)

mod6 <- gamlss(y~1, family=DLD,
               control=gamlss.control(n.cyc=500, trace=FALSE))

# Extracting the fitted values for mu
# using the inverse link function
exp(coef(mod6, what="mu"))


