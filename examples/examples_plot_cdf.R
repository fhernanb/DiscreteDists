# Example 1
# for a particular distribution

x <- 0:6
fx <- c(0, 0.19, 0.21, 0.4, 0.12, 0.05, 0.03)
plot_discrete_cdf(x, fx, las=1,
                  xlab="X",
                  ylab="Probability")

# Example 2
# for a Poisson distribution
x <- 0:10
fx <- dpois(x, lambda=3)
plot_discrete_cdf(x, fx, las=1,
                  xlab="X",
                  ylab="Probability",
                  main="CDF for Poisson")


