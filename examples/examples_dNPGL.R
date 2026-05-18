# Example 1
# Plotting the mass function for different parameter values

x_max <- 50
probs1 <- dNPGL(x=0:x_max, mu = 0.1, sigma = 2)
probs2 <- dNPGL(x=0:x_max, mu = 0.5, sigma = 5)
probs3 <- dNPGL(x=0:x_max, mu = 0.2, sigma = 6)
probs4 <- dNPGL(x=0:x_max, mu = 20, sigma = 2)

plot(x=0:x_max, y=probs1, type="h", lwd=2, col="dodgerblue", las=1,
     ylab="P(X=x)", xlab="X", main="Probability for dNPGL",
     ylim=c(0, 0.035))
legend("topright", legend="mu=0.1, sigma=2")

plot(x=0:x_max, y=probs2, type="h", lwd=2, col="tomato", las=1,
     ylab="P(X=x)", xlab="X", main="Probability for dNPGL",
     ylim=c(0, 0.1))
legend("topright", legend="mu=0.5, sigma=5")

plot(x=0:x_max, y=probs3, type="h", lwd=2, col="green4", las=1,
     ylab="P(X=x)", xlab="X", main="Probability for dNPGL",
     ylim=c(0, 0.03))
legend("topright", legend="mu=0.2, sigma=6")

plot(x=0:x_max, y=probs4, type="h", lwd=2, col="magenta", las=1,
     ylab="P(X=x)", xlab="X", main="Probability for dNPGL",
     ylim=c(0, 1))
legend("topright", legend="mu=20, sigma=2")

# Example 2
# Checking if the cumulative curves converge to 1

x_max <- 100
cumulative_probs1 <- pNPGL(q=0:x_max, mu = 0.1, sigma = 2)
cumulative_probs2 <- pNPGL(q=0:x_max, mu = 0.5, sigma = 5)
cumulative_probs3 <- pNPGL(q=0:x_max, mu = 0.2, sigma = 6)
cumulative_probs4 <- pNPGL(q=0:x_max, mu = 20, sigma = 2)

plot(x=0:x_max, y=cumulative_probs1, col="dodgerblue",
     type="o", las=1, ylim=c(0, 1),
     main="Cumulative probability for NPGL",
     xlab="X", ylab="Probability")
points(x=0:x_max, y=cumulative_probs2, type="o", col="tomato")
points(x=0:x_max, y=cumulative_probs3, type="o", col="green4")
points(x=0:x_max, y=cumulative_probs4, type="o", col="magenta")
legend("bottomright", col=c("dodgerblue", "tomato", "green4", "magenta"), lwd=3,
       legend=c("mu=0.1, sigma=2",
                "mu=0.5, sigma=5",
                "mu=0.2, sigma=6",
                "mu=20, sigma=2"))

# Example 3
# Comparing the random generator output with
# the theoretical probabilities

x_max <- 100
mu <- 0.1
sigma <- 2
probs1 <- dNPGL(x=0:x_max, mu=mu, sigma=sigma)
names(probs1) <- 0:x_max

x <- rNPGL(n=10000, mu=mu, sigma=sigma)
probs2 <- prop.table(table(x))

cn <- union(names(probs1), names(probs2))
height <- rbind(probs1[cn], probs2[cn])
nombres <- cn
mp <- barplot(height, beside = TRUE, names.arg = nombres,
              col=c("dodgerblue3","firebrick3"), las=1,
              xlab="X", ylab="Proportion")
legend("topright",
       legend=c("Theoretical", "Simulated"),
       bty="n", lwd=3,
       col=c("dodgerblue3","firebrick3"), lty=1)

# Example 4
# Checking the quantile function

mu <- 0.1
sigma <- 2
p <- seq(from=0, to=1, by=0.01)
qxx <- qNPGL(p=p, mu=mu, sigma=sigma, lower.tail=TRUE, log.p=FALSE)
plot(p, qxx, type="s", lwd=2, col="green3", ylab="quantiles",
     main="Quantiles of NPGL(mu=0.1, sigma=2)")

