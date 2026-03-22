# Example 1
# Plotting the mass function for different parameter values

x_max <- 20
probs1 <- dPTRTE(x=0:x_max, mu=0.5, sigma=0.5)
probs2 <- dPTRTE(x=0:x_max, mu=1, sigma=0.5)
probs3 <- dPTRTE(x=0:x_max, mu=2, sigma=0.5)

# To plot the first k values
plot(x=0:x_max, y=probs1, type="o", lwd=2, col="dodgerblue", las=1,
     ylab="P(X=x)", xlab="X", main="Probability for PTRTE",
     ylim=c(0, 0.80))
points(x=0:x_max, y=probs2, type="o", lwd=2, col="tomato")
points(x=0:x_max, y=probs3, type="o", lwd=2, col="green4")
legend("topright", col=c("dodgerblue", "tomato", "green4"), lwd=3,
       legend=c("mu=0.5, sigma=0.5",
                "mu=1, sigma=0.5",
                "mu=2, sigma=0.5"))

# Example 2
# Checking if the cumulative curves converge to 1

x_max <- 10

cumulative_probs1 <- pPTRTE(q=0:x_max, mu=1, sigma=0.5)
cumulative_probs2 <- pPTRTE(q=0:x_max, mu=1, sigma=0.7)
cumulative_probs3 <- pPTRTE(q=0:x_max, mu=1, sigma=0.9)

plot(x=0:x_max, y=cumulative_probs1, col="dodgerblue",
     type="o", las=1, ylim=c(0, 1),
     main="Cumulative for PTRTE",
     xlab="X", ylab="Probability")
points(x=0:x_max, y=cumulative_probs2, type="o", col="tomato")
points(x=0:x_max, y=cumulative_probs3, type="o", col="green4")
legend("bottomright", col=c("dodgerblue", "tomato", "green4"), lwd=3,
       legend=c("mu=1, sigma=0.5",
                "mu=1, sigma=0.7",
                "mu=1, sigma=0.9"))

# Example 3
# Comparing the random generator output with
# the theoretical probabilities

x_max <- 30
probs1 <- dPTRTE(x=0:x_max, mu=0.5, sigma=0.5)
names(probs1) <- 0:x_max

x <- rPTRTE(n=1000, mu=0.5, sigma=0.5)
probs2 <- prop.table(table(x))

cn <- union(names(probs1), names(probs2))
height <- rbind(probs1[cn], probs2[cn])
mp <- barplot(height, beside=TRUE, names.arg=cn,
              col=c("dodgerblue3","firebrick3"), las=1,
              xlab="X", ylab="Proportion")
legend("topright",
       legend=c("Theoretical", "Simulated"),
       bty="n", lwd=3,
       col=c("dodgerblue3","firebrick3"), lty=1)

# Example 4
# Checking the quantile function

mu <- 1
sigma <- 0.5
p <- seq(from=0, to=1, by=0.01)
qxx <- qPTRTE(p=p, mu=mu, sigma=sigma, lower.tail=TRUE, log.p=FALSE)
plot(p, qxx, type="s", lwd=2, col="green3", ylab="quantiles",
     main="Quantiles of DPTRTE(mu=1, sigma=0.5)")
