# Example 1
# Plotting the mass function for different parameter values

plot(x=0:5, y=dDBH(x=0:5, mu=0.1),
     type="h", lwd=2, col="dodgerblue", las=1,
     ylab="P(X=x)", xlab="X", ylim=c(0, 1),
     main="Probability mu=0.1")

plot(x=0:10, y=dDBH(x=0:10, mu=0.5),
     type="h", lwd=2, col="tomato", las=1,
     ylab="P(X=x)", xlab="X", ylim=c(0, 1),
     main="Probability mu=0.5")

plot(x=0:15, y=dDBH(x=0:15, mu=0.9),
     type="h", lwd=2, col="green4", las=1,
     ylab="P(X=x)", xlab="X", ylim=c(0, 1),
     main="Probability mu=0.9")

# Example 2
# Checking if the cumulative curves converge to 1

x_max <- 15
cumulative_probs1 <- pDBH(q=0:x_max, mu=0.1)
cumulative_probs2 <- pDBH(q=0:x_max, mu=0.5)
cumulative_probs3 <- pDBH(q=0:x_max, mu=0.9)

plot(x=0:x_max, y=cumulative_probs1, col="dodgerblue",
     type="o", las=1, ylim=c(0, 1),
     main="Cumulative probability for Burr-Hatke",
     xlab="X", ylab="Probability")
points(x=0:x_max, y=cumulative_probs2, type="o", col="tomato")
points(x=0:x_max, y=cumulative_probs3, type="o", col="green4")
legend("bottomright", col=c("dodgerblue", "tomato", "green4"), lwd=3,
       legend=c("mu=0.1",
                "mu=0.5",
                "mu=0.9"))

# Example 3
# Comparing the random generator output with
# the theoretical probabilities

mu <- 0.4
x_max <- 10
probs1 <- dDBH(x=0:x_max, mu=mu)
names(probs1) <- 0:x_max

x <- rDBH(n=1000, mu=mu)
probs2 <- prop.table(table(x))

cn <- union(names(probs1), names(probs2))
height <- rbind(probs1[cn], probs2[cn])
mp <- barplot(height, beside = TRUE, names.arg = cn,
              col=c("dodgerblue3","firebrick3"), las=1,
              xlab="X", ylab="Proportion")
legend("topright",
       legend=c("Theoretical", "Simulated"),
       bty="n", lwd=3,
       col=c("dodgerblue3","firebrick3"), lty=1)

# Example 4
# Checking the quantile function

mu <- 0.97
p <- seq(from=0, to=1, by = 0.01)
qxx <- qDBH(p, mu, lower.tail = TRUE, log.p = FALSE)
plot(p, qxx, type="s", lwd=2, col="green3", ylab="quantiles",
     main="Quantiles of BH(mu=0.97)")

