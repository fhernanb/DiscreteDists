# Example 1
# Plotting the mass function for diferent parameter values

x_max <- 30
probs1 <- dHYPERPO(x=0:x_max, mu=5, sigma=0.1)
probs2 <- dHYPERPO(x=0:x_max, mu=5, sigma=1.0)
probs3 <- dHYPERPO(x=0:x_max, mu=5, sigma=1.8)

# To plot the first k values
plot(x=0:x_max, y=probs1, type="o", lwd=2, col="dodgerblue", las=1,
     ylab="P(X=x)", xlab="X", main="Probability for hyper-Poisson",
     ylim=c(0, 0.20))
points(x=0:x_max, y=probs2, type="o", lwd=2, col="tomato")
points(x=0:x_max, y=probs3, type="o", lwd=2, col="green4")
legend("topright", col=c("dodgerblue", "tomato", "green4"), lwd=3,
       legend=c("mu=5, sigma=0.1",
                "mu=5, sigma=1.0",
                "mu=5, sigma=1.8"))

# Example 2
# Checking if the cumulative curves converge to 1

x_max <- 15
cumulative_probs1 <- pHYPERPO(q=0:x_max, mu=5, sigma=0.1)
cumulative_probs2 <- pHYPERPO(q=0:x_max, mu=5, sigma=1.0)
cumulative_probs3 <- pHYPERPO(q=0:x_max, mu=5, sigma=1.8)

plot(x=0:x_max, y=cumulative_probs1, col="dodgerblue",
     type="o", las=1, ylim=c(0, 1),
     main="Cumulative probability for hyper-Poisson",
     xlab="X", ylab="Probability")
points(x=0:x_max, y=cumulative_probs2, type="o", col="tomato")
points(x=0:x_max, y=cumulative_probs3, type="o", col="green4")
legend("bottomright", col=c("dodgerblue", "tomato", "green4"), lwd=3,
       legend=c("mu=5, sigma=0.1",
                "mu=5, sigma=1.0",
                "mu=5, sigma=1.8"))

# Example 3
# Comparing the random generator output with the theoretical probabilities

x_max <- 15
probs1 <- dHYPERPO(x=0:x_max, mu=3, sigma=1.1)
names(probs1) <- 0:x_max

x <- rHYPERPO(n=1000, mu=3, sigma=1.1)
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

mu <- 3
sigma <-3
p <- seq(from=0, to=1, by=0.01)
qxx <- qHYPERPO(p=p, mu=mu, sigma=sigma,
                lower.tail=TRUE, log.p=FALSE)
plot(p, qxx, type="s", lwd=2, col="green3", ylab="quantiles",
     main="Quantiles of HP(mu=3, sigma=3)")

