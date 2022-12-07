# Example 1
x_max <- 30
probs1 <- dHYPERPO2(x=0:x_max, sigma=0.01, mu=3)
probs2 <- dHYPERPO2(x=0:x_max, sigma=0.50, mu=5)
probs3 <- dHYPERPO2(x=0:x_max, sigma=1.00, mu=7)
# To plot the first k values
plot(x=0:x_max, y=probs1, type="o", lwd=2, col="dodgerblue", las=1,
     ylab="P(X=x)", xlab="X", main="Probability for hyper-Poisson",
     ylim=c(0, 0.30))
points(x=0:x_max, y=probs2, type="o", lwd=2, col="tomato")
points(x=0:x_max, y=probs3, type="o", lwd=2, col="green4")
legend("topright", col=c("dodgerblue", "tomato", "green4"), lwd=3,
       legend=c("sigma=0.01, mu=3",
                "sigma=0.50, mu=5",
                "sigma=1.00, mu=7"))

# Example 2
# Checking if curves go to 1
x_max <- 15
cumulative_probs1 <- pHYPERPO2(q=0:x_max, mu=1, sigma=1.5)
cumulative_probs2 <- pHYPERPO2(q=0:x_max, mu=3, sigma=1.5)
cumulative_probs3 <- pHYPERPO2(q=0:x_max, mu=5, sigma=1.5)

plot(x=0:x_max, y=cumulative_probs1, col="dodgerblue",
     type="o", las=1, ylim=c(0, 1),
     main="Cumulative probability for hyper-Poisson",
     xlab="X", ylab="Probability")
points(x=0:x_max, y=cumulative_probs2, type="o", col="tomato")
points(x=0:x_max, y=cumulative_probs3, type="o", col="green4")
legend("bottomright", col=c("dodgerblue", "tomato", "green4"), lwd=3,
       legend=c("sigma=1.5, mu=1",
                "sigma=1.5, mu=3",
                "sigma=1.5, mu=5"))

# Example 3
# Some tests to compare the simulator with theoretical probabilities
x_max <- 15
probs1 <- dHYPERPO2(x=0:x_max, mu=3, sigma=1.1)
names(probs1) <- 0:x_max

x <- rHYPERPO2(n=1000, mu=3, sigma=1.1)
probs2 <- prop.table(table(x))

cn <- union(names(probs1), names(probs2))
height <- rbind(probs1[cn], probs2[cn])
nombres <- cn
mp <- barplot(height, beside = TRUE, names.arg = nombres,
              col=c('dodgerblue3','firebrick3'), las=1,
              xlab='X', ylab='Proportion')
legend('topright',
       legend=c('Theoretical', 'Simulated'),
       bty='n', lwd=3,
       col=c('dodgerblue3','firebrick3'), lty=1)

# Example 4
p <- seq(from=0, to=0.99999, length.out=100)
plot(x=qHYPERPO2(p, mu=5, sigma=1.5), y=p, xlab="Quantile",
     las=1, ylab="Probability")
curve(pHYPERPO2(x, mu=5, sigma=1.5), from=0, add=TRUE, col="red")

# Example 5
p <- seq(0,1, by = 0.01)
mu <- 3
sigma<-3
qxx <- qHYPERPO2(p, mu, sigma, lower.tail = TRUE, log.p = FALSE)
plot(p, qxx, type="s", lwd=2, col="darkred", ylab="quantiles",
     main="Quantiles of HP(mu = sigma = 3)")
