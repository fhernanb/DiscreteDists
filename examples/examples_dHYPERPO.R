# Example 1
x_max <- 30
probs1 <- dHYPERPO(x=0:x_max, sigma=0.01, mu=5)
probs2 <- dHYPERPO(x=0:x_max, sigma=0.50, mu=5)
probs3 <- dHYPERPO(x=0:x_max, sigma=1.00, mu=5)
# To plot the first k values
plot(x=0:x_max, y=probs1, type="o", lwd=2, col="dodgerblue", las=1,
     ylab="P(X=x)", xlab="X", main="Probabilities for hyper-Poisson",
     ylim=c(0, 0.20))
points(x=0:x_max, y=probs2, type="o", lwd=2, col="tomato")
points(x=0:x_max, y=probs3, type="o", lwd=2, col="green4")
legend("topright", col=c("dodgerblue", "tomato", "green4"), lwd=3,
       legend=c("sigma=0.01, mu=5",
                "sigma=0.50, mu=5",
                "sigma=1.00, mu=5"))

# Example 2
# Checking if curves go to 1
x_max <- 15
cumulative_probs1 <- pHYPERPO(q=0:x_max, mu=1, sigma=1.5)
cumulative_probs2 <- pHYPERPO(q=0:x_max, mu=3, sigma=1.5)
cumulative_probs3 <- pHYPERPO(q=0:x_max, mu=5, sigma=1.5)

plot(x=0:x_max, y=cumulative_probs1, col="dodgerblue",
     type="o", las=1, ylim=c(0, 1),
     main="Cumulative probility for hyper-Poisson",
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
probs1 <- dHYPERPO(x=0:x_max, mu=3, sigma=1.1)
names(probs1) <- 0:x_max

x <- rHYPERPO(n=10000, mu=3, sigma=1.1)
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

